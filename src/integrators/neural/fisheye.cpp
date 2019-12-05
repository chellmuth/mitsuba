#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>

#include <mitsuba/core/fstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/render/gatherproc.h>
#include <mitsuba/render/photon.h>
#include <mitsuba/render/photonmap.h>

#include "photon_bundle.h"
#include "neural_frame.h"

MTS_NAMESPACE_BEGIN

class FisheyeIntegrator : public SamplingIntegrator {
public:
    typedef PointKDTree<Photon>      PhotonTree;
    typedef PhotonTree::SearchResult SearchResult;

    FisheyeIntegrator(const Properties &props) : SamplingIntegrator(props)
    {
        m_rrDepth = props.getInteger("rrDepth", 5);
        m_maxDepth = props.getInteger("maxDepth", -1);

        m_x = props.getInteger("x");
        m_y = props.getInteger("y");
        m_pdfCount = props.getSize("pdfCount");

        m_globalPhotons = props.getSize("globalPhotons", 10000);

        Log(EInfo, "Fisheye constructor (%i, %i, %i)", m_x, m_y, m_pdfCount);
    }

    FisheyeIntegrator(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager) {
        m_integrator = static_cast<SamplingIntegrator *>(manager->getInstance(stream));
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);

        manager->serialize(stream, m_integrator.get());
    }

    bool preprocess(
        const Scene *scene,
        RenderQueue *queue,
        const RenderJob *job,
        int sceneResID,
        int sensorResID,
        int samplerResID
    ) {
        if (!SamplingIntegrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID)) {
            Log(EError, "Base class error");
            return false;
        }

        if (!m_integrator->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID)) {
            Log(EError, "My integrator error");
            return false;
        }

        if (m_globalPhotonMap.get() == NULL && m_globalPhotons > 0) {
            /* Generate the global photon map */
            ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
                GatherPhotonProcess::ESurfacePhotons,
                m_globalPhotons, // count
                0, // granularity for parallelization
                m_maxDepth - 1,
                m_rrDepth,
                true, // gather locally
                true, // auto-cancel if not enough photons are generated
                job
            );

            proc->bindResource("scene", sceneResID);
            proc->bindResource("sensor", sensorResID);
            proc->bindResource("sampler", samplerResID);

            ref<Scheduler> sched = Scheduler::getInstance();
            m_proc = proc;
            sched->schedule(proc);
            sched->wait(proc);
            m_proc = NULL;

            if (proc->getReturnStatus() != ParallelProcess::ESuccess) {
                return false;
            }

            ref<PhotonMap> globalPhotonMap = proc->getPhotonMap();
            if (globalPhotonMap->isFull()) {
                Log(EDebug, "Global photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: "
                    SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

                // m_globalPhotonMap = globalPhotonMap;
                // m_globalPhotonMap->setScaleFactor(1 / (Float) proc->getShotParticles());
                // m_globalPhotonMap->build();
                // m_globalPhotonMapID = sched->registerResource(m_globalPhotonMap);
            }

            m_globalPhotonMap = globalPhotonMap;
            m_globalPhotonMap->setScaleFactor(1 / (Float) proc->getShotParticles());
            m_globalPhotonMap->build();
        }

        return true;
    }

    bool render(
        Scene *scene,
        RenderQueue *queue, const RenderJob *job,
        int sceneResID, int sensorResID, int samplerResID
    ) {
        ref<Scheduler> sched = Scheduler::getInstance();
        ref<Sensor> sensor = static_cast<Sensor *>(sched->getResource(sensorResID));
        ref<Film> film = sensor->getFilm();

        size_t nCores = sched->getCoreCount();
        const Sampler *sampler = static_cast<const Sampler *>(sched->getResource(samplerResID, 0));
        size_t sampleCount = sampler->getSampleCount();

        if (!m_integrator)
            Log(EError, "No integrator supplied to the fisheye integrator!");

        Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SIZE_T_FMT
            " %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y,
            sampleCount, sampleCount == 1 ? "sample" : "samples", nCores,
            nCores == 1 ? "core" : "cores");

        /* This is a sampling-based integrator - parallelize */
        ref<BlockedRenderProcess> proc = new BlockedRenderProcess(job, queue, scene->getBlockSize());

        // proc->setPixelFormat(
        //     m_integrators.size() > 1 ? Bitmap::EMultiSpectrumAlphaWeight : Bitmap::ESpectrumAlphaWeight,
        //     (int) (m_integrators.size() * SPECTRUM_SAMPLES + 2), false
        // );

        int integratorResID = sched->registerResource(this);
        proc->bindResource("integrator", integratorResID);
        proc->bindResource("scene", sceneResID);
        proc->bindResource("sensor", sensorResID);
        proc->bindResource("sampler", samplerResID);
        scene->bindUsedResources(proc);
        bindUsedResources(proc);
        sched->schedule(proc);

        m_process = proc;
        sched->wait(proc);
        m_process = NULL;
        sched->unregisterResource(integratorResID);

        return proc->getReturnStatus() == ParallelProcess::ESuccess;
    }

    void gatherPhotons(
        const RadianceQueryRecord &rRec,
        bool flipNormal,
        const ImageBlock *block = nullptr,
        const int identifier = -1
    ) const {
        const Intersection &its = rRec.its;

        const size_t maxPhotons = 100;
        SearchResult *results = static_cast<SearchResult *>(
            alloca((maxPhotons + 1) * sizeof(SearchResult)));

        size_t resultCount = m_globalPhotonMap->nnSearch(its.p, maxPhotons, results);
        // Log(EInfo, "Photons returned: %i", resultCount);

        // std::cout << "INTERSECTION RECORD:" << std::endl;
        // std::cout << its.p.toString() << std::endl;
        // std::cout << its.geoFrame.n.toString() << std::endl;
        // std::cout << its.wi.toString() << std::endl;

        std::ostringstream oss;
        if (identifier < 0) {
            oss << "photons_" << m_x << "_" << m_y << ".bin";
        } else {
            auto offset = Vector2i(block->getOffset());
            oss << "results/photons_" << identifier << "-block_" << offset.x << "x" << offset.y << ".bin";
        }

        ref<FileStream> fileStream = new FileStream(oss.str(), FileStream::ETruncWrite);

        float intersectionBuffer[] = {
            its.p.x,
            its.p.y,
            its.p.z,
            its.shFrame.n.x * (flipNormal ? -1.f : 1.f),
            its.shFrame.n.y * (flipNormal ? -1.f : 1.f),
            its.shFrame.n.z * (flipNormal ? -1.f : 1.f),
            its.wi.x,
            its.wi.y,
            its.wi.z,
        };
        fileStream->write(&intersectionBuffer, 9 * sizeof(float));

        uint32_t countBuffer = resultCount;
        fileStream->write(&countBuffer, sizeof(uint32_t));

        for (size_t i = 0; i < resultCount; i++) {
            const SearchResult &searchResult = results[i];
            const PhotonMap &photonMap = (*m_globalPhotonMap.get());
            const Photon &photon = photonMap[searchResult.index];

            Point position = photon.getPosition();
            Point source = photon.getSource();

            float r, g, b;
            Spectrum power = photon.getPower();
            power.toLinearRGB(r, g, b);

            float photonBuffer[] = {
                position.x,
                position.y,
                position.z,

                source.x,
                source.y,
                source.z,

                r, g, b
            };

            fileStream->write(photonBuffer, sizeof(float) * 9);

        //     std::cout << "PHOTON RECORD:" << std::endl;
        //     std::cout << photon.getPosition().toString() << std::endl;
        //     std::cout << photon.getSource().toString() << std::endl;
        //     std::cout << photon.getPower().toString() << std::endl;
        }

        fileStream->close();

        // {
        //     std::cout << "INTERSECTION:" << std::endl;
        //     std::cout << its.toString() << std::endl;

        //     Vector normal = its.shFrame.n;
        //     if (flipNormal) {
        //         normal *= -1.f;
        //     }
        //     Frame neuralFrame = constructNeuralFrame(normal, its);
        //     PhotonBundle bundle(its.p, neuralFrame, 10, 10);

        //     for (size_t i = 0; i < resultCount; i++) {
        //         const SearchResult &searchResult = results[i];
        //         const PhotonMap &photonMap = (*m_globalPhotonMap.get());
        //         const Photon &photon = photonMap[searchResult.index];

        //         bundle.splat(photon);
        //     }

        //     std::vector<Float> photonBundle = bundle.serialized();

        //     std::ostringstream oss;
        //     oss << "grid_" << m_x << "_" << m_y << ".bin";

        //     ref<FileStream> fileStream = new FileStream(oss.str(), FileStream::ETruncWrite);
        //     fileStream->write(photonBundle.data(), sizeof(float) * 100);
        //     fileStream->close();
        // }
    }

    void renderFisheye(
        const Scene *scene,
        const Sensor *sensor,
        Sampler *sampler,
        RadianceQueryRecord &rRec,
        const RayDifferential &sensorRay,
        const ImageBlock *block = nullptr,
        const int identifier = -1
    ) const {
        const int thetaSteps = 400;
        const int phiSteps = 400;
        const int spp = 32;

        Properties filmProps("HDRFilm");
        filmProps.setInteger("width", phiSteps);
        filmProps.setInteger("height", thetaSteps);
        filmProps.setBoolean("banner", false);

        ref<Film> film = static_cast<Film *>(
            PluginManager::getInstance()->createObject(MTS_CLASS(Film), filmProps)
        );

        std::ostringstream oss;
        if (identifier < 0) {
            oss << "render_" << m_x << "_" << m_y << ".exr";
        } else {
            auto offset = Vector2i(block->getOffset());
            oss << "results/pdf_" << identifier << "-block_" << offset.x << "x" << offset.y << ".exr";
        }
        film->setDestinationFile(oss.str(), 0);

        ref<Bitmap> bitmap = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, {phiSteps, thetaSteps});

        bool flippedNormal = false;
        Vector normal = rRec.its.shFrame.n;
        if (Frame::cosTheta(rRec.its.wi) < 0.f) {
            flippedNormal = true;
            normal *= -1.f;
        }

        Frame neuralFrame = constructNeuralFrame(normal, rRec.its);

        RadianceQueryRecord nestedRec(scene, sampler);
        for (int thetaStep = 0; thetaStep < thetaSteps; thetaStep++) {
            for (int phiStep = 0; phiStep < phiSteps; phiStep++) {
                Spectrum result(0.f);
                for (int sample = 0; sample < spp; sample++) {
                    const Float theta = (M_PI / 2.f) * (thetaStep + 0.5f) / thetaSteps;
                    const Float phi = (M_PI * 2.f) * (phiStep + 0.5f) / phiSteps;

                    nestedRec.newQuery(RadianceQueryRecord::ESensorRay, sensor->getMedium());

                    Float sinTheta, cosTheta;
                    Float sinPhi, cosPhi;
                    math::sincos(theta, &sinTheta, &cosTheta);
                    math::sincos(phi, &sinPhi, &cosPhi);

                    const Vector woLocal(cosPhi * sinTheta, cosTheta, sinPhi * sinTheta);
                    const Vector wo = neuralFrame.toWorld(woLocal);

                    // std::cout << "theta: " << thetaStep << " " << "phi: " << phiStep << std::endl;
                    // std::cout << woLocal.toString() << std::endl;
                    // std::cout << wo.toString() << std::endl;
                    // std::cout << rRec.its.p.toString() << std::endl;
                    // std::cout << rRec.its.shFrame.n.toString() << std::endl;
                    // std::cout << std::endl;

                    RayDifferential fisheyeRay(rRec.its.p, wo, sensorRay.time);
                    fisheyeRay.mint = Epsilon;

                    result += m_integrator->Li(fisheyeRay, nestedRec) * (1.f / spp);
                }

                bitmap->setPixel({phiStep, thetaStep}, result);
            }
        }

        film->setBitmap(bitmap);
        film->develop(scene, 0.f);

        gatherPhotons(rRec, flippedNormal, block, identifier);
    }

    void renderBlock(
        const Scene *scene,
        const Sensor *sensor, Sampler *sampler, ImageBlock *block,
        const bool &stop, const std::vector< TPoint2<uint8_t> > &points
    ) const {
        Float diffScaleFactor = 1.0f / std::sqrt((Float) sampler->getSampleCount());

        RadianceQueryRecord rRec(scene, sampler);
        Point2 apertureSample(0.5f);
        Float timeSample = 0.5f;
        RayDifferential sensorRay;

        block->clear();

        uint32_t queryType = RadianceQueryRecord::ESensorRay;

        for (size_t i = 0; i < points.size(); ++i) {
            Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
            if (stop)
                break;

            sampler->generate(offset);

            for (size_t j = 0; j < sampler->getSampleCount(); j++) {
                if (offset.x == m_x && offset.y == m_y) {}
                else { continue; }

                rRec.newQuery(queryType, sensor->getMedium());
                Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

                Spectrum spec = sensor->sampleRayDifferential(
                    sensorRay, samplePos, apertureSample, timeSample);

                sensorRay.scaleDifferential(diffScaleFactor);
                rRec.rayIntersect(sensorRay);

                if (offset.x == m_x && offset.y == m_y) {
                    renderFisheye(scene, sensor, sampler, rRec, sensorRay);
                }

                RadianceQueryRecord rRec2(rRec);
                rRec2.its = rRec.its;

                Spectrum result = spec * m_integrator->Li(sensorRay, rRec2);
                block->put(samplePos, result, rRec.alpha);
                sampler->advance();
            }
        }

        for (int i = 0; i < m_pdfCount; i++) {
            rRec.newQuery(queryType, sensor->getMedium());
            float sampleX = rRec.nextSample1D() * block->getWidth() + block->getOffset().x;
            float sampleY = rRec.nextSample1D() * block->getHeight() + block->getOffset().y;

            Point2 samplePos(sampleX, sampleY);
            sensor->sampleRayDifferential(sensorRay, samplePos, apertureSample, timeSample);
            sensorRay.scaleDifferential(diffScaleFactor);

            rRec.rayIntersect(sensorRay);

            if (rRec.its.isValid()) {
                renderFisheye(scene, sensor, sampler, rRec, sensorRay, block, i);
            }
        }

    }

    void bindUsedResources(ParallelProcess *proc) const {
        SamplingIntegrator::bindUsedResources(proc);

        m_integrator->bindUsedResources(proc);
    }

    void wakeup(ConfigurableObject *parent, std::map<std::string, SerializableObject *> &params) {
        SamplingIntegrator::wakeup(parent, params);

        m_integrator->wakeup(parent, params);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(SamplingIntegrator))) {
            SamplingIntegrator *integrator = static_cast<SamplingIntegrator *>(child);
            m_integrator = integrator;
            integrator->incRef();
        } else {
            SamplingIntegrator::addChild(name, child);
        }
    }

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        NotImplementedError("Li");
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "FisheyeIntegrator" << endl;
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    int m_maxDepth;
    int m_rrDepth;

    ref<SamplingIntegrator> m_integrator;
    int m_x;
    int m_y;
    int m_pdfCount;

    ref<ParallelProcess> m_proc;
    ref<PhotonMap> m_globalPhotonMap;
    size_t m_globalPhotons;
};

MTS_IMPLEMENT_CLASS_S(FisheyeIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(FisheyeIntegrator, "Multi-channel integrator");
MTS_NAMESPACE_END
