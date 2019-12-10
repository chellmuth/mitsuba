#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>

#include <mitsuba/core/fstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/render/gatherproc.h>
#include <mitsuba/render/photon.h>
#include <mitsuba/render/photonmap.h>

#include "neural_frame.h"
#include "photon_bundle.h"
#include "photon_helper.h"

MTS_NAMESPACE_BEGIN

class ImportanceTestingIntegrator : public SamplingIntegrator {
public:
    ImportanceTestingIntegrator(const Properties &props) : SamplingIntegrator(props)
    {
        m_rrDepth = props.getInteger("rrDepth", 5);
        m_maxDepth = props.getInteger("maxDepth", -1);
        m_strictNormals = props.getBoolean("strictNormals", false);
        m_hideEmitters = props.getBoolean("hideEmitters", false);

        m_x = props.getInteger("x");
        m_y = props.getInteger("y");
        m_pdfCount = props.getSize("pdfCount");

        Properties gtProps("path");
        gtProps.setInteger("maxDepth", m_maxDepth);
        gtProps.setBoolean("strictNormals", m_strictNormals);
        gtProps.setBoolean("hideEmitters", m_hideEmitters);

        Properties fisheyeProps("path");
        fisheyeProps.setInteger("maxDepth", m_maxDepth - 1);
        fisheyeProps.setBoolean("strictNormals", m_strictNormals);
        fisheyeProps.setBoolean("hideEmitters", m_hideEmitters);

        std::cout << gtProps.toString() << std::endl;
        std::cout << fisheyeProps.toString() << std::endl;

        m_gtIntegrator = static_cast<SamplingIntegrator *>(
            PluginManager::getInstance()->createObject(MTS_CLASS(SamplingIntegrator), gtProps)
        );

        m_fisheyeIntegrator = static_cast<SamplingIntegrator *>(
            PluginManager::getInstance()->createObject(MTS_CLASS(SamplingIntegrator), fisheyeProps)
        );

        Log(EInfo, "ImportanceTesting constructor (%i, %i, %i)", m_x, m_y, m_pdfCount);
    }

    ImportanceTestingIntegrator(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager) {
        // m_integrator = static_cast<SamplingIntegrator *>(manager->getInstance(stream));
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);

        // manager->serialize(stream, m_integrator.get());
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

        // if (!m_integrator)
        //     Log(EError, "No integrator supplied to the fisheye integrator!");

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
        const int spp = 1024;

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
        testNeuralFrame(neuralFrame, rRec.its);

        const BSDF *bsdf = rRec.its.getBSDF(sensorRay);

        Spectrum L(0.f);

        RadianceQueryRecord nestedRec(scene, sampler);
        for (int thetaStep = 0; thetaStep < thetaSteps; thetaStep++) {
            for (int phiStep = 0; phiStep < phiSteps; phiStep++) {
                const Float theta = (M_PI / 2.f) * (thetaStep + 0.5f) / thetaSteps;
                const Float phi = (M_PI * 2.f) * (phiStep + 0.5f) / phiSteps;

                Spectrum result(0.f);
                for (int sample = 0; sample < spp; sample++) {
                    nestedRec.newQuery(RadianceQueryRecord::ESensorRay, sensor->getMedium());

                    Float sinTheta, cosTheta;
                    Float sinPhi, cosPhi;
                    math::sincos(theta, &sinTheta, &cosTheta);
                    math::sincos(phi, &sinPhi, &cosPhi);

                    const Vector woLocal(cosPhi * sinTheta, cosTheta, sinPhi * sinTheta);
                    const Vector wo = neuralFrame.toWorld(woLocal);

                    RayDifferential fisheyeRay(rRec.its.p, wo, sensorRay.time);
                    fisheyeRay.mint = Epsilon;

                    BSDFSamplingRecord bRec(rRec.its, rRec.sampler, ERadiance);
                    result += m_fisheyeIntegrator->Li(fisheyeRay, nestedRec)
                        // * bsdf->eval(bRec, ESolidAngle)
                        * sinf(theta)
                        * (0.5f / M_PI)
                        * (2.f * M_PI) * (M_PI / 2.f) / (thetaSteps * phiSteps)
                        * cosTheta
                        * (1.f / spp);
                }

                bitmap->setPixel({phiStep, thetaStep}, result);

                L += result;
            }
        }

        std::cout << "RADIANCE TEST RESULT!: " << L.toString() << std::endl;

        film->setBitmap(bitmap);
        film->develop(scene, 0.f);
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
                RadianceQueryRecord rRec2(rRec);

                Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

                Spectrum spec = sensor->sampleRayDifferential(
                    sensorRay, samplePos, apertureSample, timeSample);

                sensorRay.scaleDifferential(diffScaleFactor);

                rRec.rayIntersect(sensorRay);

                if (j == 0 && offset.x == m_x && offset.y == m_y) {
                    renderFisheye(scene, sensor, sampler, rRec, sensorRay);
                }

                Spectrum result = spec * m_gtIntegrator->Li(sensorRay, rRec2);
                // std::cout << spec.toString() << " " << result.toString() << std::endl;
                block->put(samplePos, result, rRec2.alpha);
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

        // m_integrator->bindUsedResources(proc);
    }

    void wakeup(ConfigurableObject *parent, std::map<std::string, SerializableObject *> &params) {
        SamplingIntegrator::wakeup(parent, params);

        // m_integrator->wakeup(parent, params);
    }

    // void addChild(const std::string &name, ConfigurableObject *child) {
    //     if (child->getClass()->derivesFrom(MTS_CLASS(SamplingIntegrator))) {
    //         SamplingIntegrator *integrator = static_cast<SamplingIntegrator *>(child);
    //         m_integrator = integrator;
    //         integrator->incRef();
    //     } else {
    //         SamplingIntegrator::addChild(name, child);
    //     }
    // }

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        NotImplementedError("Li");
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "ImportanceTestingIntegrator" << endl;
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    int m_maxDepth;
    int m_rrDepth;
    bool m_strictNormals;
    bool m_hideEmitters;

    ref<SamplingIntegrator> m_gtIntegrator;
    ref<SamplingIntegrator> m_fisheyeIntegrator;

    // ref<SamplingIntegrator> m_integrator;
    int m_x;
    int m_y;
    int m_pdfCount;

    ref<ParallelProcess> m_proc;
};

MTS_IMPLEMENT_CLASS_S(ImportanceTestingIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(ImportanceTestingIntegrator, "Multi-channel integrator");
MTS_NAMESPACE_END
