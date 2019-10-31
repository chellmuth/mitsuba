#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>

MTS_NAMESPACE_BEGIN

class FisheyeIntegrator : public SamplingIntegrator {
public:
    FisheyeIntegrator(const Properties &props) : SamplingIntegrator(props) { }

    FisheyeIntegrator(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager) {
        m_integrator = static_cast<SamplingIntegrator *>(manager->getInstance(stream));
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);

        manager->serialize(stream, m_integrator.get());
    }

    bool preprocess(
        const Scene *scene, RenderQueue *queue,
        const RenderJob *job, int sceneResID, int sensorResID,
        int samplerResID
    ) {
        if (!SamplingIntegrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID)) {
            return false;
        }

        if (m_integrator->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID)) {
            return false;
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

    void renderBlock(
        const Scene *scene,
        const Sensor *sensor, Sampler *sampler, ImageBlock *block,
        const bool &stop, const std::vector< TPoint2<uint8_t> > &points
    ) const {

        // Float diffScaleFactor = 1.0f /
        //     std::sqrt((Float) sampler->getSampleCount());

        // bool needsApertureSample = sensor->needsApertureSample();
        // bool needsTimeSample = sensor->needsTimeSample();

        // RadianceQueryRecord rRec(scene, sampler);
        // Point2 apertureSample(0.5f);
        // Float timeSample = 0.5f;
        // RayDifferential sensorRay;

        // block->clear();

        // uint32_t queryType = RadianceQueryRecord::ESensorRay;
        // Float *temp = (Float *) alloca(sizeof(Float) * (m_integrators.size() * SPECTRUM_SAMPLES + 2));

        // for (size_t i = 0; i<points.size(); ++i) {
        //     Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
        //     if (stop)
        //         break;

        //     sampler->generate(offset);

        //     for (size_t j = 0; j<sampler->getSampleCount(); j++) {
        //         rRec.newQuery(queryType, sensor->getMedium());
        //         Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

        //         if (needsApertureSample)
        //             apertureSample = rRec.nextSample2D();

        //         if (needsTimeSample)
        //             timeSample = rRec.nextSample1D();

        //         Spectrum spec = sensor->sampleRayDifferential(
        //             sensorRay, samplePos, apertureSample, timeSample);

        //         sensorRay.scaleDifferential(diffScaleFactor);
        //         rRec.rayIntersect(sensorRay);

        //         int offset = 0;
        //         for (size_t k = 0; k<m_integrators.size(); ++k) {
        //             RadianceQueryRecord rRec2(rRec);
        //             rRec2.its = rRec.its;
        //             Spectrum result = spec * m_integrators[k]->Li(sensorRay, rRec2);
        //             for (int l = 0; l<SPECTRUM_SAMPLES; ++l)
        //                 temp[offset++] = result[l];
        //         }
        //         temp[offset++] = rRec.alpha;
        //         temp[offset] = 1.0f;
        //         block->put(samplePos, temp);
        //         sampler->advance();
        //     }
        // }
    }

    void bindUsedResources(ParallelProcess *proc) const {
        SamplingIntegrator::bindUsedResources(proc);

        m_integrator->bindUsedResources(proc);
    }

    void wakeup(ConfigurableObject *parent, std::map<std::string, SerializableObject *> &params) {
        SamplingIntegrator::wakeup(parent, params);

        m_integrator->wakeup(parent, params);
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
    ref<SamplingIntegrator> m_integrator;
};

MTS_IMPLEMENT_CLASS_S(FisheyeIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(FisheyeIntegrator, "Multi-channel integrator");
MTS_NAMESPACE_END
