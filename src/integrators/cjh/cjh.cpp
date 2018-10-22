#include <mitsuba/bidir/util.h>
#include <mitsuba/core/plugin.h>

#include "cjh.h"
#include "cjh_sampler.h"

#include <mitsuba/bidir/pathsampler.h>

MTS_NAMESPACE_BEGIN

class CJH : public Integrator {
public:
    CJH(const Properties &props) : Integrator(props) {}

    CJH(Stream *stream, InstanceManager *manager)
     : Integrator(stream, manager) {
        m_config = CJHConfiguration(stream);
        configure();
    }

    virtual ~CJH() { }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Integrator::serialize(stream, manager);
        m_config.serialize(stream);
    }

    bool preprocess(const Scene *scene, RenderQueue *queue,
            const RenderJob *job, int sceneResID, int sensorResID,
            int samplerResID) {
        Integrator::preprocess(scene, queue, job, sceneResID,
                sensorResID, samplerResID);

        return true;
    }

    void cancel() {}

    bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {

        ref<Sampler> sampler = new CJHSampler();

        ref<PathSampler> pathSampler = new PathSampler(
            PathSampler::EUnidirectional,
            scene,
            sampler, sampler, sampler, // {emitter, sensor, direct} samplers
            3, // maxDepth
            9999, // russian roulette start depth
            false, // exclude direct illumination
            true // custom sample direct logic
        );

        SplatList *current = new SplatList();
        pathSampler->sampleSplats(Point2i(-1), *current);

        std::cout << current->toString() << std::endl;

        printf("RENDERED EVERYTHING!!\n");
        return true;
    }

    MTS_DECLARE_CLASS()
private:
    CJHConfiguration m_config;
};

MTS_IMPLEMENT_CLASS_S(CJH, false, Integrator)
MTS_EXPORT_PLUGIN(CJH, "CJH Testing Plugin");
MTS_NAMESPACE_END
