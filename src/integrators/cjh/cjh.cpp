#include <mitsuba/bidir/util.h>
#include <mitsuba/core/plugin.h>

#include "cjh.h"
#include "cjh_sampler.h"

#include <mitsuba/bidir/pathsampler.h>

MTS_NAMESPACE_BEGIN

class CJH : public Integrator {
public:
    CJH(const Properties &props) : Integrator(props) {
        m_config.x = props.getFloat("x");
        m_config.y = props.getFloat("y");

        m_config.u = props.getFloat("u");
        m_config.v = props.getFloat("v");

        m_config.direct1_1 = props.getFloat("direct1_1");
        m_config.direct1_2 = props.getFloat("direct1_2");

        m_config.bsdf1_1 = props.getFloat("bsdf1_1");
        m_config.bsdf1_2 = props.getFloat("bsdf1_2");

        m_config.direct2_1 = props.getFloat("direct2_1");
        m_config.direct2_2 = props.getFloat("direct2_2");

        m_config.bsdf2_1 = props.getFloat("bsdf2_1");
        m_config.bsdf2_2 = props.getFloat("bsdf2_2");
    }

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

        ref<CJHSampler> sensorSampler = new CJHSampler("sensor");
        ref<CJHSampler> emitterSampler = new CJHSampler("emitter");
        ref<CJHSampler> directSampler = new CJHSampler("direct");

        std::vector<Float> sensorSamples = {
            m_config.u,
            m_config.v,

            m_config.x,
            m_config.y,

            // direct light
            m_config.direct1_1,
            m_config.direct1_2,

            // bsdf
            m_config.bsdf1_1,
            m_config.bsdf1_2,

            // direct light
            m_config.direct2_1,
            m_config.direct2_2,

            // bsdf
            m_config.bsdf2_1,
            m_config.bsdf2_2
        };

        std::vector<Float> emitterSamples = {};
        std::vector<Float> directSamples = {};

        sensorSampler->setSamples(sensorSamples);
        emitterSampler->setSamples(emitterSamples);
        directSampler->setSamples(directSamples);

        ref<PathSampler> pathSampler = new PathSampler(
            PathSampler::EUnidirectional,
            scene,
            sensorSampler, emitterSampler, directSampler,
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
