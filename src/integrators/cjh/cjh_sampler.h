#if !defined(__CJH_SAMPLER_H)
#define __CJH_SAMPLER_H

#include <mitsuba/render/sampler.h>

MTS_NAMESPACE_BEGIN

class CJHSampler : public Sampler {
public:
    CJHSampler();
    CJHSampler(const Properties &props);
    CJHSampler(Stream *stream, InstanceManager *manager);

    void setSamples(const std::vector<float> &samples);

    void serialize(Stream *stream, InstanceManager *manager) const;
    ref<Sampler> clone();
    std::string toString() const;

    void generate(const Point2i &pos) { Log(EError, "generate(): Unsupported!"); }
    void request1DArray(size_t size) { Log(EError, "request1DArray(): Unsupported!"); }
    void request2DArray(size_t size) { Log(EError, "request2DArray(): Unsupported!"); }

    Float next1D();
    Point2 next2D();

    MTS_DECLARE_CLASS()

private:
    ref<Random> m_random;


    std::vector<Float> m_samples;
    unsigned int m_sampleIndex;
};

MTS_NAMESPACE_END

#endif /* __CJH_SAMPLER_H */
