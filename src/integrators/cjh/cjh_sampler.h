#if !defined(__CJH_SAMPLER_H)
#define __CJH_SAMPLER_H

#include <mitsuba/render/sampler.h>

MTS_NAMESPACE_BEGIN

class CJHSampler : public Sampler {
public:
    CJHSampler();
    CJHSampler(const Properties &props);
    CJHSampler(Stream *stream, InstanceManager *manager);

    void serialize(Stream *stream, InstanceManager *manager) const;
    ref<Sampler> clone();
    void generate(const Point2i &pos);
    std::string toString() const;

    Float next1D();
    Point2 next2D();

    MTS_DECLARE_CLASS()

private:
    ref<Random> m_random;
};

MTS_NAMESPACE_END

#endif /* __CJH_SAMPLER_H */

