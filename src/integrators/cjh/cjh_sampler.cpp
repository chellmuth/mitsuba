#include "cjh_sampler.h"

MTS_NAMESPACE_BEGIN


CJHSampler::CJHSampler() : CJHSampler(Properties()) {}

CJHSampler::CJHSampler(const Properties &props) : Sampler(props) {
    m_random = new Random();
    m_sampleIndex = 0;
}

CJHSampler::CJHSampler(Stream *stream, InstanceManager *manager)
 : Sampler(stream, manager) {
    m_random = static_cast<Random *>(manager->getInstance(stream));
    m_sampleIndex = 0;
}

void CJHSampler::setSamples(const std::vector<float> &samples) {
    m_samples = samples;
    m_sampleIndex = 0;
}

void CJHSampler::serialize(Stream *stream, InstanceManager *manager) const {
    // needs updating

    Sampler::serialize(stream, manager);
    manager->serialize(stream, m_random.get());
}

ref<Sampler> CJHSampler::clone() {
    // needs updating

    ref<CJHSampler> sampler = new CJHSampler();
    sampler->m_random = new Random(m_random);
    return sampler.get();
}

Float CJHSampler::next1D() {
    printf("Request 1D\n");
    Assert(m_sampleIndex < m_samples.size());
    return m_samples[m_sampleIndex++];
}

Point2 CJHSampler::next2D() {
    printf("Request 2D\n");
    Assert(m_sampleIndex + 1 < m_samples.size());
    Float value1 = m_samples[m_sampleIndex++];
    Float value2 = m_samples[m_sampleIndex++];
    return Point2(value1, value2);
}

std::string CJHSampler::toString() const {
    std::ostringstream oss;
    oss << "CJHSampler[" << endl
        << "]";
    return oss.str();
}


MTS_IMPLEMENT_CLASS_S(CJHSampler, false, Sampler)
MTS_NAMESPACE_END
