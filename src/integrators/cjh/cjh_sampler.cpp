#include "cjh_sampler.h"

MTS_NAMESPACE_BEGIN


CJHSampler::CJHSampler() : CJHSampler(Properties()) {}

CJHSampler::CJHSampler(const Properties &props) : Sampler(props) {
    m_random = new Random();
}

CJHSampler::CJHSampler(Stream *stream, InstanceManager *manager)
 : Sampler(stream, manager) {
    m_random = static_cast<Random *>(manager->getInstance(stream));
}

void CJHSampler::serialize(Stream *stream, InstanceManager *manager) const {
    Sampler::serialize(stream, manager);
    manager->serialize(stream, m_random.get());
}

ref<Sampler> CJHSampler::clone() {
    ref<CJHSampler> sampler = new CJHSampler();
    sampler->m_random = new Random(m_random);
    return sampler.get();
}

Float CJHSampler::next1D() {
    return m_random->nextFloat();
}

Point2 CJHSampler::next2D() {
    Float value1 = m_random->nextFloat();
    Float value2 = m_random->nextFloat();
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
