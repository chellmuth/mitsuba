#include "cjh_sampler.h"

MTS_NAMESPACE_BEGIN


CJHSampler::CJHSampler() : CJHSampler(Properties()) {}

CJHSampler::CJHSampler(const Properties &props) : Sampler(props) {
    /* Number of samples per pixel when used with a sampling-based integrator */
    m_sampleCount = props.getSize("sampleCount", 4);
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
    sampler->m_sampleCount = m_sampleCount;
    sampler->m_random = new Random(m_random);
    for (size_t i=0; i<m_req1D.size(); ++i)
        sampler->request1DArray(m_req1D[i]);
    for (size_t i=0; i<m_req2D.size(); ++i)
        sampler->request2DArray(m_req2D[i]);
    return sampler.get();
}

void CJHSampler::generate(const Point2i &) {
    for (size_t i=0; i<m_req1D.size(); i++)
        for (size_t j=0; j<m_sampleCount * m_req1D[i]; ++j)
            m_sampleArrays1D[i][j] = m_random->nextFloat();
    for (size_t i=0; i<m_req2D.size(); i++)
        for (size_t j=0; j<m_sampleCount * m_req2D[i]; ++j)
            m_sampleArrays2D[i][j] = Point2(
                m_random->nextFloat(),
                m_random->nextFloat());
    m_sampleIndex = 0;
    m_dimension1DArray = m_dimension2DArray = 0;
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
        << "  sampleCount = " << m_sampleCount << endl
        << "]";
    return oss.str();
}


MTS_IMPLEMENT_CLASS_S(CJHSampler, false, Sampler)
MTS_NAMESPACE_END
