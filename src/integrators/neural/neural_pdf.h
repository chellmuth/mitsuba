#if !defined(__NEURAL_PDF_H)
#define __NEURAL_PDF_H

#include <mitsuba/mitsuba.h>

#include <vector>

MTS_NAMESPACE_BEGIN

class NeuralPDF {
public:
    bool connectToModel(int portOffset);
    bool connectToModel() { return connectToModel(0); }

    void sample(float *x, float *y, float *pdf, std::vector<float> &photonBundle) const;
    float pdf(float phi, float theta, std::vector<float> &photonBundle) const;
    std::vector<Float> batchEval(int size, std::vector<Float> photonBundle) const;

private:
    int m_socket;
};

MTS_NAMESPACE_END

#endif /* __NEURAL_PDF_H */
