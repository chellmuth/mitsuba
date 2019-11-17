#if !defined(__NEURAL_PDF_H)
#define __NEURAL_PDF_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

class NeuralPDF {
public:
    bool connectToModel(int portOffset);
    bool connectToModel() { return connectToModel(0); }

private:
    int m_socket;
};

MTS_NAMESPACE_END

#endif /* __NEURAL_PDF_H */
