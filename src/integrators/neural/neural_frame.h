#if !defined(__NEURAL_FRAME_H)
#define __NEURAL_FRAME_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/core/vector.h>

MTS_NAMESPACE_BEGIN

Frame constructNeuralFrame(const Vector &normal, const Intersection &its);

void testNeuralFrame(const Frame &neuralFrame, const Intersection &its);

MTS_NAMESPACE_END

#endif /* __NEURAL_FRAME_H */
