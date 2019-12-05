#include "neural_frame.h"

MTS_NAMESPACE_BEGIN

Frame constructNeuralFrame(const Vector &normal, const Intersection &its)
{
    Vector worldWi = its.shFrame.toWorld(its.wi);

    Vector xAxis = normalize(cross(normal, worldWi));
    Vector zAxis = normalize(cross(normal, xAxis));

    return Frame(xAxis, normal, zAxis);
}

MTS_NAMESPACE_END
