#include "neural_frame.h"

MTS_NAMESPACE_BEGIN

Frame constructNeuralFrame(Vector normal, Vector wi)
{
    Vector xAxis = normalize(cross(normal, wi));
    Vector zAxis = normalize(cross(normal, xAxis));

    return Frame(xAxis, normal, zAxis);
}

MTS_NAMESPACE_END
