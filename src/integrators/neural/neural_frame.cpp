#include "neural_frame.h"

MTS_NAMESPACE_BEGIN

static const float M_TWO_PI = M_PI * 2.f;

static void cartesianToSpherical(Vector cartesian, float *phi, float *theta)
{
    *phi = atan2f(cartesian.z, cartesian.x);
    if (*phi < 0.f) {
        *phi += 2 * M_PI;
    }
    if (*phi == M_TWO_PI) {
        *phi = 0;
    }

    *theta = acosf(cartesian.y);
}

Frame constructNeuralFrame(const Vector &normal, const Intersection &its)
{
    Vector worldWi = its.toWorld(its.wi);

    Vector xAxis = normalize(cross(normal, worldWi));
    Vector zAxis = normalize(cross(normal, xAxis));

    return Frame(xAxis, normal, zAxis);
}

void testNeuralFrame(const Frame &neuralFrame, const Intersection &its)
{
    Vector worldWi = its.toWorld(its.wi);

    Vector neuralWi = neuralFrame.toLocal(worldWi);

    Float phi, theta;
    cartesianToSpherical(neuralWi, &phi, &theta);

    std::cout << "TEST NEURAL FRAME WI:" << std::endl
              << "  LOCAL: " << neuralWi.toString() << std::endl
              << "  PHI: " << phi << " (" << phi / M_TWO_PI << ")" << std::endl
              << "  THETA: " << theta << " (" << theta / (M_PI / 2.f) << ")" << std::endl
              << "END TEST" << std::endl;
}

MTS_NAMESPACE_END
