#if !defined(__COORDINATES_H)
#define __COORDINATES_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

const float M_TWO_PI = M_PI * 2.f;

inline Vector3 sphericalToCartesian(float phi, float theta)
{
    const float y = cosf(theta);
    const float x = sinf(theta) * cosf(phi);
    const float z = sinf(theta) * sinf(phi);

    return Vector3(x, y, z);
}

inline void cartesianToSpherical(Vector3 cartesian, float *phi, float *theta)
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

MTS_NAMESPACE_END

#endif /* __COORDINATES_H */
