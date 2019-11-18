#include "photon_bundle.h"

#include <cmath>

MTS_NAMESPACE_BEGIN

static const float M_TWO_PI = M_PI * 2.f;

// (0, 1, 0) is up
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

PhotonBundle::PhotonBundle(Point point, Frame frame, int phiSteps, int thetaSteps)
    : m_point(point), m_frame(frame), m_phiSteps(phiSteps), m_thetaSteps(thetaSteps),
      m_massLookup(m_phiSteps * m_thetaSteps, 0.f), m_totalMass(0.f)
{
    std::cout << "PHOTON BUNDLE CONSTRUCTOR" << std::endl;
}

void PhotonBundle::splat(const Photon &photon)
{
    Point source = photon.getSource();
    Vector worldDirection = source - m_point;
    Vector localDirection = m_frame.toLocal(worldDirection);

    float phi, theta;
    cartesianToSpherical(localDirection, &phi, &theta);

    if (theta > M_PI / 2.f) { return; }

    const int phiStep = (int)floorf(phi / (M_TWO_PI / m_phiSteps));
    const int thetaStep = (int)floorf(theta / (M_PI / 2.f / m_thetaSteps));

    assert(phiStep >= 0);
    assert(thetaStep >= 0);
    assert(phiStep < m_phiSteps);
    assert(thetaStep < m_thetaSteps);

    float r, g, b;
    Spectrum power = photon.getPower();
    power.toLinearRGB(r, g, b);
    float mass = r + g + b;

    m_massLookup[thetaStep * m_phiSteps + phiStep] += mass;
    m_totalMass += mass;
}

std::vector<Float> PhotonBundle::serializePDF()
{
    std::vector<Float> pdf(m_phiSteps * m_thetaSteps, 0.f);

    if (m_totalMass == 0.f) {
        return pdf;
    }

    for (int i = 0; i < m_phiSteps * m_thetaSteps; i++) {
        pdf[i] = m_massLookup[i] / m_totalMass;
    }

    return pdf;
}

MTS_NAMESPACE_END
