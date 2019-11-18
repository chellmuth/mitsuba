#if !defined(__PHOTON_BUNDLE_H)
#define __PHOTON_BUNDLE_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/render/photon.h>

#include <vector>

MTS_NAMESPACE_BEGIN

class PhotonBundle {
public:
    PhotonBundle(Point point, Frame frame, int phiSteps, int thetaSteps);

    void splat(const Photon &photon);
    std::vector<Float> serializePDF();

private:
    const Point m_point;
    const Frame m_frame;

    const int m_phiSteps;
    const int m_thetaSteps;

    std::vector<Float> m_massLookup;
    Float m_totalMass;
};

MTS_NAMESPACE_END

#endif /* __PHOTON_BUNDLE_H */
