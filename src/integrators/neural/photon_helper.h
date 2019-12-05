#if !defined(__PHOTON_HELPER_H)
#define __PHOTON_HELPER_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/imageblock.h>
#include <mitsuba/render/photonmap.h>

MTS_NAMESPACE_BEGIN

void gatherPhotons(
    ref<PhotonMap> globalPhotonMap,
    int x, int y,
    const RadianceQueryRecord &rRec,
    bool flipNormal,
    const ImageBlock *block = nullptr,
    const int identifier = -1
);

MTS_NAMESPACE_END

#endif /* __PHOTON_HELPER_H */
