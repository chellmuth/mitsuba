#include "photon_helper.h"

#include <mitsuba/core/fstream.h>
#include <mitsuba/core/kdtree.h>
#include <mitsuba/render/gatherproc.h>
#include <mitsuba/render/photon.h>

MTS_NAMESPACE_BEGIN

void gatherPhotons(
    const std::string &prefix,
    ref<PhotonMap> globalPhotonMap,
    int x, int y,
    const Intersection &its,
    bool flipNormal,
    const ImageBlock *block,
    const int identifier
) {
    typedef PointKDTree<Photon>      PhotonTree;
    typedef PhotonTree::SearchResult SearchResult;

    const size_t maxPhotons = 10;
    SearchResult *results = static_cast<SearchResult *>(
        alloca((maxPhotons + 1) * sizeof(SearchResult)));

    size_t resultCount = globalPhotonMap->nnSearch(its.p, maxPhotons, results);
    // Log(EInfo, "Photons returned: %i", resultCount);

    std::cout << "INTERSECTION RECORD:" << std::endl;
    std::cout << its.p.toString() << std::endl;
    // std::cout << its.geoFrame.n.toString() << std::endl;
    std::cout << "SH FRAME.N " << its.shFrame.n.toString() << std::endl;
    std::cout << "TO WORLD SH FRAME.N " << its.toWorld(its.shFrame.n).toString() << std::endl;
    std::cout << "TO LOCAL SH FRAME.N " << its.toLocal(its.shFrame.n).toString() << std::endl;
    std::cout << "WI LOCAL: " << its.wi.toString() << " WI WORLD: " << its.toWorld(its.wi).toString() << std::endl;

    std::ostringstream oss;
    if (identifier < 0) {
        oss << prefix << "_" << x << "_" << y << ".bin";
    } else {
        auto offset = Vector2i(block->getOffset());
        oss << "results/" << prefix << "_" << identifier << "-block_" << offset.x << "x" << offset.y << ".bin";
    }

    ref<FileStream> fileStream = new FileStream(oss.str(), FileStream::ETruncWrite);

    float intersectionBuffer[] = {
        its.p.x,
        its.p.y,
        its.p.z,
        its.shFrame.n.x * (flipNormal ? -1.f : 1.f),
        its.shFrame.n.y * (flipNormal ? -1.f : 1.f),
        its.shFrame.n.z * (flipNormal ? -1.f : 1.f),
        its.toWorld(its.wi).x,
        its.toWorld(its.wi).y,
        its.toWorld(its.wi).z,
    };
    fileStream->write(&intersectionBuffer, 9 * sizeof(float));

    uint32_t countBuffer = resultCount;
    fileStream->write(&countBuffer, sizeof(uint32_t));

    for (size_t i = 0; i < resultCount; i++) {
        const SearchResult &searchResult = results[i];
        const PhotonMap &photonMap = (*globalPhotonMap.get());
        const Photon &photon = photonMap[searchResult.index];

        Point position = photon.getPosition();
        Point source = photon.getSource();

        float r, g, b;
        Spectrum power = photon.getPower();
        power.toLinearRGB(r, g, b);

        float photonBuffer[] = {
            position.x,
            position.y,
            position.z,

            source.x,
            source.y,
            source.z,

            r, g, b
        };

        fileStream->write(photonBuffer, sizeof(float) * 9);
    }

    fileStream->close();
}

MTS_NAMESPACE_END
