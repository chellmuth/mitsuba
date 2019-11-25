#include "neural_frame.h"
#include "neural_pdf.h"
#include "neural_util.h"
#include "photon_bundle.h"

#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>
#include <mitsuba/core/statistics.h>

#include <mitsuba/core/fstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/gatherproc.h>
#include <mitsuba/render/photon.h>
#include <mitsuba/render/photonmap.h>

#include <cmath>

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength("Path tracer", "Average path length", EAverage);

static Vector sphericalToCartesian(float phi, float theta)
{
    const float y = cosf(theta);
    const float x = sinf(theta) * cosf(phi);
    const float z = sinf(theta) * sinf(phi);

    return Vector3(x, y, z);
}

class NeuralIntegrator : public MonteCarloIntegrator {
public:
    typedef PointKDTree<Photon>      PhotonTree;
    typedef PhotonTree::SearchResult SearchResult;

    NeuralIntegrator(const Properties &props)
        : MonteCarloIntegrator(props), m_neuralPDF()
        {
            m_x = props.getInteger("x");
            m_y = props.getInteger("y");

            if (m_x >= 0 && m_y >= 0) {
                Log(EInfo, "Hunting pixel (%i, %i)", m_x, m_y);
            } else {
                Log(EInfo, "Full render");
            }

            m_globalPhotons = props.getSize("globalPhotons", 10000);
        }

    /// Unserialize from a binary data stream
    NeuralIntegrator(Stream *stream, InstanceManager *manager)
        : MonteCarloIntegrator(stream, manager) { }

    bool preprocess(
        const Scene *scene,
        RenderQueue *queue,
        const RenderJob *job,
        int sceneResID,
        int sensorResID,
        int samplerResID
    ) {
        if (!MonteCarloIntegrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID)) {
            Log(EError, "Base class error");
            return false;
        }

        m_neuralPDF.connectToModel();

        if (m_globalPhotonMap.get() == NULL && m_globalPhotons > 0) {
            /* Generate the global photon map */
            ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
                GatherPhotonProcess::ESurfacePhotons,
                m_globalPhotons, // count
                0, // granularity for parallelization
                m_maxDepth - 1,
                m_rrDepth,
                true, // gather locally
                true, // auto-cancel if not enough photons are generated
                job
            );

            proc->bindResource("scene", sceneResID);
            proc->bindResource("sensor", sensorResID);
            proc->bindResource("sampler", samplerResID);

            ref<Scheduler> sched = Scheduler::getInstance();
            m_proc = proc;
            sched->schedule(proc);
            sched->wait(proc);
            m_proc = NULL;

            if (proc->getReturnStatus() != ParallelProcess::ESuccess) {
                return false;
            }

            ref<PhotonMap> globalPhotonMap = proc->getPhotonMap();
            if (globalPhotonMap->isFull()) {
                Log(EDebug, "Global photon map full. Shot " SIZE_T_FMT " particles, excess photons due to parallelism: "
                    SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

                // m_globalPhotonMap = globalPhotonMap;
                // m_globalPhotonMap->setScaleFactor(1 / (Float) proc->getShotParticles());
                // m_globalPhotonMap->build();
                // m_globalPhotonMapID = sched->registerResource(m_globalPhotonMap);
            }

            m_globalPhotonMap = globalPhotonMap;
            m_globalPhotonMap->setScaleFactor(1 / (Float) proc->getShotParticles());
            m_globalPhotonMap->build();
        }

        return true;
    }

    Spectrum neuralSample(const BSDF *bsdf, BSDFSamplingRecord &bRec, Float &pdf, bool debugPixel) const {
        const Intersection &its = bRec.its;

        // std::cout << "INTERSECTION:" << std::endl;
        // std::cout << its.toString() << std::endl;

        const size_t maxPhotons = 100;
        SearchResult *results = static_cast<SearchResult *>(
            alloca((maxPhotons + 1) * sizeof(SearchResult)));

        size_t resultCount = m_globalPhotonMap->nnSearch(its.p, maxPhotons, results);

        bool flippedNormal = false;
        Vector normal = its.shFrame.n;
        if (Frame::cosTheta(bRec.wi) < 0.f) {
            flippedNormal = true;
            normal *= -1.f;
        }
        Frame neuralFrame = constructNeuralFrame(normal, its.wi);
        PhotonBundle bundle(its.p, neuralFrame, 10, 10);

        for (size_t i = 0; i < resultCount; i++) {
            const SearchResult &searchResult = results[i];
            const PhotonMap &photonMap = (*m_globalPhotonMap.get());
            const Photon &photon = photonMap[searchResult.index];

            bundle.splat(photon);
        }

        std::vector<Float> photonBundle = bundle.serialized();
        float phi, theta, pdf2;
        m_neuralPDF.sample(&phi, &theta, &pdf2, photonBundle);

        if (debugPixel) {
            {
                std::ostringstream oss;
                oss << "batch_" << m_x << "_" << m_y << ".exr";

                std::vector<Float> batchedResult = m_neuralPDF.batchEval(400, photonBundle);
                imageFromVector(oss.str(), batchedResult);
            }
            {
                std::ostringstream oss;
                oss << "grid_" << m_x << "_" << m_y << ".bin";

                ref<FileStream> fileStream = new FileStream(oss.str(), FileStream::ETruncWrite);
                fileStream->write(photonBundle.data(), sizeof(float) * 100);
                fileStream->close();
            }
        }

        Vector localDirection = sphericalToCartesian(phi, theta);
        // std::cout << "PHI: " << phi << " THETA: " << theta << std::endl;
        // std::cout << "MY FORMAT: " << localDirection.toString() << std::endl;
        localDirection = Vector(localDirection.x, localDirection.z, localDirection.y);
        // std::cout << "MITSUBA FORMAT: " << localDirection.toString() << std::endl;

        if (flippedNormal) {
            // std::cout << "FLIPPING NORMAL" << std::endl;
            localDirection.z *= -1.f;
        }

        // std::cout << "FLIPPED MITSUBA: " << localDirection.toString() << std::endl;
        // std::cout << "PDF: " << pdf2 << std::endl;

        bRec.wo = localDirection;
        bRec.eta = 1.f;
        bRec.sampledComponent = 0;
        bRec.sampledType = BSDF::EDiffuseReflection;

        pdf = pdf2;
        return bsdf->eval(bRec);
    }

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        return Li(r, rRec, false);
    }

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec, bool debugPixel) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        RayDifferential ray(r);
        Spectrum Li(0.0f);
        bool scattered = false;

        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        rRec.rayIntersect(ray);
        ray.mint = Epsilon;

        Spectrum throughput(1.0f);
        Float eta = 1.0f;

        while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
            if (!its.isValid()) {
                /* If no intersection could be found, potentially return
                   radiance from a environment luminaire if it exists */
                if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                    && (!m_hideEmitters || scattered))
                    Li += throughput * scene->evalEnvironment(ray);
                break;
            }

            const BSDF *bsdf = its.getBSDF(ray);

            /* Possibly include emitted radiance if requested */
            if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                && (!m_hideEmitters || scattered))
                Li += throughput * its.Le(-ray.d);

            /* Include radiance from a subsurface scattering model if requested */
            if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
                Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

            if ((rRec.depth >= m_maxDepth && m_maxDepth > 0)
                || (m_strictNormals && dot(ray.d, its.geoFrame.n)
                    * Frame::cosTheta(its.wi) >= 0)) {

                /* Only continue if:
                   1. The current path length is below the specifed maximum
                   2. If 'strictNormals'=true, when the geometric and shading
                      normals classify the incident direction to the same side */
                break;
            }

            /* ==================================================================== */
            /*                     Direct illumination sampling                     */
            /* ==================================================================== */

            /* Estimate the direct illumination if this is requested */
            DirectSamplingRecord dRec(its);

            if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
                (bsdf->getType() & BSDF::ESmooth)) {
                Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D());
                if (!value.isZero()) {
                    const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

                    /* Allocate a record for querying the BSDF */
                    BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

                    /* Evaluate BSDF * cos(theta) */
                    const Spectrum bsdfVal = bsdf->eval(bRec);

                    /* Prevent light leaks due to the use of shading normals */
                    if (!bsdfVal.isZero() && (!m_strictNormals
                            || dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {

                        /* Calculate prob. of having generated that direction
                           using BSDF sampling */
                        Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                            ? bsdf->pdf(bRec) : 0;

                        /* Weight using the power heuristic */
                        Float weight = miWeight(dRec.pdf, bsdfPdf);
                        Li += throughput * value * bsdfVal * weight;

                        // std::cout << "EMITTER HIT!" << std::endl;
                    }
                }
            }

            /* ==================================================================== */
            /*                            BSDF sampling                             */
            /* ==================================================================== */

            /* Sample BSDF * cos(theta) */
            Float bsdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);

            // Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
            Spectrum bsdfWeight = neuralSample(bsdf, bRec, bsdfPdf, debugPixel && rRec.depth == 1);
            // for (int i = 0; i < 20; i++) {
            //      bsdfWeight = neuralSample(bsdf, bRec, bsdfPdf, debugPixel && rRec.depth == 1);
            //      const Vector wo = its.toWorld(bRec.wo);
            //      ray = Ray(its.p, wo, ray.time);
            //      if (scene->rayIntersect(ray, its)) {
            //          std::cout << "HIT!" << std::endl;
            //      } else {
            //          std::cout << "MISS!" << std::endl;
            //      }
            // }

            if (bsdfWeight.isZero())
                break;

            scattered |= bRec.sampledType != BSDF::ENull;

            /* Prevent light leaks due to the use of shading normals */
            const Vector wo = its.toWorld(bRec.wo);
            Float woDotGeoN = dot(its.geoFrame.n, wo);
            if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                break;

            bool hitEmitter = false;
            Spectrum value;

            /* Trace a ray in this direction */
            ray = Ray(its.p, wo, ray.time);
            if (scene->rayIntersect(ray, its)) {
                // std::cout << "HIT!" << std::endl;
                /* Intersected something - check if it was a luminaire */
                if (its.isEmitter()) {
                    value = its.Le(-ray.d);
                    dRec.setQuery(ray, its);
                    hitEmitter = true;
                }
            } else {
                // std::cout << "MISSED INTERSECTION!" << std::endl;
                /* Intersected nothing -- perhaps there is an environment map? */
                const Emitter *env = scene->getEnvironmentEmitter();

                if (env) {
                    if (m_hideEmitters && !scattered)
                        break;

                    value = env->evalEnvironment(ray);
                    if (!env->fillDirectSamplingRecord(dRec, ray))
                        break;
                    hitEmitter = true;
                } else {
                    break;
                }
            }

            /* Keep track of the throughput and relative
               refractive index along the path */
            throughput *= bsdfWeight;
            eta *= bRec.eta;

            /* If a luminaire was hit, estimate the local illumination and
               weight using the power heuristic */
            if (hitEmitter &&
                (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
                /* Compute the prob. of generating that direction using the
                   implemented direct illumination sampling technique */
                const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
                    scene->pdfEmitterDirect(dRec) : 0;
                Li += throughput * value * miWeight(bsdfPdf, lumPdf);
            }

            /* ==================================================================== */
            /*                         Indirect illumination                        */
            /* ==================================================================== */

            /* Set the recursive query type. Stop if no surface was hit by the
               BSDF sample or if indirect illumination was not requested */
            if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
                break;
            rRec.type = RadianceQueryRecord::ERadianceNoEmission;

            if (rRec.depth++ >= m_rrDepth) {
                /* Russian roulette: try to keep path weights equal to one,
                   while accounting for the solid angle compression at refractive
                   index boundaries. Stop with at least some probability to avoid
                   getting stuck (e.g. due to total internal reflection) */

                Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
                if (rRec.nextSample1D() >= q)
                    break;
                throughput /= q;
            }
        }

        /* Store statistics */
        avgPathLength.incrementBase();
        avgPathLength += rRec.depth;

        return Li;
    }

    inline Float miWeight(Float pdfA, Float pdfB) const {
        pdfA *= pdfA;
        pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }

    void renderBlock(const Scene *scene,
            const Sensor *sensor, Sampler *sampler, ImageBlock *block,
            const bool &stop, const std::vector< TPoint2<uint8_t> > &points) const {

        Float diffScaleFactor = 1.0f /
            std::sqrt((Float) sampler->getSampleCount());

        bool needsApertureSample = sensor->needsApertureSample();
        bool needsTimeSample = sensor->needsTimeSample();

        RadianceQueryRecord rRec(scene, sampler);
        Point2 apertureSample(0.5f);
        Float timeSample = 0.5f;
        RayDifferential sensorRay;

        block->clear();

        uint32_t queryType = RadianceQueryRecord::ESensorRay;

        if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we don't have to */
            queryType &= ~RadianceQueryRecord::EOpacity;

        for (size_t i = 0; i<points.size(); ++i) {
            Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());

            bool debugPixel = false;
            if (m_x == offset.x && m_y == offset.y) {
                debugPixel = true;
                Log(EInfo, "Rendering Pixel (%i, %i)", m_x, m_y);
            } else if (m_x != -1 && m_y != -1) {
                continue;
            }

            if (stop)
                break;

            sampler->generate(offset);

            for (size_t j = 0; j<sampler->getSampleCount(); j++) {
                rRec.newQuery(queryType, sensor->getMedium());
                Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

                if (needsApertureSample)
                    apertureSample = rRec.nextSample2D();
                if (needsTimeSample)
                    timeSample = rRec.nextSample1D();

                Spectrum spec = sensor->sampleRayDifferential(
                    sensorRay, samplePos, apertureSample, timeSample);

                sensorRay.scaleDifferential(diffScaleFactor);

                spec *= Li(sensorRay, rRec, debugPixel);
                block->put(samplePos, spec, rRec.alpha);
                sampler->advance();
            }
        }
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        MonteCarloIntegrator::serialize(stream, manager);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "NeuralIntegrator[" << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  strictNormals = " << m_strictNormals << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()

private:
    int m_x, m_y;

    ref<ParallelProcess> m_proc;
    ref<PhotonMap> m_globalPhotonMap;
    size_t m_globalPhotons;

    NeuralPDF m_neuralPDF;
};

MTS_IMPLEMENT_CLASS_S(NeuralIntegrator, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(NeuralIntegrator, "Neural integrator");
MTS_NAMESPACE_END
