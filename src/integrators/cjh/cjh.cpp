#include <mitsuba/bidir/util.h>
#include <mitsuba/core/plugin.h>

#include "cjh.h"
#include "cjh_sampler.h"

#include <mitsuba/core/sstream.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/plugin.h>


#ifdef __WINDOWS__
#include <io.h>
#include <ws2tcpip.h>
#include <mitsuba/core/getopt.h>
#else
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>
#define INVALID_SOCKET -1
#define SOCKET int
#endif


static SOCKET sock = INVALID_SOCKET;

MTS_NAMESPACE_BEGIN

class CJH : public Integrator {
public:
    CJH(const Properties &props) : Integrator(props) {
        m_config.x = props.getFloat("x");
        m_config.y = props.getFloat("y");

        m_config.u = props.getFloat("u");
        m_config.v = props.getFloat("v");

        m_config.direct1_1 = props.getFloat("direct1_1");
        m_config.direct1_2 = props.getFloat("direct1_2");

        m_config.bsdf1_1 = props.getFloat("bsdf1_1");
        m_config.bsdf1_2 = props.getFloat("bsdf1_2");

        m_config.direct2_1 = props.getFloat("direct2_1");
        m_config.direct2_2 = props.getFloat("direct2_2");

        m_config.bsdf2_1 = props.getFloat("bsdf2_1");
        m_config.bsdf2_2 = props.getFloat("bsdf2_2");
    }

    CJH(Stream *stream, InstanceManager *manager)
     : Integrator(stream, manager) {
        m_config = CJHConfiguration(stream);
        configure();
    }

    virtual ~CJH() { }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Integrator::serialize(stream, manager);
        m_config.serialize(stream);
    }

    bool preprocess(const Scene *scene, RenderQueue *queue,
            const RenderJob *job, int sceneResID, int sensorResID,
            int samplerResID) {
        Integrator::preprocess(scene, queue, job, sceneResID,
                sensorResID, samplerResID);

        return true;
    }

    void cancel() {}

    void runSample(Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {

        ref<FileStream> input = new FileStream("randoms.bin", FileStream::EReadOnly);
        ref<FileStream> output = new FileStream("luminance.bin", FileStream::ETruncWrite);

        unsigned int samples = input->readUInt();
        // printf("found %u samples\n", samples);

        Vector2i size = scene->getFilm()->getSize();

        for (unsigned int i = 0; i < samples; i++) {
            ref<CJHSampler> sensorSampler = new CJHSampler("sensor");
            ref<CJHSampler> emitterSampler = new CJHSampler("emitter");
            ref<CJHSampler> directSampler = new CJHSampler("direct");

            std::vector<Float> sensorSamples = {
                // t
                input->readFloat(),

                // u,v
                input->readFloat(),
                input->readFloat(),

                // x, y
                input->readFloat(),
                input->readFloat(),

                // direct 1
                input->readFloat(),
                input->readFloat(),

                // bsdf 1
                input->readFloat(),
                input->readFloat(),

                // direct 2
                input->readFloat(),
                input->readFloat(),

                // bsdf 2
                input->readFloat(),
                input->readFloat(),
            };

            std::vector<Float> emitterSamples = {};
            std::vector<Float> directSamples = {};

            sensorSampler->setSamples(sensorSamples);
            emitterSampler->setSamples(emitterSamples);
            directSampler->setSamples(directSamples);

            Properties props("path");
            props.setInteger("maxDepth", 3);
            props.setInteger("rrDepth", 9999);

            ref<SamplingIntegrator> pathTracer = static_cast<SamplingIntegrator *>(
                PluginManager::getInstance()->createObject(
                    MTS_CLASS(SamplingIntegrator),
                    props
                )
            );

            Sensor *sensor = scene->getSensor();

            Float timeSample = sensorSampler->next1D();
            Point2 apertureSample = sensorSampler->next2D();
            Point2 samplePos = sensorSampler->next2D();
            samplePos.x *= size[0];
            samplePos.y *= size[1];

            RayDifferential sensorRay;
            Spectrum value = sensor->sampleRayDifferential(
                sensorRay, samplePos, apertureSample, timeSample);

            RadianceQueryRecord rRec(scene, sensorSampler);

            uint32_t queryType = RadianceQueryRecord::ESensorRay \
                & ~RadianceQueryRecord::EOpacity;

            rRec.newQuery(queryType, sensor->getMedium());

            value *= pathTracer->Li(sensorRay, rRec);

            output->writeFloat(value.getLuminance());
        }

        // printf("RENDERED EVERYTHING!!\n");
    }

    bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {


        int listenPort = 65432;
        struct addrinfo hints, *servinfo, *p = NULL;
        memset(&hints, 0, sizeof(struct addrinfo));
        hints.ai_family = AF_UNSPEC;
        hints.ai_flags = AI_PASSIVE;
        hints.ai_socktype = SOCK_STREAM;
        char portName[8];
        sock = INVALID_SOCKET;

        snprintf(portName, sizeof(portName), "%i", listenPort);

        getaddrinfo(NULL, portName, &hints, &servinfo);
        p = servinfo;
        sock = socket(p->ai_family, p->ai_socktype, p->ai_protocol);
        int bindResult = bind(sock, p->ai_addr, (socklen_t) p->ai_addrlen);
        printf("bind: %d\n", bindResult);

        freeaddrinfo(servinfo);

        int listenResult = listen(sock, 1);
        printf("listen: %d\n", listenResult);

        int count = 1;
        while (true) {
            socklen_t addrlen = sizeof(sockaddr_storage);
            struct sockaddr_storage sockaddr;
            memset(&sockaddr, 0, addrlen);

            SOCKET newSocket = accept(sock, (struct sockaddr *) &sockaddr, &addrlen);
            if (newSocket == INVALID_SOCKET) {
                printf("INVALID!\n");
                break;

#if defined(__WINDOWS__)
                if (!running)
                    break;
#else
                if (errno == EINTR)
                    continue;
#endif
                SocketStream::handleError("none", "accept", EWarn);
                continue;
            }

            // printf("count: %d\n", count);
            runSample(scene, queue, job, sceneResID, sensorResID, samplerResID);
            send(newSocket, "!", 1, 0);
            count += 1;
            close(newSocket);
        }

        return true;
    }

    MTS_DECLARE_CLASS()
private:
    CJHConfiguration m_config;
};



MTS_IMPLEMENT_CLASS_S(CJH, false, Integrator)
MTS_EXPORT_PLUGIN(CJH, "CJH Testing Plugin");
MTS_NAMESPACE_END
