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

    float runSample(Scene *scene, float *samples, int sampleCount) {
        // ref<FileStream> output = new FileStream("luminance.bin", FileStream::ETruncWrite);

        Vector2i size = scene->getFilm()->getSize();

        for (unsigned int i = 0; i < sampleCount; i++) {
            ref<CJHSampler> sensorSampler = new CJHSampler("sensor");
            ref<CJHSampler> emitterSampler = new CJHSampler("emitter");
            ref<CJHSampler> directSampler = new CJHSampler("direct");

            int sampleIndex = 0;
            std::vector<Float> sensorSamples = {
                // t
                samples[sampleIndex++],

                // u,v
                samples[sampleIndex++],
                samples[sampleIndex++],

                // x, y
                samples[sampleIndex++],
                samples[sampleIndex++],

                // direct 1
                samples[sampleIndex++],
                samples[sampleIndex++],

                // bsdf 1
                samples[sampleIndex++],
                samples[sampleIndex++],

                // direct 2
                samples[sampleIndex++],
                samples[sampleIndex++],

                // bsdf 2
                samples[sampleIndex++],
                samples[sampleIndex++]
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

            // output->writeFloat(value.getLuminance());
            return value.getLuminance();
        }

        // printf("RENDERED EVERYTHING!!\n");
    }

    bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {

#if defined(__WINDOWS__)
        WSADATA wsaData;
        if (WSAStartup(MAKEWORD(2,2), &wsaData))
            SLog(EError, "Could not initialize WinSock2!");
        if (LOBYTE(wsaData.wVersion) != 2 || HIBYTE(wsaData.wVersion) != 2)
            SLog(EError, "Could not find the required version of winsock.dll!");
#endif

        int listenPort = 65432;
        struct addrinfo hints, *servinfo, *p = NULL;
        memset(&hints, 0, sizeof(struct addrinfo));
        hints.ai_family = AF_INET;
        hints.ai_flags = AI_PASSIVE;
        hints.ai_socktype = SOCK_STREAM;
        // hints.ai_protocol = IPPROTO_TCP;
        char portName[8];
        sock = INVALID_SOCKET;

        snprintf(portName, sizeof(portName), "%i", listenPort);

        getaddrinfo(NULL, portName, &hints, &servinfo);
        p = servinfo;
        sock = socket(p->ai_family, p->ai_socktype, p->ai_protocol);

        int enable = 1;
        if (setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, (const char *)&enable, sizeof(int)) < 0) {
            printf("setsockopt(SO_REUSEADDR) failed");
            return true;
        }
        // int disable = 0;
        // if (setsockopt(sock, SOL_SOCKET, SO_LINGER, (const char *)&disable, sizeof(int)) < 0) {
        //     printf("setsockopt(SO_LINGER) failed");
        //     return true;
        // }

        int bindResult = bind(sock, p->ai_addr, (socklen_t) p->ai_addrlen);
        printf("bind: %d\n", bindResult);

        freeaddrinfo(servinfo);

        int listenResult = listen(sock, 1);
        printf("listen: %d\n", listenResult);
        printf("OKAY TO RUN MATLAB FUNCTIONS\n");

        int count = 1;
        while (true) {
            socklen_t addrlen = sizeof(sockaddr_storage);
            struct sockaddr_storage sockaddr;
            memset(&sockaddr, 0, addrlen);

            SOCKET newSocket = accept(sock, (struct sockaddr *) &sockaddr, &addrlen);
            if (newSocket == INVALID_SOCKET) {
                printf("INVALID!\n");
                break;

                SocketStream::handleError("none", "accept", EWarn);
                continue;
            }

            ssize_t bytesRead;
            int incomingSamples;
            // printf("about to read\n");
            // bytesRead = read(newSocket, &incomingSamples, sizeof(incomingSamples));
            bytesRead = recv(newSocket, (char *)&incomingSamples, sizeof(incomingSamples), 0);

            // printf("initial read done (%d)\n", bytesRead);
            // printf("read: %d\n", incomingSamples);
            for (int i = 0; i < incomingSamples; i++) {
                float samples[13];
                // bytesRead = read(newSocket, &samples, sizeof(samples), 0);
                bytesRead = recv(newSocket, (char *)&samples, sizeof(samples), 0);
                // printf("bytes read (%d/%d)\n", bytesRead, sizeof(samples));
                if (bytesRead != sizeof(samples)) {
                    printf("Only read %d bytes\n", bytesRead);
                    break;
                }
                // printf("sample1: %f, sample5: %f, sample13: %f\n", samples[0], samples[4], samples[12]);
                float result = runSample(scene, samples, incomingSamples);

                // printf("result: %f (size %d)\n", result, sizeof(result));
                send(newSocket, (char *)&result, 4, 0);
                count += 1;
            }

            // shutdown(newSocket, 2);
#if defined(__WINDOWS__)
            int closeCode = closesocket(newSocket);
            if (closeCode != 0) {
                printf("BAD CLOSE CODE %d\n", closeCode);
                break;
            }
#else
            close(newSocket);
#endif
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
