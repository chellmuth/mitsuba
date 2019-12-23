#include "neural_pdf.h"

#include <assert.h>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <string.h>

MTS_NAMESPACE_BEGIN

static const int PORT = 65432;
static const float M_TWO_PI = M_PI * 2.f;

bool NeuralPDF::connectToModel(int portOffset)
{
    struct sockaddr_in serv_addr;

    if ((m_socket = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    {
        printf("\n Socket creation error \n");
        return false;
    }

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(PORT + portOffset);

    // Convert IPv4 and IPv6 addresses from text to binary form
    if(inet_pton(AF_INET, "127.0.0.1", &serv_addr.sin_addr)<=0)
    {
        printf("\nInvalid address/ Address not supported \n");
        return false;
    }

    if (connect(m_socket, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
    {
        printf("\nConnection Failed \n");
        return false;
    }

    return true;
}

void NeuralPDF::sample(float *phi, float *theta, float *pdf, std::vector<float> &photonBundle) const
{
    const int count = 1;
    int hello[] = { 0, count };
    send(m_socket, &hello, sizeof(int) * 2, 0);

    float *photonData = photonBundle.data();
    send(m_socket, photonData, sizeof(float) * photonBundle.size(), 0);

    float buffer[3] = {0.f, 0.f, 0.f};
    int bytesRead = recv(m_socket, buffer, sizeof(buffer), 0);

    assert(bytesRead == sizeof(float) * 3);

    *phi = buffer[0] * M_TWO_PI;
    // *theta = (1.f - buffer[1]) * (M_PI / 2.f);
    *theta = buffer[1] * (M_PI / 2.f);
    // std::cout << "PDF: " << buffer[2] << std::endl;
    *pdf = buffer[2] / sinf(*theta) / (M_TWO_PI * M_PI / 2.f);

    if (*phi < 0 || *phi > M_TWO_PI || *theta < 0 || *theta > (M_PI/2.f) || *pdf < 0) {
        std::cout << "UHOH: " << *phi << " " << *theta << " " << *pdf
                  << " " << buffer[0] << " " << buffer[1] << " " << buffer[2]
                  << std::endl;

        sample(phi, theta, pdf, photonBundle);
    }
    // printf("(%f %f %f) (%f %f)\n", *phi, *theta, *pdf, buffer[0], buffer[1]);
}

float NeuralPDF::pdf(float phi, float theta, std::vector<float> &photonBundle) const
{
    const int count = 1;
    int hello[] = { 2, count };
    send(m_socket, &hello, sizeof(int) * 2, 0);

    float coordinates[] = { phi, theta };
    send(m_socket, &coordinates, sizeof(float) * 2, 0);

    float *photonData = photonBundle.data();
    send(m_socket, photonData, sizeof(float) * photonBundle.size(), 0);

    float buffer[1] = {0.f};
    int bytesRead = recv(m_socket, buffer, sizeof(buffer), 0);

    assert(bytesRead == sizeof(float) * 1);

    return buffer[0] / sinf(theta) / (M_TWO_PI * M_PI / 2.f);
}

std::vector<Float> NeuralPDF::batchEval(int size, std::vector<Float> photonBundle) const
{
    const int count = size * size;

    std::vector<Float> pdfs(count);

    int hello[] = { 1, count };
    send(m_socket, &hello, sizeof(int) * 2, 0);

    float *photonData = photonBundle.data();
    send(m_socket, photonData, sizeof(float) * photonBundle.size(), 0);

    float buffer[count];
    int bytesRead = recv(m_socket, buffer, sizeof(buffer), MSG_WAITALL);
    assert(bytesRead == 4 * count);

    const int thetaSteps = (int)sqrtf(count);
    const int phiSteps = (int)sqrtf(count);

    for (int i = 0; i < count; i++) {
        const int thetaStep = (int)floorf(i / phiSteps);
        const int phiStep = i % phiSteps;
        const int sourceIndex = thetaStep * phiSteps + phiStep;

        const float theta = (M_PI / 2.f) * (thetaStep + 0.5f) / thetaSteps;
        pdfs[i] = buffer[sourceIndex] / sinf(theta) / (M_TWO_PI * M_PI / 2.f);
    }

    return pdfs;
}

MTS_NAMESPACE_END
