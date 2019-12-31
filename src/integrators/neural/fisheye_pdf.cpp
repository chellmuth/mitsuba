#include "fisheye_pdf.h"

#include "coordinates.h"

#include <iostream>

MTS_NAMESPACE_BEGIN

FisheyePDF::FisheyePDF(ref<Bitmap> bitmap)
    : m_bitmap(bitmap), m_sum(0.f)
{
    const int width = m_bitmap->getWidth();
    const int height = m_bitmap->getHeight();
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            Point2i point(col, row);
            m_sum += m_bitmap->getPixel(point).average();
        }
    }
}

void FisheyePDF::sample(float *x, float *y, float *pdf) const
{
    float cdf = 0.f;
    float lastPixel = 0.f;
    float currentPixel = 0.f;

    const int width = m_bitmap->getWidth();
    const int height = m_bitmap->getHeight();
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            lastPixel = currentPixel;

            Point2i point(col, row);
            currentPixel = m_bitmap->getPixel(point).average();

            cdf += currentPixel / m_sum;

            if (xi <= cdf) {
                *x = col;
                *y = row;
                *pdf = (currentPixel - lastPixel) / m_sum;
                return;
            }
        }
    }
}

float FisheyePDF::pdf(float phi, float theta) const
{
    const int x = (int)floorf(phi * M_TWO_PI / m_bitmap->getWidth());
    const int y = (int)floorf(theta * (M_PI / 2.f) / m_bitmap->getHeight());

    Point2i point(x, y);
    Spectrum result = m_bitmap->getPixel(point);
    return result.average() / m_sum;
}


MTS_NAMESPACE_END
