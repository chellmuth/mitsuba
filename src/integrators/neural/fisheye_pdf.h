#if !defined(__FISHEYE_PDF_H)
#define __FISHEYE_PDF_H

#include <mitsuba/core/bitmap.h>
#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

class FisheyePDF {
public:
    FisheyePDF(ref<Bitmap> bitmap);
    void sample(float *x, float *y, float *pdf) const;
    float pdf(float phi, float theta) const;

private:
    float m_sum;
    ref<Bitmap> m_bitmap;
};

MTS_NAMESPACE_END

#endif /* __FISHEYE_PDF_H */
