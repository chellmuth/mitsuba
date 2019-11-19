#include "neural_util.h"

#include <mitsuba/core/bitmap.h>

#include <cmath>

MTS_NAMESPACE_BEGIN

void imageFromVector(const fs::path &filename, std::vector<Float> values)
{
    const int totalSize = (int)values.size();
    const int side = (int)floorf(sqrtf(totalSize));

    ref<Bitmap> bitmap = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, {side, side});
    for (int i = 0; i < totalSize; i++) {
        const int x = i % side;
        const int y = (int)floorf(i / side);

        Spectrum value(values[i]);
        bitmap->setPixel({x, y}, value);
    }

    bitmap->write(Bitmap::EOpenEXR, filename);
}

MTS_NAMESPACE_END
