#if !defined(__CJH_H)
#define __CJH_H

/* #include <mitsuba/bidir/pathsampler.h> */
/* #include <mitsuba/core/bitmap.h> */

MTS_NAMESPACE_BEGIN

struct CJHConfiguration {
    Float x;
    Float y;

    Float direct1_1;
    Float direct1_2;

    Float bsdf1_1;
    Float bsdf1_2;

    inline CJHConfiguration() { }

    inline CJHConfiguration(Stream *stream) {}

    inline void serialize(Stream *stream) const {}
};

MTS_NAMESPACE_END

#endif /* __CJH_H */
