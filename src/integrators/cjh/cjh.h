#if !defined(__CJH_H)
#define __CJH_H

/* #include <mitsuba/bidir/pathsampler.h> */
/* #include <mitsuba/core/bitmap.h> */

MTS_NAMESPACE_BEGIN

struct CJHConfiguration {
    inline CJHConfiguration() { }

    inline CJHConfiguration(Stream *stream) {}

    inline void serialize(Stream *stream) const {}
};

MTS_NAMESPACE_END

#endif /* __CJH_H */
