#if !defined(__NEURAL_UTIL_H)
#define __NEURAL_UTIL_H

#include <mitsuba/mitsuba.h>

#include <string>
#include <vector>

MTS_NAMESPACE_BEGIN

void imageFromVector(const fs::path &filename, std::vector<Float> values);

MTS_NAMESPACE_END

#endif /* __NEURAL_UTIL_H */
