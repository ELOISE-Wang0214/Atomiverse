/// @file core/real.h
/// @brief Compile-time floating-point precision, set by CMake -DATV_REAL=float|double.
#pragma once

namespace ATV {

#ifndef ATV_REAL_TYPE
/// @brief Preprocessor-selected floating-point type used by Atomiverse when
/// CMake does not override it.
#define ATV_REAL_TYPE double
#endif

/// @brief Project-wide floating-point type.
///
/// Controlled at configure time through `-DATV_REAL=float` or
/// `-DATV_REAL=double`.
using Real = ATV_REAL_TYPE;

} // namespace ATV
