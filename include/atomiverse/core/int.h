/// @file core/int.h
/// @brief Compile-time integer precision, set by CMake -DATV_INT=int|int64_t.
#pragma once

#include <cstdint>

namespace ATV {

#ifndef ATV_INT_TYPE
/// @brief Preprocessor-selected integer type used by Atomiverse when CMake does
/// not override it.
#define ATV_INT_TYPE int
#endif

/// @brief Project-wide signed integer type.
///
/// Controlled at configure time through `-DATV_INT=int` or
/// `-DATV_INT=int64_t`.
using Int = ATV_INT_TYPE;

} // namespace ATV
