/// @file core/real.h
/// @brief Legacy project-wide real type. Prefer the role-specific tiers in
/// core/precision.h (::ATV::StorageReal / ::ATV::WorkReal / ::ATV::AccReal).
#pragma once

#include "atomiverse/core/precision.h"

namespace ATV {

/// @brief Project-wide floating-point type (legacy alias).
///
/// Kept for backward compatibility: it aliases ::ATV::AccReal (always
/// `double`), the accumulator-precision tier. New code should choose an
/// explicit tier from core/precision.h — ::ATV::StorageReal for per-atom
/// stored fields, ::ATV::WorkReal for transient kernel math, ::ATV::AccReal for
/// accumulators — instead of this catch-all. Precision is no longer a single
/// global `-DATV_REAL` switch; see core/precision.h for the tiered controls.
using Real = AccReal;

} // namespace ATV
