/// @file core/int.h
/// @brief Legacy project-wide integer type. Prefer the role-specific tiers in
/// core/precision.h (::ATV::LocalInt / ::ATV::TagInt / ::ATV::ImageInt /
/// ::ATV::BigInt).
#pragma once

#include "atomiverse/core/precision.h"

namespace ATV {

/// @brief Project-wide signed integer type (legacy alias).
///
/// Kept for backward compatibility: it aliases ::ATV::LocalInt (32-bit), the
/// on-rank count/index tier. New code should choose an explicit tier from
/// core/precision.h — ::ATV::LocalInt for local counts, ::ATV::TagInt for atom
/// IDs, ::ATV::ImageInt for image flags, ::ATV::BigInt for global totals — so
/// 64-bit width is spent only where a count can overflow. Precision is no
/// longer a single global `-DATV_INT` switch; see core/precision.h.
using Int = LocalInt;

} // namespace ATV
