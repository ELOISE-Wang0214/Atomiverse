/// @file core/constants.h
/// @brief Deprecated — included only for back-compat. New code should
///        include `atomiverse/core/units.h` directly.
///
/// Historically this file held `kPI`, `kDegToRad`, `kAngstromToBohr`,
/// while `units.h` held the atomistic conversion factors. The split
/// was arbitrary — both kinds of constants described unit conversions.
/// `units.h` is now the single registry; this file just forwards.
#pragma once

#include "atomiverse/core/units.h"
