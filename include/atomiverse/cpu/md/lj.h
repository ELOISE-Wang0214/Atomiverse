/// @file cpu/md/lj.h
/// @brief Backwards-compatibility alias header.
///
/// `ATV::cpu::md::LJ` is now an alias for `ATV::cpu::PairLJCut` after the
/// 2026-05-19 LAMMPS-naming refactor (`pe-potential-abstraction.md`).
/// New code should use `cpu::PairLJCut` from
/// `atomiverse/cpu/potentials/pair_lj_cut.h` directly. This header is
/// retained so existing call sites (`sim.set_pair(cpu::md::LJ{...})`)
/// compile unchanged.
#pragma once

#include "atomiverse/cpu/potentials/pair_lj_cut.h"

namespace ATV {
namespace cpu {
namespace md {

using LJ = ATV::cpu::PairLJCut;

} // namespace md
} // namespace cpu
} // namespace ATV
