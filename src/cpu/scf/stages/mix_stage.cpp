/// @file cpu/scf/stages/mix_stage.cpp
/// @brief MixStage::launch — DIIS error + Fock-space extrapolation + level
///        shift, byte-for-byte the shipping SCF loop body (steps 2 + 2b).

#include "atomiverse/cpu/scf/stages/mix_stage.h"

#include <cmath>

#include "atomiverse/cpu/scf/internal/level_shift.h"

namespace ATV {
namespace cpu {
namespace scf {

void MixStage::launch(const cpu::LaunchContext& ctx) {
    if (!ctx.at(ATV::scf::Phase::mix))
        return;

    // Step 2: DIIS error from the UNSHIFTED Fock (scf_solver.cpp:246-249).
    err_alpha_ = internal::diis_error(F_alpha_, P_alpha_, S_, X_);
    if (!restricted_)
        err_beta_ = internal::diis_error(F_beta_, P_beta_, S_, X_);

    // Step 2: DIIS extrapolation in Fock space (scf_solver.cpp:251-263).
    if (use_diis_ && iter_ >= diis_start_iter_) {
        if (restricted_) {
            F_alpha_ = diis_.step(F_alpha_, err_alpha_);
            F_beta_ = F_alpha_;
        } else {
            auto [Fa_diis, Fb_diis] =
                diis_.step(F_alpha_, err_alpha_, F_beta_, err_beta_);
            F_alpha_ = Fa_diis;
            F_beta_ = Fb_diis;
        }
    }

    // Step 2b: Saunders-Hillier level shift, strictly after DIIS (which built
    // its error and extrapolated the unshifted Fock) and strictly before
    // diagonalization (scf_solver.cpp:270-282). Restricted D_occ = 0.5*P
    // (occupancy-1 convention); unrestricted P_sigma is already S-idempotent.
    const bool shift_active =
        std::abs(level_shift_) > internal::kLevelShiftSkipThreshold;
    if (shift_active) {
        if (restricted_) {
            RealMatrix D_occ = 0.5 * P_alpha_;
            F_alpha_ =
                internal::apply_level_shift(F_alpha_, S_, D_occ, level_shift_);
            F_beta_ = F_alpha_;
        } else {
            F_alpha_ = internal::apply_level_shift(F_alpha_, S_, P_alpha_,
                                                   level_shift_);
            F_beta_ =
                internal::apply_level_shift(F_beta_, S_, P_beta_, level_shift_);
        }
    }
}

} // namespace scf
} // namespace cpu
} // namespace ATV
