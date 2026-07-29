/// @file cpu/component/qm_engine.h
/// @brief Shared QM/MM embedding site types (protocol classes land in
///        esmm-qmmm-02).
#pragma once

namespace ATV {
namespace cpu {

/// @brief One external MM-like embedding charge (public currency).
///
/// Positions in Å, charge in e. Used by mean-field `point_charges` and by
/// the private embedding gradient path. Not a Config bag — a plain value type.
struct EmbeddingSite {
    double x = 0; ///< Position x [Å].
    double y = 0; ///< Position y [Å].
    double z = 0; ///< Position z [Å].
    double q = 0; ///< Charge [e].
};

/// @brief Force on one embedding site (public currency).
///
/// Components in eV/Å. Filled by the analytic embedding gradient when
/// `point_charges` is non-empty; empty otherwise.
struct EmbeddingForce {
    double fx = 0; ///< Force x [eV/Å].
    double fy = 0; ///< Force y [eV/Å].
    double fz = 0; ///< Force z [eV/Å].
};

} // namespace cpu
} // namespace ATV
