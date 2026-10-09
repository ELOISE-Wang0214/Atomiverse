/// @file cpu/scf/internal/dhf.h
/// @brief Four-component Dirac–Hartree–Fock driver (private).
///
/// Transplants `pyscf.scf.dhf` (get_hcore, get_ovlp, get_occ,
/// energy_elec, _call_veff_llll, _call_veff_ssll) on l ≤ 1 shells.
/// The j-adapted coefficients are `pyscf.symm.sph.sph2spinor` for l = 0, 1.
#pragma once

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace internal {

/// @brief Converged (or last) Dirac–Hartree–Fock total energy.
struct DhfEnergy {
    double hartree = 0.0;
    bool converged = false;
    int iterations = 0;
};

/// @brief Dirac–Coulomb SCF. Energy is Hartree.
/// @throws std::runtime_error naming `dhf_ao_lgt1` when a shell has l > 1.
DhfEnergy dhf_energy(const Frame& frame, BasisSet basis, int charge, double tol,
                     int max_iter);

} // namespace internal
} // namespace scf
} // namespace cpu
} // namespace ATV
