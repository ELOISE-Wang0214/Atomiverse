/// @file cpu/scf/internal/complex_ghf.h
/// @brief Complex two-component Hartree–Fock in the spin-orbital basis.
///
/// Transplants, contracted AO basis, `xuncontract=False`, `approx='1E'`:
///   pyscf.scf.ghf.get_hcore / get_ovlp / get_jk / get_occ / get_veff
///   pyscf.scf.hf.energy_elec / make_rdm1
///   pyscf.x2c.x2c.SpinOrbitalX2CHelper.get_hcore
///   pyscf.x2c.x2c._x2c1e_get_hcore / _sigma_dot / _block_diag
/// The modified-Dirac blocks inside `_x2c1e_get_hcore` are the
/// spin-orbital form of pyscf.scf.dhf.get_hcore (large–small kinetic
/// coupling, small-component `W/(4c²) − T`). Not transplanted:
/// pyscf.hessian.rhf.hess_elec (RHF CPHF; GHF.nuc_grad_method and
/// x2c.SCF.nuc_grad_method raise NotImplementedError) and
/// pyscf.df.df_jk.density_fit (JK selector; this path stores the dense ERI).
#pragma once

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/dense_matrix.h"
#include "atomiverse/cpu/scf/internal/basis_workspace.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace internal {

using ComplexMatrix = xt::xtensor<std::complex<double>, 2>;

/// @brief Cartesian `int1e_spnucsp` quaternion (PySCF component order).
///
/// `s` is ⟨p|V|p⟩ (`int1e_pnucp`, symmetric). `x,y,z` are the
/// antisymmetric cross terms ⟨p|V|p⟩_vector = p × (V p).
struct Spnucsp {
    RealMatrix x;
    RealMatrix y;
    RealMatrix z;
    RealMatrix s;
};

/// @brief ⟨σ·p | V_nuc | σ·p⟩ on real GTOs, four real matrices.
/// @param shells Cartesian shells (Bohr). l > 3 is rejected by the
///        nuclear primitive (`complex_ghf_ao_lgt3`).
/// @param atoms Nuclei (Bohr, charge Z).
Spnucsp int1e_spnucsp(const ShellList& shells, const AtomList& atoms);

/// @brief Spin-orbital embedding of a real AO matrix, `block_diag(A, A)`.
ComplexMatrix block_diag_real(const RealMatrix& a);

/// @brief `x2c._sigma_dot`: (i σ)·w_vec + w_scalar on the 2n spin-orbital
///        pairs. `wx,wy,wz,ws` are n×n.
ComplexMatrix sigma_dot(const RealMatrix& wx, const RealMatrix& wy,
                        const RealMatrix& wz, const RealMatrix& ws);

/// @brief Foldy–Wouthuysen X2C-1e core, `x2c._x2c1e_get_hcore`.
///
/// `t,v,w,s` are already spin-orbital (2 n_ao). `w` is `sigma_dot`
/// of `int1e_spnucsp`. Returns the 2 n_ao electronic hcore [Hartree].
/// @param c Speed of light [a.u.].
ComplexMatrix x2c1e_get_hcore(const ComplexMatrix& t, const ComplexMatrix& v,
                              const ComplexMatrix& w, const ComplexMatrix& s,
                              double c);

/// @brief Converged (or last) complex-GHF total energy.
struct ComplexGhfEnergy {
    double hartree = 0.0;
    bool converged = false;
    int iterations = 0;
};

/// @brief Complex GHF SCF. `x2c` selects `GHF.x2c1e` (`true`) or
///        non-relativistic `GHF.get_hcore` (`false`).
ComplexGhfEnergy complex_ghf_energy(const Frame& frame, BasisSet basis,
                                    int charge, double tol, int max_iter,
                                    bool x2c);

} // namespace internal
} // namespace scf
} // namespace cpu
} // namespace ATV
