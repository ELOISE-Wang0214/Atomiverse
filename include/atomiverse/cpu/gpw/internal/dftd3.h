/// @file cpu/gpw/internal/dftd3.h
/// @brief Grimme DFT-D3 dispersion correction (zero / "original" damping) for
///        the PBE reference functional, matching CP2K's `TYPE DFTD3`.
///
/// A self-contained additive correction (energy + analytic forces) on top of
/// the GPW DFT, evaluated as a pairwise sum over the (periodically imaged)
/// atoms. Method: Grimme, Antony, Ehrlich, Krieg, J. Chem. Phys. 132, 154104
/// (2010), zero damping; bit-for-bit faithful to CP2K's `qs_dispersion_d3`:
///
///   E_disp = - Σ_{A<B} Σ_{n=6,8} s_n (C_n^AB / R_AB^n) f_{d,n}(R_AB)
///   f_{d,n}(R) = 1 / (1 + 6 (R / (s_{r,n} R0_AB))^{-α_n}),  α6=14, α8=16
///   C_8^AB = 3 C_6^AB r2r4_A r2r4_B,  R0_AB from CP2K's r0ab table [Bohr]
///   C_6^AB = CN-interpolated reference C6 (k3 = -4)
///   CN_A   = Σ_{B≠A} 1/(1+exp(-k1((rcov_A+rcov_B)/R_AB - 1))) · ½(1-tanh(R-2 rcov_AB))
///
/// PBE parameters: s6 = 1.0, s8 = 0.722, sr6 = 1.217, sr8 = 1.0. A smooth
/// fifth-order cutoff (CP2K `damping_d3`) zeroes each term over [rcut-1, rcut].
/// Periodic (orthorhombic, Γ-point) lattice sums use the same image enumeration
/// as the GPW integrals; the three-body (C9 / ATM) term is not included (CP2K's
/// default with `TYPE DFTD3` for the ZnO_MD target leaves it off).
#pragma once

#include <array>
#include <string>
#include <vector>

#include "atomiverse/cpu/gpw/internal/dftd3_params.h"

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Result of a DFT-D3 evaluation.
struct D3Result {
    double energy = 0.0; ///< Dispersion energy [Hartree].
    /// @brief Forces F_A = -dE/dR_A [Hartree/Bohr], one xyz triple per atom.
    std::vector<std::array<double, 3>> forces;
};

/// @brief PBE-D3 (zero damping) dispersion energy and forces.
///
/// @param Z        Atomic numbers, length N (1-based element identity).
/// @param pos_bohr Atom positions [Bohr], length N.
/// @param cell_bohr Orthorhombic box lengths a,b,c [Bohr] (used iff periodic).
/// @param periodic When true, sum over lattice images within `rcut_bohr`; when
///        false, a single isolated cluster (no images).
/// @param params   D3 reference parameters (`D3Params{}`, in-tree tables).
/// @param rcut_bohr Pair cutoff [Bohr]; CP2K default is 2*R_CUTOFF = 40.
/// @return Energy [Hartree] and forces [Hartree/Bohr].
D3Result dftd3_pbe(const std::vector<int>& Z,
                   const std::vector<std::array<double, 3>>& pos_bohr,
                   const double cell_bohr[3], bool periodic,
                   const D3Params& params, double rcut_bohr = 40.0);

/// @brief PBE-D3 with Becke–Johnson (rational) damping — molecular/PySCF path.
///
///   E = −Σ_{A<B} [ s6 C6/(R⁶ + f⁶) + s8 C8/(R⁸ + f⁸) ]
///   f  = a1 √(C8/C6) + a2
///
/// PBE-D3(BJ) defaults (Grimme et al., J. Comput. Chem. 32, 1456 (2011)):
/// s6 = 1.0, s8 = 0.7875, a1 = 0.4289, a2 = 4.4407 Bohr. Coordination numbers
/// use the original Grimme form (no CP2K long-range tanh factor), matching
/// simple-dftd3 / PySCF `disp='d3bj'`. Three-body ATM is off.
///
/// @param Z         Atomic numbers, length N.
/// @param pos_bohr  Atom positions [Bohr], length N.
/// @param cell_bohr Orthorhombic box lengths a,b,c [Bohr] (used iff periodic).
/// @param periodic  Lattice images within `rcut_bohr` when true; isolated else.
/// @param params    D3 reference parameters (`D3Params{}`).
/// @param rcut_bohr Pair cutoff [Bohr] (default 94.2 ≈ 50 Å).
/// @param s6        C6 scale (default 1.0 for PBE).
/// @param s8        C8 scale (default 0.7875 for PBE-BJ).
/// @param a1        BJ steepness (default 0.4289 for PBE).
/// @param a2        BJ offset [Bohr] (default 4.4407 for PBE).
/// @return Energy [Hartree] and forces [Hartree/Bohr].
D3Result dftd3_pbe_bj(const std::vector<int>& Z,
                      const std::vector<std::array<double, 3>>& pos_bohr,
                      const double cell_bohr[3], bool periodic,
                      const D3Params& params, double rcut_bohr = 94.2,
                      double s6 = 1.0, double s8 = 0.7875, double a1 = 0.4289,
                      double a2 = 4.4407);

} // namespace gpw
} // namespace cpu
} // namespace ATV
