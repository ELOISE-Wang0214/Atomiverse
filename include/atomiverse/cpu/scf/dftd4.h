/// @file cpu/scf/dftd4.h
/// @brief Public molecular DFT-D4 with Becke–Johnson damping (PBE defaults).
///
/// Geometry-only additive dispersion correction for the Gaussian-basis path
/// (Caldeweyher et al., J. Chem. Phys. 150, 154122 (2019)). Independent of the
/// electronic wavefunction: evaluate on a Frame of atomic numbers + positions
/// [Å]. PySCF twin: `disp='d4'` (default rational damping, two-body term;
/// three-body ATM is not included — for H₂O it is ≪ 1e-10 Ha).
///
/// @code
///   cpu::scf::D4 d4;                   // PBE-D4 defaults
///   double e_disp = d4.energy(frame);  // [Hartree]
///
///   // Or fold into RKS total energy:
///   cpu::scf::RKS rks(cpu::scf::XCType::PBE);
///   rks.d4 = true;
///   rks.run(frame);
///   double ed = rks.energy_disp();     // [Hartree]
/// @endcode
#pragma once

#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Molecular Grimme DFT-D4 with Becke–Johnson (rational) damping.
///
/// Default parameters are the PBE-D4 set of Caldeweyher et al. (SAW190103
/// BJ+ATM fit; two-body only here): s6 = 1, s8 = 0.95948085, a1 = 0.38574991,
/// a2 = 4.80688534 Bohr. Charge model is EEQ2019; ga = 3, gc = 2, wf = 6.
/// Tables are in-tree (Z ≤ 94). Three-body ATM is not included.
class D4 {
  public:
    double s6 = 1.0;          ///< C6 scale (dimensionless; PBE = 1).
    double s8 = 0.95948085;   ///< C8 scale (dimensionless; PBE-D4).
    double a1 = 0.38574991;   ///< BJ steepness (dimensionless; PBE).
    double a2 = 4.80688534;   ///< BJ offset [Bohr]; PBE.
    double ga = 3.0;          ///< Charge-scaling height.
    double gc = 2.0;          ///< Charge-scaling steepness.
    double wf = 6.0;          ///< CN Gaussian weight factor.
    double charge = 0.0;      ///< Molecular charge for EEQ (electrons).

    /// @brief Dispersion energy for the Frame geometry.
    ///
    /// Reads `kAtom_type` (atomic number Z) and `kAtom_x/y/z` [Å]. Isolated
    /// (non-periodic) evaluation — the Frame box is ignored. Unsupported Z
    /// (outside the in-tree D4 table, Z > 94 or no reference states) throws.
    ///
    /// @param frame System geometry (Å, Z).
    /// @return E_disp [Hartree] (negative when dispersion is attractive).
    double energy(const Frame& frame) const;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
