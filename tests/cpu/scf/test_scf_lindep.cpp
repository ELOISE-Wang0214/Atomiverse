/// @file tests/cpu/scf/test_scf_lindep.cpp
/// @brief SCF linear-dependence handling via canonical
///        orthogonalization (hf-robust-02, ac-002/ac-003/ac-004).
///
/// Two systems:
///  - clean H2 / 6-31G (no linear dependence): the default
///    SCFConfig::lindep_threshold = 1e-8 must leave the dsygvd fast path
///    untouched — energy pinned bit-identical to the pre-change baseline;
///    lindep_threshold = 0 (disabled) must give the same bits.
///  - H2 + H+ / 6-31G with a third proton at (0, 0, 1e-5 Å), nearly
///    coincident with the first H: S eigenvalues {3.5e-12, 1.1e-10, 0.126,
///    ...} force the canonical-orthogonalization path, dropping exactly 2
///    columns (n_mo = 4 of n_ao = 6).
///
/// PySCF 2.13.1 reference (cart=True, conv_tol=1e-12), regenerated for this
/// exact geometry/charge on 2026-07-02 with a custom eig hook implementing
/// raw-S canonical orthogonalization at tau = 1e-8
/// (X = U_kept / sqrt(s_kept)):
///   S eigenvalues: 3.531647e-12 1.093954e-10 1.259887e-01 6.175935e-01
///                  9.063536e-01 4.350064e+00
///   n_dropped = 2, n_mo = 4
///   E_total      =  52914.892606617992 Ha  (V_nn = 52919.148609614618 Ha)
///   E_electronic =     -4.256002996626 Ha
/// The total energy is dominated by the nuclear repulsion of the two
/// near-coincident protons (~5.3e4 Ha), which is sensitive to the Å→Bohr
/// constant at the 1e-6-Ha scale, so the cross-code comparison pins the
/// ELECTRONIC energy E - V_nn instead.
///
/// Tutorial regression: the ch02 section "Why HF becomes a matrix equation —
/// and when the basis goes bad" quotes s_min ≈ 3.5e-12, 2 dropped columns,
/// and the clean-H2 energy below; those literals live HERE and in the
/// tutorial in lockstep.

#ifdef ATV_HAS_RHF

#include <cmath>

#include <gtest/gtest.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/scf/scf_solver.h"

using namespace ATV;

namespace {

// Pre-change baseline: single-thread probe at base 2c98906, clean H2/6-31G,
// rhf_config, 17 significant digits (ac-002; dsygvd fast path).
constexpr double kH2Rhf631gPreChange = -1.1267339671159133;

// PySCF 2.13.1 electronic energy for the near-degenerate system (see @file).
constexpr double kPyscfElectronicRef = -4.256002996626;

// PySCF 2.13.1 smallest raw S eigenvalue (dimensionless; tutorial ch02
// quotes 3.5e-12).
constexpr double kPyscfSmin = 3.531647e-12;

/// Clean H2 at R = 0.7414 Å.
cpu::Frame make_h2() {
    cpu::Frame f(2);
    f.atoms.x = {0.0, 0.0};
    f.atoms.y = {0.0, 0.0};
    f.atoms.z = {0.0, 0.7414};
    f.atoms.type = {1, 1};
    f.atoms.mass = {1.008, 1.008};
    f.box = Mat3::cubic(50.0);
    return f;
}

/// H2 + H+: third proton nearly coincident with the first center
/// (0, 0, 1e-5 Å). charge = +1 keeps 2 electrons.
cpu::Frame make_h2_hplus() {
    cpu::Frame f(3);
    f.atoms.x = {0.0, 0.0, 0.0};
    f.atoms.y = {0.0, 0.0, 0.0};
    f.atoms.z = {0.0, 0.7414, 1e-5};
    f.atoms.type = {1, 1, 1};
    f.atoms.mass = {1.008, 1.008, 1.008};
    f.box = Mat3::cubic(50.0);
    return f;
}

} // namespace

// ---------------------------------------------------------------------------
// ac-002 — healthy-basis no-regression (dsygvd fast path untouched).
// ---------------------------------------------------------------------------

// tutorial ch02 "when the basis goes bad" worked example (clean-H2 energy);
// regression anchor test_cpu_scf_lindep (ac-002, ac-004). Single OMP thread
// pins the reduction order so the bit-identity comparisons are well-defined.
TEST(SCFLindep, CleanH2DefaultThresholdBitIdenticalToBaseline) {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    auto frame = make_h2();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631G);
    auto cfg = cpu::scf::rhf_config(ints); // lindep_threshold default 1e-8

    const auto result = cpu::scf::scf_solve(ints, cfg);
    ASSERT_TRUE(result.converged);
    // Healthy basis => n_dropped == 0 => dsygvd path — value must not move
    // at all (1e-12 spec bar; asserted bit-identical, which is stronger).
    EXPECT_DOUBLE_EQ(result.energy, kH2Rhf631gPreChange);
    // No columns dropped: full n_mo == n_ao.
    EXPECT_EQ(result.eigenvalues_alpha.size(), 4u);
}

// ac-002: lindep_threshold = 0 disables filtering entirely and must
// reproduce the current dsygvd result bitwise.
TEST(SCFLindep, CleanH2ThresholdZeroBitwiseIdentical) {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    auto frame = make_h2();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631G);
    auto cfg = cpu::scf::rhf_config(ints);
    cfg.lindep_threshold = 0.0;

    const auto result = cpu::scf::scf_solve(ints, cfg);
    ASSERT_TRUE(result.converged);
    EXPECT_DOUBLE_EQ(result.energy, kH2Rhf631gPreChange);
}

// ---------------------------------------------------------------------------
// ac-003 — near-degenerate H2 + H+ converges on the filtered subspace.
// ---------------------------------------------------------------------------

// tutorial ch02 "when the basis goes bad" worked example (s_min, dropped
// columns); regression anchor test_cpu_scf_lindep (ac-003, ac-004).
TEST(SCFLindep, NearDegenerateOverlapSpectrumMatchesPyscf) {
    auto frame = make_h2_hplus();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631G,
                                          /*charge=*/+1);
    ASSERT_EQ(ints.n_electrons, 2);
    ASSERT_EQ(ints.n_basis, 6);

    cpu::linalg::CanonicalOrth ortho(ints.S, /*tau=*/1e-8);
    EXPECT_EQ(ortho.n_dropped, 2)
        << "two near-null overlap directions below tau = 1e-8";
    EXPECT_EQ(ortho.X.shape(0), 6u);
    EXPECT_EQ(ortho.X.shape(1), 4u);
    // Raw smallest eigenvalue ~ 3.5e-12 (PySCF cross-check, 10% slop for
    // integral last-digit differences between codes).
    EXPECT_NEAR(ortho.s_min, kPyscfSmin, 0.35e-12);
}

// ac-003: RHF on the near-degenerate system converges with exactly 2
// dropped columns (n_mo == 4 of n_ao == 6) and matches the regenerated
// PySCF electronic energy to < 1e-6 Ha.
TEST(SCFLindep, NearDegenerateRhfMatchesPyscfElectronicEnergy) {
    auto frame = make_h2_hplus();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631G,
                                          /*charge=*/+1);
    auto cfg = cpu::scf::rhf_config(ints); // lindep_threshold default 1e-8

    const auto result = cpu::scf::scf_solve(ints, cfg);
    ASSERT_TRUE(result.converged)
        << "canonical orthogonalization must rescue the near-singular basis";
    // Rectangular spectrum: 4 of 6 columns survive.
    EXPECT_EQ(result.eigenvalues_alpha.size(), 4u);
    EXPECT_EQ(result.C_alpha.shape(0), 6u);
    EXPECT_EQ(result.C_alpha.shape(1), 4u);
    // Density stays n x n.
    EXPECT_EQ(result.P_alpha.shape(0), 6u);
    EXPECT_EQ(result.P_alpha.shape(1), 6u);

    const double e_electronic = result.energy - ints.V_nuc;
    EXPECT_NEAR(e_electronic, kPyscfElectronicRef, 1e-6)
        << "electronic energy vs PySCF 2.13.1 (raw-S tau=1e-8 reference)";
}

#endif // ATV_HAS_RHF
