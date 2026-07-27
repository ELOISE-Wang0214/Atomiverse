/// @file tests/cpu/scf/test_scf_cpmd.cpp
/// @brief esmm-cpmd-04-ao-cpmd — Car–Parrinello MD on the Gaussian SCF
///        methods: the refusal lattice, the dt-protocol seam, the fixed-state
///        gradient primitive, and the production-step signature.
///
/// Mirrors two production units:
///   - `include/atomiverse/cpu/scf/mean_field_method.h` — the `cpmd*` members,
///     the `TimestepConsumer` implementation, and the entry-point gates
///     (`MeanFieldMethodCpmd` suite);
///   - `include/atomiverse/cpu/scf/internal/gradient.h` —
///     `ScfGradient::compute_at_density` (`ScfGradientAtDensity` suite).
///
/// **Scope.** Unit only. The Car–Parrinello *physics* (conserved-energy drift,
/// adiabaticity, μ/dt tuning) is not testable in seconds and is pinned by
/// `regressions/dynamics/{cpmd-rhf-h2,cpmd-rks-water,cpmd-uhf-o2}` against
/// Car & Parrinello, PRL **55**, 2471 (1985) and Marx & Hutter, *Ab Initio
/// Molecular Dynamics* (CUP 2009) §2.4–2.5. What is here is the smallest
/// system that exercises each contract: H2/STO-3G (2 AOs, 2 electrons) and,
/// for the open-shell refusals, the H2⁺ doublet. Every refusal row throws
/// *before* any SCF work, so those tests cost a driver build.
///
/// **The propagator itself** (`electronic::Cpmd` — RATTLE-constrained
/// extended-Lagrangian advance) is method-agnostic dense algebra pinned by
/// tests/cpu/electronic/test_cpmd.cpp; nothing here re-tests it.
///
/// **Timestep.** Under `cpmd = true` the electronic propagation shares the
/// nuclear clock: dt arrives only through `cpu::TimestepConsumer::set_timestep`
/// (esmm-cpmd-02), pushed by `cpu::md::MD::build()`. There is no dt member to
/// set. The push walk's own rules live in
/// tests/cpu/dynamics/test_timestep_protocol.cpp; what is pinned here is the
/// SCF side of that seam.
///
/// **Throw-type convention** (mirrors the GPW CPMD twin,
/// src/cpu/gpw/gpw_evaluator.cpp, and is the contract these tests fix):
///   - a *static configuration combination*, knowable before any compute
///     (ROHF, RI, x2c, d3bj, …) → `std::invalid_argument`;
///   - *state discovered while driving* (no dt was ever pushed; near-singular
///     S at the seed) → `std::runtime_error`;
///   - the *wrong entry point* — standalone `run()` / `forces()` with
///     `cpmd = true`, i.e. CPMD outside an MD loop → `std::logic_error`
///     (spec Design; `std::invalid_argument` also satisfies it, being a
///     `std::logic_error`).
///
/// **Refusal lattice coverage** (ac-004). Rows pinned here: ROHF, ROKS, GHF,
/// `predictor()`, `ri != None`, `x2c`, `solvent_eps > 0`, `mom != Off`,
/// `soscf`, non-empty `point_charges`, nonzero `external_field`, `d3bj`, `d4`,
/// `n_dropped_ > 0` at seed, plus both standalone entries. The `xl_bomd` row is
/// pinned **from the XL side**: esmm-cpmd-05 added the member together with
/// `MeanFieldMethodXlBomd.RefusesCpmd` (tests/cpu/scf/test_scf_xl_bomd.cpp) and
/// the symmetric `validate_cpmd_config_` throw (`cpmd` + `xl_bomd` →
/// `std::invalid_argument`). Both flags live on one method object, so whichever
/// validator runs first sees the pair and one guard test covers both
/// directions; this file does not duplicate it.
///
/// **Determinism.** Fixed geometries and velocities, no RNG, no wall clock, no
/// filesystem, no network, no third-party oracle. The gradient parity bound is
/// an algebraic identity (below), so its 1e-10 tolerance is a hard-coded
/// property of the two formulas, not a fitted number.

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/dynamics/velocity_verlet.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/linalg/matrix.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/gradient_result.h"
#include "atomiverse/cpu/scf/internal/gradient.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/methods.h"
#include "atomiverse/cpu/scf/scf_enums.h"
#include "atomiverse/cpu/scf/slater_determinant.h"

#include "cpu/driver_test_helpers.h"

using namespace ATV;

using cpu::linalg::RealMatrix;
using cpu::linalg::RealVector;
using cpu::scf::BasisSet;
using cpu::scf::GradientResult;
using cpu::scf::IntegralPack;
using cpu::scf::SlaterDeterminant;
using Spin = cpu::scf::SlaterDeterminant::Spin;

namespace {

/// Nuclear MD timestep [fs] — the CPMD-appropriate small step; published by
/// the integrator and pushed into the method by `MD::build()`.
constexpr double kDtFs = 0.1;

/// Fictitious electronic mass [a.u.] (Car–Parrinello μ). The AO residual is
/// stiffer than a plane-wave one; 400 is the member default and is only ever
/// *carried* here — no adiabaticity claim is made by a unit test.
constexpr double kMu = 400.0;

/// Element-wise force-parity bound [eV/Å] between `compute_at_density` and
/// `compute` at a converged determinant (ac-001). This is not a fitted
/// tolerance: with C the canonical MOs of F, CᵀSC = I, the multiplier weight
/// W = occ·C_occ (C_occᵀ F C_occ) C_occᵀ equals occ·Σ_i ε_i c_i c_iᵀ *exactly*
/// (C_occᵀ F C_occ = diag(ε_occ)), so the two assemblies differ only by
/// floating-point summation order.
constexpr double kParityTolEvPerA = 1e-10;

bool message_mentions(const std::exception& e, const std::string& needle) {
    return std::string(e.what()).find(needle) != std::string::npos;
}

/// @brief H2 at the STO-3G-relevant equilibrium separation, 100 Å vacuum box.
/// @param bond H–H separation [Å].
/// @return Two-atom Frame (Z = 1, 1; masses in amu; atoms at rest).
cpu::Frame make_h2(double bond = 0.74) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_x)[1] = bond;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    f.span(kAtom_mass)[0] = 1.008;
    f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_type)[1] = 1;
    f.box = cpu::Box::cubic(100.0);
    return f;
}

/// @brief H2 with the two nuclei approaching along x at ±0.01 Å/fs.
///
/// A CPMD production step propagates C under the force of the *moved* nuclei;
/// with the nuclei frozen at the seed geometry the electronic force is the
/// (vanishing) SCF residual and T_e would be indistinguishable from zero.
/// @return Two-atom Frame with a real nuclear velocity [Å/fs].
cpu::Frame make_moving_h2() {
    cpu::Frame f = make_h2();
    f.span(kAtom_vx)[0] = 0.01;
    f.span(kAtom_vx)[1] = -0.01;
    return f;
}

/// @brief HeH⁺/STO-3G (Z = 2, 1) at r = 0.93 Å, nuclei approaching along x at
///        the same ±0.01 Å/fs as `make_moving_h2()` — used only by
///        `PublishesElectronicEnergyMetaKeys`.
///
/// H2 is homonuclear, so its AO Fock is persymmetric at every bond length:
/// the occupied sigma_g vector is an *exact* S-metric eigenvector throughout
/// the whole moving-nuclei trajectory, the CP band force
/// -(I - S C Cᵀ)(2 occ H C) is identically zero, and T_e stays exactly 0
/// under any correct implementation — not a bug, a symmetry identity.
/// HeH⁺ (Z = 2 ≠ 1) has no such symmetry, so the carried-forward
/// coefficients are only approximately the new eigenvector once the
/// geometry moves and T_e is genuinely positive, which is what this test
/// needs to observe. Charge is a property of the *method*
/// (`RHF::charge`, set by the caller), not the Frame.
/// @return Two-atom Frame (Z = 2, 1; masses in amu) with a real nuclear
///         velocity [Å/fs].
cpu::Frame make_moving_heh_plus() {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_x)[1] = 0.93;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    f.span(kAtom_mass)[0] = 4.0026;
    f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_type)[0] = 2;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_vx)[0] = 0.01;
    f.span(kAtom_vx)[1] = -0.01;
    f.box = cpu::Box::cubic(100.0);
    return f;
}

/// @brief RHF/STO-3G with CPMD switched on and nothing else unusual.
/// @return The method by value (moved into the driver by the caller).
cpu::scf::RHF cpmd_rhf() {
    cpu::scf::RHF rhf;
    rhf.basis = BasisSet::STO3G;
    rhf.tol = 1e-10;
    rhf.cpmd = true;
    rhf.cpmd_mu = kMu;
    return rhf;
}

/// @brief Bind `method` at the MD force phase and take exactly one step.
///
/// The uniform driver for every refusal row: `cpmd` validation fires at the
/// nested entry (or at the nested `build_on`), and both sites are inside this
/// call, so a row's `EXPECT_THROW` does not depend on which of the two the
/// implementation picks. `VelocityVerlet` is the `TimestepProvider`, so the dt
/// gate is satisfied and each row fails on its own gate.
/// @param method The concrete mean-field method under test (moved in).
/// @param frame  Geometry driven for one step.
template <typename Method>
void drive_nested_md_step(Method method, cpu::Frame& frame) {
    cpu::md::MD sim;
    sim.set_integrator(cpu::md::VelocityVerlet{kDtFs});
    sim.on(md::Phase::force).add(std::move(method));
    sim.build();
    sim.run(frame, 1);
}

cpu::scf::internal::ScfGradient make_grad(cpu::scf::SpinType spin) {
    return cpu::scf::internal::ScfGradient{cpu::scf::MeanFieldType::HF,
                                           spin,
                                           cpu::scf::XCType::LSDA,
                                           cpu::scf::GridLevel::SG0,
                                           0,
                                           0,
                                           false};
}

/// @brief Dense C = A·B for the tiny (2×2) AO matrices in this file.
template <typename A, typename B>
RealMatrix matmul(const A& a, const B& b) {
    const std::size_t n = a.shape()[0];
    const std::size_t k = a.shape()[1];
    const std::size_t m = b.shape()[1];
    RealMatrix out = xt::zeros<double>({n, m});
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t p = 0; p < k; ++p)
            for (std::size_t j = 0; j < m; ++j)
                out(i, j) += a(i, p) * b(p, j);
    return out;
}

RealMatrix transpose_of(const RealMatrix& a) {
    RealMatrix out = xt::zeros<double>({a.shape()[1], a.shape()[0]});
    for (std::size_t i = 0; i < a.shape()[0]; ++i)
        for (std::size_t j = 0; j < a.shape()[1]; ++j)
            out(j, i) = a(i, j);
    return out;
}

/// @brief The AO Fock matrix of a converged canonical determinant.
///
/// F C = S C ε with CᵀSC = I ⇒ C⁻¹ = CᵀS ⇒ **F = S C diag(ε) Cᵀ S**. Exact,
/// and it needs no Fock rebuild — the method's own `F_alpha_` is private
/// workspace, and reconstructing it from the published (C, ε, S) is what makes
/// the parity check a statement about the *formulas* rather than about the
/// solver's internals.
/// @param S   AO overlap (dimensionless).
/// @param C   MO coefficients of the channel (square, no dropped directions).
/// @param eps Orbital energies [Hartree].
/// @return F^σ in the AO basis [Hartree].
RealMatrix canonical_fock(const xt::xtensor<double, 2>& S, const RealMatrix& C,
                          const RealVector& eps) {
    RealMatrix C_eps = C;
    for (std::size_t i = 0; i < C_eps.shape()[0]; ++i)
        for (std::size_t j = 0; j < C_eps.shape()[1]; ++j)
            C_eps(i, j) *= eps(j);
    return matmul(matmul(S, C_eps), matmul(transpose_of(C), S));
}

/// @brief The leading `n_occ` (occupied) columns of a coefficient matrix.
RealMatrix occupied_columns(const RealMatrix& C, int n_occ) {
    RealMatrix out =
        xt::zeros<double>({C.shape()[0], static_cast<std::size_t>(n_occ)});
    for (std::size_t i = 0; i < C.shape()[0]; ++i)
        for (int j = 0; j < n_occ; ++j)
            out(i, static_cast<std::size_t>(j)) =
                C(i, static_cast<std::size_t>(j));
    return out;
}

void expect_forces_equal(const GradientResult& got, const GradientResult& want,
                         double tol) {
    ASSERT_EQ(got.forces.shape()[0], want.forces.shape()[0]);
    ASSERT_EQ(got.forces.shape()[1], want.forces.shape()[1]);
    for (std::size_t a = 0; a < want.forces.shape()[0]; ++a)
        for (std::size_t d = 0; d < want.forces.shape()[1]; ++d)
            EXPECT_NEAR(got.forces(a, d), want.forces(a, d), tol)
                << "atom " << a << " axis " << d;
}

} // namespace

// ─── the CPMD member surface (mirrors GPW, minus dt) ───────────────────

/// Parameters are **members with defaults**, not a `CpmdConfig` bag and not
/// constructor arguments — a BOMD user must be able to ignore all of them, and
/// a CPMD user must get the same μ / thermostat defaults the GPW evaluator
/// ships. `cpmd_dt` is deliberately absent: the clock comes from the
/// integrator (ac-005). The readbacks start at zero, before any electron has
/// moved.
TEST(MeanFieldMethodCpmd, MemberSurfaceDefaults) {
    cpu::scf::RHF rhf;
    EXPECT_FALSE(rhf.cpmd);
    EXPECT_DOUBLE_EQ(rhf.cpmd_mu, 400.0); // [a.u.] fictitious mass
    EXPECT_FALSE(rhf.cpmd_thermostat);
    EXPECT_DOUBLE_EQ(rhf.cpmd_ekin_target_ha, 1e-5); // [Hartree]
    EXPECT_DOUBLE_EQ(rhf.cpmd_tau_e, 50.0);          // [fs] coupling time
    EXPECT_DOUBLE_EQ(rhf.last_electronic_ke_ev(), 0.0);
    EXPECT_DOUBLE_EQ(rhf.last_electronic_thermostat_ev(), 0.0);
}

// ─── ac-004: refusal lattice — spin shapes that are not productized ─────

/// ROHF shares one set of spatial orbitals between the channels, so there are
/// no independent α/β coefficient matrices to propagate — the extended
/// Lagrangian is a different constrained problem, not a wiring detail.
TEST(MeanFieldMethodCpmd, RefusesRohf) {
    auto f = make_h2();
    cpu::scf::ROHF rohf{/*n_alpha=*/1, /*n_beta=*/0};
    rohf.charge = 1; // H2⁺ doublet — the natural restricted-open case
    rohf.basis = BasisSet::STO3G;
    rohf.cpmd = true;
    EXPECT_THROW(drive_nested_md_step(std::move(rohf), f),
                 std::invalid_argument);
}

TEST(MeanFieldMethodCpmd, RefusesRoks) {
    auto f = make_h2();
    cpu::scf::ROKS roks{/*n_alpha=*/1, /*n_beta=*/0};
    roks.charge = 1;
    roks.basis = BasisSet::STO3G;
    roks.cpmd = true;
    EXPECT_THROW(drive_nested_md_step(std::move(roks), f),
                 std::invalid_argument);
}

/// GHF propagates two-component spinors: Car–Parrinello on a spinor manifold
/// is new physics (the RATTLE metric is the spinor overlap), not RHF with
/// twice the columns.
TEST(MeanFieldMethodCpmd, RefusesGhf) {
    auto f = make_h2();
    cpu::scf::GHF ghf{/*n_electrons=*/1};
    ghf.charge = 1;
    ghf.basis = BasisSet::STO3G;
    ghf.cpmd = true;
    EXPECT_THROW(drive_nested_md_step(std::move(ghf), f),
                 std::invalid_argument);
}

// ─── ac-004: refusal lattice — competing density propagators ────────────

/// ASPC extrapolates the density across nuclear steps; CPMD propagates the
/// coefficients. Running both means two answers for "what is P now" — GPW
/// refuses the same pair.
TEST(MeanFieldMethodCpmd, RefusesPredictor) {
    auto f = make_h2();
    auto m = cpmd_rhf();
    m.predictor(3);
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::invalid_argument);
}

// ─── ac-004: refusal lattice — Fock builds CPMD's gradient does not cover ──

TEST(MeanFieldMethodCpmd, RefusesDensityFitting) {
    auto f = make_h2();
    auto m = cpmd_rhf();
    m.ri = cpu::scf::RI::JK;
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::invalid_argument);
}

TEST(MeanFieldMethodCpmd, RefusesX2c) {
    auto f = make_h2();
    auto m = cpmd_rhf();
    m.x2c = true;
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::invalid_argument);
}

TEST(MeanFieldMethodCpmd, RefusesSolvent) {
    auto f = make_h2();
    auto m = cpmd_rhf();
    m.solvent_eps = 78.4; // water
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::invalid_argument);
}

/// Maximum-overlap occupation re-selects the occupied subspace every
/// iteration; CPMD never diagonalizes, so there is no occupation to re-select.
TEST(MeanFieldMethodCpmd, RefusesMom) {
    auto f = make_h2();
    auto m = cpmd_rhf();
    m.mom = cpu::scf::MomStyle::IMOM;
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::invalid_argument);
}

/// SOSCF is a second-order *convergence* accelerator for the SCF loop the
/// production CPMD step deliberately does not run.
TEST(MeanFieldMethodCpmd, RefusesSoscf) {
    auto f = make_h2();
    auto m = cpmd_rhf();
    m.soscf = true;
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::invalid_argument);
}

/// The QM/MM × CPMD seam. Lifting this one gate later is the whole embedding
/// integration, so it must be shut loudly rather than silently half-working.
TEST(MeanFieldMethodCpmd, RefusesPointCharges) {
    auto f = make_h2();
    auto m = cpmd_rhf();
    m.point_charges.push_back({/*x=*/3.0, /*y=*/0.0, /*z=*/0.0, /*q=*/-0.5});
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::invalid_argument);
}

TEST(MeanFieldMethodCpmd, RefusesExternalField) {
    auto f = make_h2();
    auto m = cpmd_rhf();
    m.external_field = {{1e-4, 0.0, 0.0}}; // [a.u.]
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::invalid_argument);
}

// ─── ac-004: refusal lattice — energy terms with no force counterpart ───
//
// `ScfGradient` carries no dispersion force terms while `e_total_` includes
// E_disp (verified 2026-07-27, spec task 1). Under BOMD that is a
// pre-existing non-conservative force debt; under CPMD it would silently
// destroy the E_cons gate that makes the method verifiable, so both flags are
// refused outright.

TEST(MeanFieldMethodCpmd, RefusesD3Bj) {
    auto f = make_h2();
    auto m = cpmd_rhf();
    m.d3bj = true;
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::invalid_argument);
}

TEST(MeanFieldMethodCpmd, RefusesD4) {
    auto f = make_h2();
    auto m = cpmd_rhf();
    m.d4 = true;
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::invalid_argument);
}

// ─── ac-004: refusal lattice — ill-conditioned RATTLE metric at the seed ──

/// A near-singular overlap makes the holonomic constraint CᵀSC = I
/// ill-conditioned, so the propagated coefficients drift off the manifold.
/// `lindep_threshold = 0.5` drops the antibonding direction of H2/STO-3G
/// (S eigenvalues 1 ± S₁₂ ≈ {1.66, 0.34}) — a cheap, deterministic way to
/// reach `n_dropped_ > 0` without a diffuse basis. Discovered at the seed,
/// hence `std::runtime_error` (see the throw-type convention above).
TEST(MeanFieldMethodCpmd, RefusesDroppedOverlapDirectionsAtSeed) {
    auto f = make_h2();
    auto m = cpmd_rhf();
    m.lindep_threshold = 0.5;
    EXPECT_THROW(drive_nested_md_step(std::move(m), f), std::runtime_error);
}

// ─── ac-004: CPMD outside an MD loop is a category error, not a config one ──

/// `run()` is the standalone single-point entry: there is no nuclear
/// trajectory for the electrons to follow, so a CP propagation is meaningless
/// rather than merely unsupported. The dt is pushed first so this fails on the
/// entry-point gate and not on the missing-clock gate.
TEST(MeanFieldMethodCpmd, StandaloneRunRefusesCpmd) {
    auto f = make_h2();
    auto m = cpmd_rhf();
    m.set_timestep(kDtFs);
    EXPECT_THROW(m.run(f), std::logic_error);
}

TEST(MeanFieldMethodCpmd, StandaloneForcesRefusesCpmd) {
    auto f = make_h2();
    auto m = cpmd_rhf();
    m.set_timestep(kDtFs);
    EXPECT_THROW(m.forces(f), std::logic_error);
}

// ─── ac-005: dt arrives only through the TimestepConsumer protocol ──────

/// `timestep_required()` is the gate MD's push walk reads: it must track
/// `cpmd` exactly, so a plain BOMD method never demands an integrator and a
/// CPMD one always does.
TEST(MeanFieldMethodCpmd, TimestepRequiredTracksCpmdFlag) {
    cpu::scf::RHF rhf;
    EXPECT_FALSE(rhf.timestep_required());
    rhf.cpmd = true;
    EXPECT_TRUE(rhf.timestep_required());
}

TEST(MeanFieldMethodCpmd, SetTimestepRejectsNonPositive) {
    auto m = cpmd_rhf();
    EXPECT_THROW(m.set_timestep(0.0), std::invalid_argument);
}

TEST(MeanFieldMethodCpmd, SetTimestepRejectsNonFinite) {
    auto m = cpmd_rhf();
    EXPECT_THROW(m.set_timestep(std::numeric_limits<double>::infinity()),
                 std::invalid_argument);
}

/// The MD-side half of the seam: a nested method driver is reached by the push
/// walk through the `DriverNode` cross-cast, so a CPMD method in a loop whose
/// integrator publishes no clock is refused at `MD::build()` — naming the
/// protocol, because "add an integrator that provides it" is the fix.
TEST(MeanFieldMethodCpmd, NestedCpmdRequiresTimestepProvider) {
    auto f = make_h2();
    atv_test::StubIntegrator integrator; // an Integrator, but NO provider
    auto rhf = cpmd_rhf();
    cpu::md::MD sim; // declared last: the driver destructs first

    sim.on(md::Phase::initial_integrate).add_ref(integrator);
    sim.on(md::Phase::force).add_ref(rhf);

    try {
        sim.build();
        FAIL() << "expected MD::build() to refuse a cpmd method with no "
                  "TimestepProvider";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "TimestepProvider")) << e.what();
    }
}

/// The method-side half: a loop with no integrator at all (Evaluate) pushes
/// nothing, so the missing clock has to surface from the method's own build —
/// otherwise `cpmd_dt_fs_ = 0` would silently propagate the electrons by zero.
TEST(MeanFieldMethodCpmd, NestedBuildWithoutPushedTimestepThrows) {
    auto f = make_h2();
    auto rhf = cpmd_rhf();
    cpu::Evaluate eval; // declared last: the driver destructs first

    eval.on(evaluate::Phase::evaluate).add_ref(rhf);

    try {
        eval.build();
        eval.run(f);
        FAIL() << "expected the cpmd method to refuse to build with no pushed "
                  "timestep";
    } catch (const std::runtime_error& e) {
        EXPECT_TRUE(message_mentions(e, "timestep")) << e.what();
    }
}

// ─── ac-002: the seed is one full SCF, production steps are not ─────────

/// The cold seed must be a converged wavefunction: CP propagates *from* a
/// stationary state, and starting off it injects fictitious kinetic energy
/// that never leaves.
TEST(MeanFieldMethodCpmd, SeedStepRunsAFullScf) {
    auto f = make_moving_h2();
    auto rhf = cpmd_rhf();
    cpu::md::MD sim;

    sim.set_integrator(cpu::md::VelocityVerlet{kDtFs});
    sim.on(md::Phase::force).add_ref(rhf);
    sim.build();
    sim.run(f, /*n_steps=*/1);

    EXPECT_GT(rhf.iterations(), 0) << "the seed step must iterate the SCF";
    EXPECT_TRUE(rhf.converged()) << "an unconverged seed must throw, not seed";
    EXPECT_LT(rhf.energy(), 0.0);
}

/// The defining CPMD signature: after the seed, every step is a fixed-density
/// evaluation plus one electronic velocity-Verlet advance — zero SCF
/// iterations, forever.
TEST(MeanFieldMethodCpmd, ProductionStepsRunZeroScfIterations) {
    auto f = make_moving_h2();
    auto rhf = cpmd_rhf();
    cpu::md::MD sim;

    sim.set_integrator(cpu::md::VelocityVerlet{kDtFs});
    sim.on(md::Phase::force).add_ref(rhf);
    sim.build();
    sim.run(f, /*n_steps=*/3); // 1 seed + 2 production steps

    EXPECT_EQ(rhf.iterations(), 0)
        << "a production CP step must be a fixed-density eval, not an SCF";
    // Non-vacuity: the run really did produce a wavefunction and an energy.
    EXPECT_TRUE(rhf.has_result());
    EXPECT_LT(rhf.energy(), 0.0);
    EXPECT_TRUE(std::isfinite(static_cast<double>(f.meta.get(kEnergy))));
}

// ─── ac-002 / observability: the electronic energies reach the Frame ────

/// `CpmdConserved` (esmm-cpmd-03) throws on an absent key rather than silently
/// dropping T_e from E_cons, so publication is a contract, not a convenience —
/// and the published value must be the same number the readback reports.
///
/// System: HeH⁺, not H2 — see `make_moving_heh_plus()`. Homonuclear H2's
/// persymmetric AO Fock keeps the occupied MO an exact S-metric eigenvector
/// at every bond length, so its CP band force is identically zero and T_e
/// stays exactly 0.0 under any correct implementation; `EXPECT_GT(te, 0.0)`
/// below would fail on physical grounds, not an implementation bug. HeH⁺
/// breaks that symmetry, so T_e is genuinely positive.
TEST(MeanFieldMethodCpmd, PublishesElectronicEnergyMetaKeys) {
    auto f = make_moving_heh_plus();
    auto rhf = cpmd_rhf();
    rhf.charge = 1; // He(Z=2) + H(Z=1) - 1 = 2 electrons, closed-shell RHF
    cpu::md::MD sim;

    sim.set_integrator(cpu::md::VelocityVerlet{kDtFs});
    sim.on(md::Phase::force).add_ref(rhf);
    sim.build();

    // Open-schema keys: `Frame(int)` does not pre-insert them and nothing has
    // published them yet, so "present afterwards" is a real observation.
    ASSERT_FALSE(f.meta.contains(kEnergyElectronKinetic.key));
    ASSERT_FALSE(f.meta.contains(kEnergyElectronThermostat.key));

    sim.run(f, /*n_steps=*/3); // 1 seed + 2 production steps

    ASSERT_TRUE(f.meta.contains(kEnergyElectronKinetic.key));
    ASSERT_TRUE(f.meta.contains(kEnergyElectronThermostat.key));
    const double te = f.meta.get(kEnergyElectronKinetic);
    const double etherm = f.meta.get(kEnergyElectronThermostat);
    EXPECT_TRUE(std::isfinite(te)) << te;
    EXPECT_TRUE(std::isfinite(etherm)) << etherm;
    EXPECT_EQ(te, rhf.last_electronic_ke_ev());
    EXPECT_EQ(etherm, rhf.last_electronic_thermostat_ev());
    EXPECT_GT(te, 0.0) << "the electrons have followed two nuclear steps; a "
                          "T_e of exactly zero would make the equality vacuous";
    EXPECT_EQ(etherm, 0.0) << "electrons are NVE here (no cpmd_thermostat)";
}

// ─── ac-001: fixed-state gradient parity (the physics core) ─────────────

/// At a converged determinant the multiplier form of the Pulay weight,
/// W = occ·C_occ (C_occᵀ F C_occ) C_occᵀ, reduces to the energy-weighted
/// density occ·Σ_i ε_i c_i c_iᵀ that `compute` uses — because C_occᵀ F C_occ
/// = diag(ε_occ) exactly when C diagonalizes F in the S metric. That identity
/// is the ground truth for the non-eigenstate formula CPMD needs, so the two
/// must agree to floating-point noise (Pulay, Mol. Phys. **17**, 197 (1969);
/// Marx & Hutter §2.5).
///
/// Restricted leg: RHF/STO-3G H2, occ = 2 on the single α channel; the β
/// arguments are the empty matrices (a restricted method has one channel).
TEST(ScfGradientAtDensity, RestrictedMatchesComputeAtConvergedDeterminant) {
    auto f = make_h2();
    cpu::scf::RHF rhf;
    rhf.basis = BasisSet::STO3G;
    rhf.tol = 1e-10;
    rhf.max_iter = 200;
    rhf.run(f);
    ASSERT_TRUE(rhf.converged());

    const IntegralPack ints{f, BasisSet::STO3G, /*charge=*/0,
                            /*compute_eri=*/false};
    const SlaterDeterminant& det = rhf.determinant();
    ASSERT_TRUE(det.restricted());

    const auto grad = make_grad(cpu::scf::SpinType::Restricted);
    const GradientResult reference = grad.compute(f, ints, det, rhf.energy());

    const RealMatrix& C = det.coefficients(Spin::alpha);
    const RealMatrix F =
        canonical_fock(ints.S, C, det.orbital_energies(Spin::alpha));
    const RealMatrix C_occ = occupied_columns(C, det.n_electrons(Spin::alpha));

    const GradientResult at_density = grad.compute_at_density(
        f, ints, C_occ, F, RealMatrix{}, RealMatrix{}, rhf.energy());

    expect_forces_equal(at_density, reference, kParityTolEvPerA);
    EXPECT_NEAR(at_density.energy, reference.energy, kParityTolEvPerA);
}

/// Unrestricted leg: the β channel must carry its own multiplier weight with
/// occ = 1, not a doubled α. UHF H2 (n_α = n_β = 1) is the cheapest
/// two-channel determinant that converges from the shipped guess.
TEST(ScfGradientAtDensity, UnrestrictedMatchesComputeAtConvergedDeterminant) {
    auto f = make_h2();
    cpu::scf::UHF uhf{/*n_alpha=*/1, /*n_beta=*/1};
    uhf.basis = BasisSet::STO3G;
    uhf.tol = 1e-10;
    uhf.max_iter = 200;
    uhf.run(f);
    ASSERT_TRUE(uhf.converged());

    const IntegralPack ints{f, BasisSet::STO3G, /*charge=*/0,
                            /*compute_eri=*/false};
    const SlaterDeterminant& det = uhf.determinant();
    ASSERT_FALSE(det.restricted());

    const auto grad = make_grad(cpu::scf::SpinType::Unrestricted);
    const GradientResult reference = grad.compute(f, ints, det, uhf.energy());

    const RealMatrix& C_a = det.coefficients(Spin::alpha);
    const RealMatrix& C_b = det.coefficients(Spin::beta);
    const RealMatrix F_a =
        canonical_fock(ints.S, C_a, det.orbital_energies(Spin::alpha));
    const RealMatrix F_b =
        canonical_fock(ints.S, C_b, det.orbital_energies(Spin::beta));
    const RealMatrix C_occ_a =
        occupied_columns(C_a, det.n_electrons(Spin::alpha));
    const RealMatrix C_occ_b =
        occupied_columns(C_b, det.n_electrons(Spin::beta));

    const GradientResult at_density = grad.compute_at_density(
        f, ints, C_occ_a, F_a, C_occ_b, F_b, uhf.energy());

    expect_forces_equal(at_density, reference, kParityTolEvPerA);
    EXPECT_NEAR(at_density.energy, reference.energy, kParityTolEvPerA);
}
