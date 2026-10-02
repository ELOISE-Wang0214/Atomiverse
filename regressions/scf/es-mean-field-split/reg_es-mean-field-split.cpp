/// @file regressions/scf/es-mean-field-split/reg_es-mean-field-split.cpp
/// @brief Mean-field split on classic systems, pinned to offline PySCF.
///
/// H₂/STO-3G: `RHF::run` and `RHF::forces` stay on the method. MP2 and the
/// CPHF nuclear Hessian do not (`MP2<RHF>`, `Hessian<RHF>`).
/// H₂O/6-31G: RHF energy stays on the method. C-PCM and static
/// polarizability attach beside it (`WithSolvent<RHF>`,
/// `Polarizability<RHF>`).
///
/// Out-of-library consumer: public headers + link atomiverse only.
/// Pins and gates: REFERENCE.md. Exit non-zero on mismatch. ≤ 60 s.

#include <cmath>
#include <cstdio>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"

namespace {

ATV::cpu::Frame h2() {
    ATV::cpu::Frame frame(2);
    frame.span(ATV::kAtom_type)[0] = 1;
    frame.span(ATV::kAtom_x)[0] = 0.0;
    frame.span(ATV::kAtom_y)[0] = 0.0;
    frame.span(ATV::kAtom_z)[0] = 0.0;
    frame.span(ATV::kAtom_type)[1] = 1;
    frame.span(ATV::kAtom_x)[1] = 0.0;
    frame.span(ATV::kAtom_y)[1] = 0.0;
    frame.span(ATV::kAtom_z)[1] = 0.74;
    frame.box = ATV::cpu::Box::cubic(20.0);
    return frame;
}

ATV::cpu::Frame h2o() {
    ATV::cpu::Frame frame(3);
    frame.span(ATV::kAtom_type)[0] = 8;
    frame.span(ATV::kAtom_x)[0] = 0.0;
    frame.span(ATV::kAtom_y)[0] = 0.0;
    frame.span(ATV::kAtom_z)[0] = 0.1173;
    frame.span(ATV::kAtom_type)[1] = 1;
    frame.span(ATV::kAtom_x)[1] = 0.0;
    frame.span(ATV::kAtom_y)[1] = 0.7572;
    frame.span(ATV::kAtom_z)[1] = -0.4692;
    frame.span(ATV::kAtom_type)[2] = 1;
    frame.span(ATV::kAtom_x)[2] = 0.0;
    frame.span(ATV::kAtom_y)[2] = -0.7572;
    frame.span(ATV::kAtom_z)[2] = -0.4692;
    frame.box = ATV::cpu::Box::cubic(20.0);
    return frame;
}

int fail(const char* msg) {
    std::printf("FAIL: %s\n", msg);
    return 1;
}

} // namespace

int main() {
    using namespace ATV;

    // ── H₂ / STO-3G: run and forces on RHF; MP2 and CPHF Hessian beside ──
    auto h2_frame = h2();
    cpu::scf::RHF rhf{cpu::scf::BasisSet::STO3G, /*charge=*/0, /*tol=*/1e-12};
    rhf.max_iter = 100;
    rhf.run(h2_frame);
    if (!rhf.converged())
        return fail("H2 RHF did not converge");
    constexpr double kH2Energy = -1.1167593073964255;
    if (!(std::abs(rhf.energy() - kH2Energy) <= 1e-8)) {
        std::printf("FAIL: H2 E_RHF got=%.12f pin=%.12f\n", rhf.energy(),
                    kH2Energy);
        return 1;
    }

    const cpu::scf::GradientResult g = rhf.forces(h2_frame);
    auto grad = [&](int atom, int axis) {
        return -units::convert(g.forces(atom, axis),
                               units::Force::eVPerAngstrom,
                               units::Force::HartreePerBohr);
    };
    const double want_g[2] = {-0.02767960070679465, 0.02767960070679487};
    double max_g = 0.0;
    for (int atom = 0; atom < 2; ++atom) {
        for (int axis = 0; axis < 3; ++axis) {
            const double got = grad(atom, axis);
            const double want = (axis == 2) ? want_g[atom] : 0.0;
            max_g = std::max(max_g, std::abs(got - want));
        }
    }
    std::printf("H2 |grad err| = %.3e Ha/Bohr\n", max_g);
    if (!(max_g <= 1e-6))
        return fail("H2 RHF gradient off PySCF");

    const cpu::scf::MP2Result mp2 =
        cpu::scf::MP2<cpu::scf::RHF>{rhf}.run(h2_frame);
    constexpr double kMp2Tot = -1.1298973809859585;
    constexpr double kMp2Corr = -0.013138073589533008;
    std::printf("H2 MP2 e_tot=%.12f e_corr=%.12f\n", mp2.e_tot, mp2.e_corr);
    if (!(std::abs(mp2.e_tot - kMp2Tot) <= 1e-6))
        return fail("H2 MP2 total off PySCF");
    if (!(std::abs(mp2.e_corr - kMp2Corr) <= 1e-7))
        return fail("H2 MP2 correlation off PySCF");

    const cpu::vib::HessianResult hess =
        cpu::scf::Hessian<cpu::scf::RHF>{rhf}.run(h2_frame);
    // Atom-major 6×6 from pyscf.hessian.RHF.kernel, eV/Å².
    constexpr double kHxx = 1.923435534077;
    constexpr double kHzz = 47.077486422795;
    const double want_h[6][6] = {
        {kHxx, 0, 0, -kHxx, 0, 0}, {0, kHxx, 0, 0, -kHxx, 0},
        {0, 0, kHzz, 0, 0, -kHzz}, {-kHxx, 0, 0, kHxx, 0, 0},
        {0, -kHxx, 0, 0, kHxx, 0}, {0, 0, -kHzz, 0, 0, kHzz},
    };
    if (hess.n_atoms != 2)
        return fail("H2 Hessian atom count");
    double max_h = 0.0;
    for (int r = 0; r < 6; ++r)
        for (int c = 0; c < 6; ++c)
            max_h =
                std::max(max_h, std::abs(hess.hessian(r, c) - want_h[r][c]));
    std::printf("H2 |Hessian err| = %.3e eV/Ang^2\n", max_h);
    if (!(max_h <= 1e-5))
        return fail("H2 CPHF Hessian off PySCF");

    // ── H₂O / 6-31G: energy on RHF; solvent and polarizability beside ──
    auto water = h2o();
    cpu::scf::RHF water_rhf{cpu::scf::BasisSet::Pople631G, 0, 1e-12};
    water_rhf.max_iter = 100;
    water_rhf.run(water);
    if (!water_rhf.converged())
        return fail("H2O RHF did not converge");
    constexpr double kWaterEnergy = -75.98397447272193;
    if (!(std::abs(water_rhf.energy() - kWaterEnergy) <= 1e-8)) {
        std::printf("FAIL: H2O E_RHF got=%.12f pin=%.12f\n", water_rhf.energy(),
                    kWaterEnergy);
        return 1;
    }

    const auto alpha = cpu::scf::Polarizability<cpu::scf::RHF>{water_rhf}.run();
    const double want_a[3] = {1.394213713033, 6.648784943839, 4.410031211301};
    for (int a = 0; a < 3; ++a) {
        const double err = std::abs(alpha(a, a) - want_a[a]);
        std::printf("alpha[%d%d] got=%.10f pin=%.10f |err|=%.3e\n", a, a,
                    alpha(a, a), want_a[a], err);
        if (!(err <= 1e-5))
            return fail("H2O polarizability diagonal off PySCF");
        for (int b = 0; b < a; ++b) {
            if (!(std::abs(alpha(a, b)) <= 1e-4))
                return fail("H2O polarizability off-diagonal");
        }
    }

    cpu::scf::RHF sol_mf{cpu::scf::BasisSet::Pople631G, 0, 1e-10};
    sol_mf.max_iter = 100;
    cpu::scf::ContinuumSolvent pcm;
    pcm.eps = 78.4;
    pcm.radii_scale = 1.2;
    pcm.x = 0.0;
    pcm.model = cpu::scf::SolventModel::CPCM;
    cpu::scf::WithSolvent<cpu::scf::RHF> sol{std::move(sol_mf), pcm};
    sol.run(water);
    if (!sol.mf.converged())
        return fail("H2O C-PCM did not converge");
    constexpr double kHaInKcal = 627.5094740631;
    constexpr double kPinSolvHa = -0.01614471877573211;
    const double dg_kcal = sol.e_solvent() * kHaInKcal;
    const double pin_kcal = kPinSolvHa * kHaInKcal;
    std::printf("H2O dG_el = %.6f kcal/mol (pin %.6f)\n", dg_kcal, pin_kcal);
    if (!(std::abs(dg_kcal - pin_kcal) <= 0.1))
        return fail("H2O C-PCM free energy off PySCF");
    if (!(sol.mf.energy() < water_rhf.energy()))
        return fail("H2O C-PCM total not below gas phase");

    std::printf("PASS\n");
    return 0;
}
