/// @file regressions/scf/mean-field-method-split/reg_mean-field-method-split.cpp
/// @brief Mean-field split on classic systems, pinned to offline PySCF.
///
/// H₂/STO-3G: `RHF::run` stays. TDHF, CIS(D), and G0W0 do not
/// (`TDHF<RHF>`, `CISD<RHF>`, `G0W0<RHF>`).
/// H₂O/6-31G: `RHF::forces` stays on the mean-field type.
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
    std::printf("H2 E_RHF = %.12f Ha\n", rhf.energy());

    const cpu::scf::ExcitedStatesResult tdhf =
        cpu::scf::TDHF<cpu::scf::RHF>{rhf}.run(1);
    constexpr double kTdhfEv = 25.33200837852095;
    constexpr double kTdhfF = 0.8857897216216963;
    if (tdhf.energy_ev.size() != 1 || tdhf.tda)
        return fail("H2 TDHF root count or TDA flag");
    std::printf("H2 TDHF w=%.8f eV f=%.8f\n", tdhf.energy_ev[0],
                tdhf.oscillator_strength[0]);
    if (!(std::abs(tdhf.energy_ev[0] - kTdhfEv) <= 1e-4))
        return fail("H2 TDHF energy off PySCF");
    if (!(std::abs(tdhf.oscillator_strength[0] - kTdhfF) <= 1e-3))
        return fail("H2 TDHF oscillator strength off PySCF");

    const cpu::scf::CISDResult cisd = cpu::scf::CISD<cpu::scf::RHF>{rhf}.run(1);
    constexpr double kCisEv = 25.807465778338745;
    constexpr double kCisdEv = 26.16497097331174;
    if (cisd.energy_ev.size() != 1)
        return fail("H2 CIS(D) root count");
    std::printf("H2 CIS w=%.8f eV CIS(D)=%.8f eV\n", cisd.cis_energy_ev[0],
                cisd.energy_ev[0]);
    if (!(std::abs(cisd.cis_energy_ev[0] - kCisEv) <= 1e-4))
        return fail("H2 CIS energy off PySCF");
    if (!(std::abs(cisd.energy_ev[0] - kCisdEv) <= 1e-4))
        return fail("H2 CIS(D) energy off PySCF");

    const cpu::scf::G0W0Result gw = cpu::scf::G0W0<cpu::scf::RHF>{rhf}.run();
    constexpr double kGapEv = 35.007583892679726;
    constexpr double kIpEv = 16.2440458258512;
    constexpr double kEaEv = 18.763538066828527;
    constexpr double kQpHomo = -0.5969576734903094;
    constexpr double kQpLumo = 0.6895473055730482;
    if (!gw.converged || gw.qp_energy_ha.size() != 2)
        return fail("H2 G0W0 did not return two QP energies");
    std::printf("H2 G0W0 gap=%.8f eV IP=%.8f EA=%.8f\n", gw.gap_ev, gw.ip_ev,
                gw.ea_ev);
    if (!(std::abs(gw.gap_ev - kGapEv) <= 1e-4))
        return fail("H2 G0W0 gap off PySCF");
    if (!(std::abs(gw.ip_ev - kIpEv) <= 1e-4))
        return fail("H2 G0W0 IP off PySCF");
    if (!(std::abs(gw.ea_ev - kEaEv) <= 1e-4))
        return fail("H2 G0W0 EA off PySCF");
    if (!(std::abs(gw.qp_energy_ha[0] - kQpHomo) <= 1e-5))
        return fail("H2 QP HOMO off PySCF");
    if (!(std::abs(gw.qp_energy_ha[1] - kQpLumo) <= 1e-5))
        return fail("H2 QP LUMO off PySCF");

    auto water = h2o();
    cpu::scf::RHF water_rhf{cpu::scf::BasisSet::Pople631G, 0, 1e-12};
    water_rhf.max_iter = 100;
    const cpu::scf::GradientResult wg = water_rhf.forces(water);
    if (!water_rhf.converged())
        return fail("H2O RHF did not converge");
    auto wgrad = [&](int atom, int axis) {
        return -units::convert(wg.forces(atom, axis),
                               units::Force::eVPerAngstrom,
                               units::Force::HartreePerBohr);
    };
    // PySCF rhf --grad, Ha/Bohr, atom-major (O, H, H).
    constexpr double kWGrad[3][3] = {
        {-7.996236897155793e-18, -5.316927453868914e-16, 0.02381868594901171},
        {1.192691453658455e-16, -0.004413880349737553, -0.011909342974508963},
        {-1.1127290846869028e-16, 0.004413880349738886, -0.011909342974507409},
    };
    double max_wg = 0.0;
    for (int atom = 0; atom < 3; ++atom)
        for (int axis = 0; axis < 3; ++axis)
            max_wg = std::max(max_wg,
                              std::abs(wgrad(atom, axis) - kWGrad[atom][axis]));
    std::printf("H2O |grad err| = %.3e Ha/Bohr\n", max_wg);
    if (!(max_wg <= 1e-6))
        return fail("H2O RHF gradient off PySCF");

    std::printf("PASS\n");
    return 0;
}
