/// @file regressions/scf/hf-complex-spinor/reg_hf-complex-spinor.cpp
/// @brief H2 and H2O STO-3G complex spin-orbital HF vs PySCF GHF.x2c1e.
///
/// Public API only. Pins: REFERENCE.md (PySCF 2.14.0, xuncontract=False).

#include <cmath>
#include <cstdio>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"

namespace {

ATV::cpu::Frame h2() {
    ATV::cpu::Frame frame(2);
    frame.span(ATV::kAtom_type)[0] = 1;
    frame.span(ATV::kAtom_type)[1] = 1;
    frame.span(ATV::kAtom_z)[1] = 0.74;
    frame.box = ATV::cpu::Box::cubic(30.0);
    return frame;
}

ATV::cpu::Frame h2o() {
    constexpr double kPi = 3.14159265358979323846;
    const double half = 0.5 * 104.5 * kPi / 180.0;
    const double hx = 0.95 * std::sin(half);
    const double hz = 0.95 * std::cos(half);
    ATV::cpu::Frame frame(3);
    frame.span(ATV::kAtom_type)[0] = 8;
    frame.span(ATV::kAtom_type)[1] = 1;
    frame.span(ATV::kAtom_type)[2] = 1;
    frame.span(ATV::kAtom_x)[1] = hx;
    frame.span(ATV::kAtom_z)[1] = hz;
    frame.span(ATV::kAtom_x)[2] = -hx;
    frame.span(ATV::kAtom_z)[2] = hz;
    frame.box = ATV::cpu::Box::cubic(30.0);
    return frame;
}

bool check(const char* label, ATV::cpu::Frame frame, bool x2c, double pin) {
    ATV::cpu::scf::ComplexGHF ghf{ATV::cpu::scf::BasisSet::STO3G, 0, 1e-12, 80};
    ghf.x2c = x2c;
    ghf.run(frame);
    const double e = ghf.energy();
    const double diff = std::abs(e - pin);
    std::printf("%s: E = %.12f Ha (pin %.12f, |dE| = %.3e, iters %d)\n", label,
                e, pin, diff, ghf.iterations());
    if (!ghf.converged()) {
        std::printf("FAIL: %s did not converge\n", label);
        return false;
    }
    const double e_ev = frame.meta.get(ATV::kEnergy) / ATV::units::kHartreeInEv;
    if (std::abs(e_ev - e) > 1e-12) {
        std::printf("FAIL: %s kEnergy does not match energy()\n", label);
        return false;
    }
    if (!(diff <= 1e-8)) {
        std::printf("FAIL: %s off PySCF pin\n", label);
        return false;
    }
    return true;
}

} // namespace

int main() {
    bool ok = true;
    ok = check("H2 x2c", h2(), true, -1.11677137468639) && ok;
    ok = check("H2O x2c", h2o(), true, -75.00458277099555) && ok;
    ok = check("H2O GHF", h2o(), false, -74.96175405132141) && ok;
    std::printf("hf-complex-spinor: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
