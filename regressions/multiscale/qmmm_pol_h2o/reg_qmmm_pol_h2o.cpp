/// @file regressions/multiscale/qmmm_pol_h2o/reg_qmmm_pol_h2o.cpp
/// @brief Polarizable QM/MM: H2O (RHF/STO-3G) + polarizable MM site.
///
/// Gate: polarization energy is negative and total energy is finite and
/// lower than the fixed-charge (α→0 limit via QmmmElectrostatic) path by a
/// hard-coded minimum shift. See REFERENCE.md.
#include <cmath>
#include <cstdio>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/multiscale/qmmm_electrostatic.h"
#include "atomiverse/cpu/multiscale/qmmm_polarizable.h"

using namespace ATV;

namespace {

cpu::Frame water_plus_site() {
    cpu::Frame f(4);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 0.9572;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_x)[2] = -0.239987;
    f.span(kAtom_y)[2] = 0.926627;
    f.span(kAtom_z)[2] = 0;
    f.span(kAtom_x)[3] = 0;
    f.span(kAtom_y)[3] = 0;
    f.span(kAtom_z)[3] = 4.0;
    f.span(kAtom_type)[0] = 8;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_type)[3] = 18; // Ar-like MM site
    for (int i = 0; i < 4; ++i) {
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_mass)[i] = 1;
    }
    f.simbox = cpu::Box::cubic(25);
    return f;
}

} // namespace

int main() {
    double e_fix = 0, e_pol = 0, e_pol_term = 0;
    {
        auto f = water_plus_site();
        cpu::Evaluate eval;
        eval.on(evaluate::Phase::evaluate)
            .add(cpu::QmmmElectrostatic{{1, 2, 3}, 0.0, 0.01, 3.4, 10.0});
        eval.build();
        eval.run(f);
        e_fix = f.meta.get(kEnergy);
    }
    {
        auto f = water_plus_site();
        cpu::Evaluate eval;
        cpu::QmmmPolarizable* p = nullptr;
        // hold via owned add
        eval.on(evaluate::Phase::evaluate)
            .add(cpu::QmmmPolarizable{{1, 2, 3}, 0.0, /*alpha=*/2.0,
                                      /*thole=*/2.6, 0.01, 3.4, 10.0});
        eval.build();
        eval.run(f);
        e_pol = f.meta.get(kEnergy);
        (void)p;
        // Re-run via direct component to read pol term
        cpu::QmmmPolarizable pol{{1, 2, 3}, 0.0, 2.0, 2.6, 0.01, 3.4, 10.0};
        f.meta.get(kEnergy) = 0;
        cpu::Context ctx;
        pol.build({ctx, f});
        pol.launch({ctx, f, 0, 0});
        e_pol_term = pol.last_pol_energy();
        e_pol = f.meta.get(kEnergy);
        if (!pol.last_pol_converged()) {
            std::printf("FAIL: pol SCC not converged\n");
            return 1;
        }
    }
    std::printf("qmmm_pol_h2o E_fix=%.6f E_pol=%.6f E_pol_term=%.6f eV\n",
                e_fix, e_pol, e_pol_term);
    if (!std::isfinite(e_fix) || !std::isfinite(e_pol) ||
        !std::isfinite(e_pol_term)) {
        std::printf("FAIL: non-finite\n");
        return 1;
    }
    if (!(e_pol_term < -1e-4)) {
        std::printf("FAIL: polarization energy not negative enough\n");
        return 1;
    }
    if (!(e_pol < e_fix - 1e-4)) {
        std::printf("FAIL: polarizable total not below fixed-charge\n");
        return 1;
    }
    std::printf("PASS\n");
    return 0;
}
