/// @file regressions/md/pim-induced-dipole-bulk/reg_pim-induced-dipole-bulk.cpp
/// @brief END-TO-END: bulk PIM induced dipoles vs offline MetalWalls-form pin.
///
/// See REFERENCE.md for geometry and pin values.
/// Public headers + link atomiverse only. Flat main(); exits non-zero on fail.

#include <cmath>
#include <cstdio>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/potentials/pim_induced_dipole.h"

using namespace ATV;

int main() {
    // Offline pin — Applequist/Thole bulk PIM (MetalWalls E4 form).
    // Geometry: q=+1 α=1.0 at origin; q=-1 α=1.5 at (3,0,0); L=50; thole_a=0.39.
    // See regressions/md/pim-induced-dipole-bulk/REFERENCE.md.
    const double r = 3.0;
    const double L = 50.0;
    const double thole_a = 0.39;
    const double tol_e = 1e-9;
    const double tol_f = 1e-8;

    const double E_pin = -0.12040946414707571; // [eV]
    const double f0x_pin = 0.040341072868572643;
    const double f1x_pin = -0.040341072868572643;

    cpu::Frame frame(2);
    frame.span(kAtom_x)[0] = 0.0;
    frame.span(kAtom_y)[0] = 0.0;
    frame.span(kAtom_z)[0] = 0.0;
    frame.span(kAtom_x)[1] = r;
    frame.span(kAtom_y)[1] = 0.0;
    frame.span(kAtom_z)[1] = 0.0;
    frame.span(kAtom_charge)[0] = 1.0;
    frame.span(kAtom_charge)[1] = -1.0;
    frame.span(kAtom_type)[0] = 1;
    frame.span(kAtom_type)[1] = 2;
    frame.span(kAtom_mass)[0] = frame.span(kAtom_mass)[1] = 1.0;
    frame.span(kAtom_tag)[0] = 1;
    frame.span(kAtom_tag)[1] = 2;
    frame.span(kAtom_fx)[0] = frame.span(kAtom_fy)[0] =
        frame.span(kAtom_fz)[0] = 0;
    frame.span(kAtom_fx)[1] = frame.span(kAtom_fy)[1] =
        frame.span(kAtom_fz)[1] = 0;
    frame.simbox = cpu::Box::cubic(L);
    frame.meta.get(kEnergy) = 0;

    cpu::PimInducedDipole pim(20.0, /*alpha_by_type=*/{0.0, 1.0, 1.5}, thole_a,
                              1e-12, 200);

    cpu::Evaluate eval;
    eval.on(evaluate::Phase::evaluate).add_ref(pim);
    eval.build();
    eval.run(frame);

    const double E = frame.meta.get(kEnergy);
    const double f0x = frame.span(kAtom_fx)[0];
    const double f0y = frame.span(kAtom_fy)[0];
    const double f0z = frame.span(kAtom_fz)[0];
    const double f1x = frame.span(kAtom_fx)[1];
    const double f1y = frame.span(kAtom_fy)[1];
    const double f1z = frame.span(kAtom_fz)[1];

    std::printf("PIM: E=%.16e want=%.16e\n", E, E_pin);
    std::printf("  F0=(%.16e, %.16e, %.16e)\n", f0x, f0y, f0z);
    std::printf("  F1=(%.16e, %.16e, %.16e)\n", f1x, f1y, f1z);
    std::printf("  e_pol=%.16e iters=%d\n", pim.last_pol_energy(),
                pim.last_pol_iters());

    bool ok = true;
    if (!(std::fabs(E - E_pin) <= tol_e)) {
        std::printf("FAIL energy got=%.16e want=%.16e tol=%.1e\n", E, E_pin,
                    tol_e);
        ok = false;
    }
    if (!(std::fabs(pim.last_pol_energy() - E_pin) <= tol_e)) {
        std::printf("FAIL last_pol_energy got=%.16e want=%.16e\n",
                    pim.last_pol_energy(), E_pin);
        ok = false;
    }
    if (!(std::fabs(f0x - f0x_pin) <= tol_f
          && std::fabs(f1x - f1x_pin) <= tol_f)) {
        std::printf(
            "FAIL force got F0x=%.16e F1x=%.16e want %.16e / %.16e tol=%.1e\n",
            f0x, f1x, f0x_pin, f1x_pin, tol_f);
        ok = false;
    }
    if (!(std::fabs(f0y) < 1e-12 && std::fabs(f0z) < 1e-12
          && std::fabs(f1y) < 1e-12 && std::fabs(f1z) < 1e-12)) {
        std::printf("FAIL transverse forces non-zero\n");
        ok = false;
    }
    if (!(std::fabs(f0x + f1x) < 1e-12)) {
        std::printf("FAIL action-reaction\n");
        ok = false;
    }
    if (!pim.last_pol_converged()) {
        std::printf("FAIL SCC not converged\n");
        ok = false;
    }
    if (!(E < 0.0)) {
        std::printf("FAIL E_pol should be negative\n");
        ok = false;
    }

    if (!ok)
        return 1;
    std::printf(
        "PASS: pim-induced-dipole-bulk vs offline MetalWalls-form pin\n");
    return 0;
}
