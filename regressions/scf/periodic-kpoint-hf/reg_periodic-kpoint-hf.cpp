/// @file regressions/scf/periodic-kpoint-hf/reg_periodic-kpoint-hf.cpp
/// @brief H₂ / STO-3G KRHF, 4 Bohr cube, Monkhorst–Pack 1×1×2.
///
/// Public API only: `KRHF{STO3G, k_mesh}.run(frame)`.
/// Oracle: PySCF 2.14.0 `pbc.scf.KRHF` + default FFTDF / exxdiv='ewald'
/// (tools/pyscf_oracle/pin_scf.py krhf). See REFERENCE.md.

#include <cmath>
#include <cstdio>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/krhf.h"

int main() {
    using namespace ATV;

    // PySCF cell.unit = 'Bohr': lattice 4 Bohr, bond 1.4 Bohr.
    const double L_A = 4.0 * units::kBohrInAngstrom;
    const double z_A = 1.4 * units::kBohrInAngstrom;

    cpu::Frame frame(2);
    frame.span(kAtom_type)[0] = frame.span(kAtom_type)[1] = 1;
    frame.span(kAtom_x)[0] = frame.span(kAtom_y)[0] = frame.span(kAtom_z)[0] =
        0.0;
    frame.span(kAtom_x)[1] = frame.span(kAtom_y)[1] = 0.0;
    frame.span(kAtom_z)[1] = z_A;
    frame.box = cpu::Box::cubic(L_A);

    cpu::scf::KRHF krhf{cpu::scf::BasisSet::STO3G,
                        {1, 1, 2},
                        /*charge=*/0,
                        /*tol=*/1e-10,
                        /*max_iter=*/40};
    krhf.run(frame);

    if (!krhf.converged()) {
        std::fprintf(stderr, "FAIL: KRHF did not converge\n");
        return 1;
    }

    // PySCF 2.14.0 pbc.scf.KRHF.kernel on this cell.
    constexpr double kEpinHa = -1.3541438798330374;
    constexpr double kTolHa = 1.0e-6;

    const double e_ha = krhf.energy();
    const double diff = std::fabs(e_ha - kEpinHa);
    std::printf("E_KRHF(H2/STO-3G, 4 Bohr, 1x1x2) got=%.16f want=%.16f "
                "|diff|=%.3e tol=%.0e iters=%d\n",
                e_ha, kEpinHa, diff, kTolHa, krhf.iterations());
    if (!std::isfinite(e_ha) || diff > kTolHa) {
        std::fprintf(stderr, "FAIL: KRHF energy off pin\n");
        return 1;
    }
    const double e_ev = frame.meta.get(kEnergy);
    if (std::fabs(e_ev - e_ha * units::kHartreeInEv) > 1e-8) {
        std::fprintf(stderr, "FAIL: kEnergy [eV] does not match energy()\n");
        return 1;
    }
    std::printf("PASS\n");
    return 0;
}
