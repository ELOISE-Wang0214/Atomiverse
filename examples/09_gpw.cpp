/// @file examples/09_gpw.cpp
/// @brief GPW DFT single-point energy for an isolated water molecule.
///
/// The third electronic-structure family in these tutorials: GPW (Gaussian and
/// Plane Waves). The Kohn-Sham problem is solved in a contracted Gaussian basis
/// (like cpu::scf) while the Hartree + XC potentials are evaluated on a uniform
/// real-space grid via FFT (like cpu::pw). This is the method behind CP2K's
/// Quickstep, so it reproduces CP2K's numbers for the same input.
///
/// This mirrors CP2K's getting-started h2o.inp: isolated H2O, PBE,
/// DZVP-MOLOPT-GTH basis + GTH-PBE pseudopotentials, Martyna-Tuckerman isolated
/// electrostatics. The user-facing component is cpu::gpw::GPWEvaluator (role
/// Evaluator, capability Energy), run through cpu::single_point exactly like
/// cpu::pw::PWEvaluator and cpu::scf::SCFEvaluator.
///
/// Build:
///   cmake -B build -DATV_BUILD_RHF=ON -DATV_BUILD_PW=ON -DATV_BUILD_GPW=ON \
///                  -DATV_BUILD_EXAMPLES=ON                 # needs FFTW3
///   cmake --build build --target 09_gpw
///
/// Run (point at a CP2K data/ directory with BASIS_MOLOPT + GTH_POTENTIALS):
///   ./build/examples/09_gpw /path/to/cp2k/data
///
/// CP2K reference total energy for this input: -17.21233714 Hartree.

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "atomiverse/cpu/driver/single_point.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"

using namespace ATV;

namespace {

/// @brief Water molecule from the CP2K h2o.inp tutorial (coordinates in Å).
cpu::Frame make_water() {
    cpu::Frame f(3);
    const int z[3] = {8, 1, 1};
    const double x[3] = {5.0000, 5.7586, 4.2414};
    const double y[3] = {5.0000, 5.0000, 5.0000};
    const double zc[3] = {5.0000, 5.5043, 5.5043};
    for (int i = 0; i < 3; ++i) {
        f.atoms.type[i] = z[i];
        f.atoms.x[i] = x[i];
        f.atoms.y[i] = y[i];
        f.atoms.z[i] = zc[i];
    }
    f.box = cpu::Box::cubic(10.0); // 10 Å box
    return f;
}

} // namespace

int main(int argc, char** argv) {
    const std::string data = argc > 1 ? argv[1]
                                      : "/nobackup/proj/disk/teoroo/personal/"
                                        "jicli594/work/cp2k/data";

    cpu::Frame frame = make_water();

    cpu::gpw::GpwElement O, H;
    O.z = 8;
    O.basis = cpu::gpw::parse_gth_basis(data + "/BASIS_MOLOPT", "O", "DZVP-MOLOPT-GTH");
    O.pot = cpu::gpw::parse_gth_potential(data + "/GTH_POTENTIALS", "O", "GTH-PBE-q6");
    H.z = 1;
    H.basis = cpu::gpw::parse_gth_basis(data + "/BASIS_MOLOPT", "H", "DZVP-MOLOPT-GTH");
    H.pot = cpu::gpw::parse_gth_potential(data + "/GTH_POTENTIALS", "H", "GTH-PBE-q1");

    cpu::gpw::GpwEvaluatorConfig cfg;
    cfg.cutoff_ry = 400.0;
    cfg.xc = cpu::gpw::XcKind::PBE;
    cfg.nonperiodic = true; // isolated molecule (Martyna-Tuckerman)
    cfg.elements = {O, H};

    cpu::gpw::GPWEvaluator eval(cfg);
    cpu::single_point(frame, eval);

    const double e_ha = eval.last_energy_hartree();
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "GPW H2O single-point (PBE, DZVP-MOLOPT-GTH, MT)\n"
              << "  converged : " << (eval.last_result().converged ? "yes" : "no")
              << " in " << eval.last_result().iters << " iters\n"
              << "  E_total   : " << e_ha << " Hartree   (" << frame.scalars.energy
              << " eV)\n"
              << "  CP2K ref  : -17.21233714 Hartree\n"
              << "  |dE|      : " << std::abs(e_ha + 17.21233714458516) << " Hartree\n";
    return 0;
}
