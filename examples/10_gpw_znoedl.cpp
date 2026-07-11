/// @file examples/10_gpw_znoedl.cpp
/// @brief Full-system GPW single point for the ZnO electric-double-layer (EDL)
///        slab — the gpw-edl-aimd capstone.
///
/// Loads the 1214-atom ZnO-slab + water + NaCl electrolyte frame
/// (`final_frame.xyz`, extended XYZ with a `Lattice=` header) through molrs via
/// `cpu::load_extxyz`, configures GPW for all five elements (Zn/O/H/Na/Cl) with
/// GTH-PBE pseudopotentials + the DZVP-MOLOPT-SR-GTH basis (exactly CP2K's
/// `ZnO_MD.inp` electronic structure), and runs a periodic PBE single point.
///
/// This is the production target of the gpw-edl-aimd spec. With analytic forces
/// + CSVR NVT (shown at the bottom, not run here) it becomes the ZnO_MD AIMD:
///
///   cpu::md::MD aimd;
///   aimd.set_nvt(cpu::md::NVT{0.5, 330.0, 20.0, cpu::md::NVT::Thermostat::CSVR});
///   cpu::gpw::GpwEvaluatorConfig fcfg = cfg; fcfg.compute_forces = true;
///   aimd.on(ATV::md::Phase::force).add(cpu::gpw::GPWEvaluator{fcfg});
///   aimd.build(); aimd.run(frame, nsteps);
///
/// Build:
///   cmake -B build -DATV_BUILD_RHF=ON -DATV_BUILD_PW=ON -DATV_BUILD_GPW=ON \
///                  -DATV_BUILD_EXAMPLES=ON
///   cmake --build build --target 10_gpw_znoedl
///
/// Run (args: [cutoff_ry] [z_lo z_hi] — the optional z-window in Å
/// cuts a slab fragment for scaling studies; default = all 1214 atoms):
///   ./build/examples/10_gpw_znoedl 280

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "atomiverse/cpu/driver/single_point.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/basis_data/gpw_data.h"
#include "atomiverse/cpu/io/extxyz.h"

using namespace ATV;

namespace {

const char* kDefaultFrame =
    "/nobackup/proj/disk/teoroo/personal/jicli594/work/molpack-paper/figures/"
    "final_frame.xyz";

/// @brief Keep only atoms whose z (Å) is in [z_lo, z_hi]; preserves the box.
cpu::Frame z_slab(const cpu::Frame& src, double z_lo, double z_hi) {
    std::vector<int> keep;
    for (int i = 0; i < src.atoms.size(); ++i)
        if (src.atoms.z[i] >= z_lo && src.atoms.z[i] <= z_hi) keep.push_back(i);
    cpu::Frame f(static_cast<int>(keep.size()));
    for (std::size_t k = 0; k < keep.size(); ++k) {
        const int i = keep[k];
        f.atoms.type[k] = src.atoms.type[i];
        f.atoms.x[k] = src.atoms.x[i];
        f.atoms.y[k] = src.atoms.y[i];
        f.atoms.z[k] = src.atoms.z[i];
    }
    f.box = src.box;
    return f;
}

/// @brief One GpwElement (Z, DZVP-MOLOPT-SR-GTH basis, GTH-PBE-q* potential).
cpu::gpw::GpwElement element(int z, const char* sym, const char* pot) {
    cpu::gpw::GpwElement e;
    e.z = z;
    e.basis = cpu::gpw::builtin_gth_basis(sym, "DZVP-MOLOPT-SR-GTH");
    e.pot = cpu::gpw::builtin_gth_potential(sym, pot);
    return e;
}

} // namespace

int main(int argc, char** argv) {
    // No data-directory argument: the GTH/MOLOPT tables are in-tree.
    const double cutoff = argc > 1 ? std::atof(argv[1]) : 280.0;

    cpu::Frame frame = cpu::load_extxyz(kDefaultFrame);
    if (argc > 3) frame = z_slab(frame, std::atof(argv[2]), std::atof(argv[3]));

    std::cout << "ZnO-EDL GPW single point\n"
              << "  atoms   : " << frame.atoms.size() << "\n"
              << "  box [A] : " << frame.box.length(0) << " x " << frame.box.length(1)
              << " x " << frame.box.length(2) << "\n"
              << "  cutoff  : " << cutoff << " Ry (CP2K ZnO_MD.inp uses 320)\n"
              << std::flush;

    cpu::gpw::GpwEvaluatorConfig cfg;
    cfg.cutoff_ry = cutoff;
    cfg.xc = cpu::gpw::XcKind::PBE;
    cfg.nonperiodic = false; // periodic XYZ (slab)
    cfg.max_iter = 100;
    cfg.tol = 1e-6;
    cfg.elements = {
        element(30, "Zn", "GTH-PBE-q12"), element(8, "O", "GTH-PBE-q6"),
        element(1, "H", "GTH-PBE-q1"),    element(11, "Na", "GTH-PBE-q9"),
        element(17, "Cl", "GTH-PBE-q7"),
    };

    cpu::gpw::GPWEvaluator eval(cfg);

    const auto t0 = std::chrono::steady_clock::now();
    cpu::single_point(frame, eval);
    const auto t1 = std::chrono::steady_clock::now();
    const double secs =
        std::chrono::duration<double>(t1 - t0).count();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  converged : " << (eval.last_result().converged ? "yes" : "no")
              << " in " << eval.last_result().iters << " SCF iters\n"
              << "  E_total   : " << eval.last_energy_hartree() << " Hartree   ("
              << frame.scalars.energy << " eV)\n"
              << "  wall-time : " << secs << " s\n";
    return 0;
}
