/// @file benchmarks/md/lj/consistency/cpu_reference.cpp
/// @brief CPU reference half of the cpu↔cuda LJ backend-consistency
///        benchmark. Runs the canonical FCC-argon LJ force field on the CPU
///        (the reference backend, all-double) and dumps energy + per-atom
///        forces to a golden file for `cuda_compare.cu` to check against.
///
/// This realizes the manual-system §6 policy — "CPU is the reference; every
/// physical quantity has a CPU golden value" — as a committed parity case.
/// It is a plain executable, NOT a ctest gate (backend parity compares two
/// engines, so it lives in benchmarks/, like every other cross-engine check).
///
/// Usage:
///   bench_lj_consistency_cpu [--golden PATH] [--steps N] [--dt FS]
///
///   --golden PATH   output golden file (default: lj_consistency_golden.txt)
///   --steps N       integration steps before the snapshot (default: 0 =
///                   post-setup single force field, the tightest comparison)
///   --dt FS         timestep [fs] when --steps > 0 (default: 0.5)
///
/// Exit code 0 = golden written, 1 = failure (I/O, or a suspiciously trivial
/// zero force field — see the guard below).

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "consistency_system.h"

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"

using namespace ATV;
using namespace atv_consistency;

int main(int argc, char** argv) {
    std::string golden_path = "lj_consistency_golden.txt";
    ConsistencyParams p;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&]() -> std::string {
            return (i + 1 < argc) ? argv[++i] : "";
        };
        if (arg == "--golden")
            golden_path = next();
        else if (arg == "--steps")
            p.steps = std::atoi(next().c_str());
        else if (arg == "--dt")
            p.dt = std::atof(next().c_str());
        else {
            std::fprintf(stderr, "unknown argument: %s\n", arg.c_str());
            return 1;
        }
    }

    auto frame = make_consistency_system();

    cpu::md::MD sim;
    sim.set_integrator(cpu::md::VelocityVerlet{p.dt});
    sim.set_pair(cpu::PairLJCut{kEpsilon, kSigma, kCutoff});
    sim.set_neighbor(cpu::md::NeighborList<>{kCutoff});
    sim.build();
    sim.run(frame, p.steps);

    // Snapshot the reference energy + per-atom forces.
    Golden g;
    g.energy = frame.scalars.energy;
    const int n = frame.atoms.size();
    g.fx.reserve(n);
    g.fy.reserve(n);
    g.fz.reserve(n);
    double fmax = 0.0;
    for (int i = 0; i < n; ++i) {
        g.fx.push_back(frame.atoms.fx[i]);
        g.fy.push_back(frame.atoms.fy[i]);
        g.fz.push_back(frame.atoms.fz[i]);
        fmax = std::fmax(fmax, std::fabs(frame.atoms.fx[i]));
        fmax = std::fmax(fmax, std::fabs(frame.atoms.fy[i]));
        fmax = std::fmax(fmax, std::fabs(frame.atoms.fz[i]));
    }

    // Guard against a silent false-pass: if the force field is trivially
    // zero the driver did not compute forces at this phase point (e.g.
    // --steps 0 on a build that does not force-evaluate at setup). A golden
    // of all zeros would make the GPU compare pass vacuously. Fail loudly.
    if (fmax == 0.0) {
        std::fprintf(stderr,
                     "FAIL: reference force field is identically zero at "
                     "steps=%d — the CPU driver did not force-evaluate at "
                     "this phase point. Re-run with --steps 1.\n",
                     p.steps);
        return 1;
    }

    if (!write_golden(golden_path, g)) {
        std::fprintf(stderr, "FAIL: could not write golden file %s\n",
                     golden_path.c_str());
        return 1;
    }

    std::printf("# cpu-reference LJ consistency golden\n");
    std::printf("  atoms   = %d\n", n);
    std::printf("  steps   = %d (dt = %.3f fs)\n", p.steps, p.dt);
    std::printf("  energy  = %.10g eV\n", g.energy);
    std::printf("  |F|max  = %.10g eV/Angstrom\n", fmax);
    std::printf("  golden  -> %s\n", golden_path.c_str());
    return 0;
}
