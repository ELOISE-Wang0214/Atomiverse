/// @file benchmarks/md/lj/consistency/cuda_compare.cu
/// @brief GPU half of the cpu↔cuda LJ backend-consistency benchmark. Replays
///        the IDENTICAL FCC-argon system on `cuda::md::MD`, reads the CPU
///        golden written by `cpu_reference.cpp`, and asserts per-quantity
///        mixed-precision tolerances.
///
/// **GPU-only / nightly.** The login nodes have no CUDA toolkit, so this
/// target builds and runs only on the GH200 node (`ATV_BUILD_CUDA=ON`). It is
/// deliberately excluded from the CPU CI matrix — wire it into the self-hosted
/// GPU/nightly job, not `ci.yml`.
///
/// **Why not bitwise.** The CUDA backend stores per-atom state at `real_t`
/// (float — see cuda/frame/types.cuh) and does the pair math at `work_t`, so
/// positions carry ~1e-7 relative float-storage error that the r^-12 LJ core
/// amplifies. The comparison is therefore to the *storage-precision floor*, not
/// bit-for-bit (manual-system §5 `backend-dependent` result class, §6 policy).
/// Accumulators (energy) are AccReal=double on both backends, but the CUDA
/// energy sums float-stored contributions, so it too carries the float floor.
///
/// Tolerance policy (per quantity; defaults are the initial float-storage
/// floor — CALIBRATE on the first real GPU run, then only ever TIGHTEN):
///   - energy:  |ΔE| ≤ etol_rel · |E_cpu|                 (default 2e-4)
///   - force:   |ΔFc| ≤ ftol_rel · |Fc_cpu| + ftol_abs    per component
///              (default ftol_rel 2e-3; ftol_abs auto = 1e-3 · F_rms, so a
///               near-zero component is bounded relative to the force scale)
///
/// Usage:
///   bench_lj_consistency_cuda [--golden PATH] [--steps N] [--dt FS]
///                             [--etol-rel R] [--ftol-rel R] [--ftol-abs A]
/// Exit code 0 = PASS (within tolerance), 1 = FAIL / golden missing.

#ifdef ATV_BUILD_CUDA

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "consistency_system.h"

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cuda/md/md.cuh"
#include "atomiverse/cuda/md/velocity_verlet.cuh"
#include "atomiverse/cuda/neighbor/neighbor_list.cuh"
#include "atomiverse/cuda/potentials/pair_lj_cut.cuh"

using namespace ATV;
using namespace atv_consistency;

int main(int argc, char** argv) {
    std::string golden_path = "lj_consistency_golden.txt";
    ConsistencyParams p;
    double etol_rel = 2e-4;
    double ftol_rel = 2e-3;
    double ftol_abs = -1.0; // <0 → auto = 1e-3 · F_rms

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
        else if (arg == "--etol-rel")
            etol_rel = std::atof(next().c_str());
        else if (arg == "--ftol-rel")
            ftol_rel = std::atof(next().c_str());
        else if (arg == "--ftol-abs")
            ftol_abs = std::atof(next().c_str());
        else {
            std::fprintf(stderr, "unknown argument: %s\n", arg.c_str());
            return 1;
        }
    }

    // ── Load the CPU golden ───────────────────────────────────────────
    Golden ref;
    if (!read_golden(golden_path, ref)) {
        std::fprintf(stderr,
                     "FAIL: could not read golden %s — run "
                     "bench_lj_consistency_cpu first (same --steps/--dt).\n",
                     golden_path.c_str());
        return 1;
    }

    // ── Replay the identical system on the CUDA backend ───────────────
    // `cuda::md::MD::run(cpu::Frame&)` materializes to device and writes the
    // result back to the host frame (cf. bench_md_step.cu). NeighborList /
    // PairLJCut live in `cuda::` (not `cuda::md::`); VelocityVerlet in
    // `cuda::md::`. Bind explicitly to avoid the `md` namespace ambiguity.
    auto frame = make_consistency_system();
    cuda::md::MD sim;
    sim.add(cuda::md::VelocityVerlet{p.dt});
    sim.on(ATV::md::Phase::force)
        .add(cuda::NeighborList{kCutoff}, -10)
        .add(cuda::PairLJCut{kEpsilon, kSigma, kCutoff});
    sim.build();
    sim.run(frame, p.steps);
    cudaDeviceSynchronize();

    const int n = frame.atoms.size();
    if (n != ref.n()) {
        std::fprintf(stderr, "FAIL: atom count mismatch (cuda %d vs golden %d)\n",
                     n, ref.n());
        return 1;
    }

    // ── Resolve the absolute force floor from the reference RMS ───────
    double f2sum = 0.0;
    for (int i = 0; i < n; ++i)
        f2sum += ref.fx[i] * ref.fx[i] + ref.fy[i] * ref.fy[i]
                 + ref.fz[i] * ref.fz[i];
    const double f_rms = std::sqrt(f2sum / (3.0 * n));
    if (ftol_abs < 0.0)
        ftol_abs = 1e-3 * f_rms;

    // ── Compare ───────────────────────────────────────────────────────
    double e_cuda = frame.scalars.energy;
    const double e_abs = std::fabs(e_cuda - ref.energy);
    const double e_rel =
        (ref.energy != 0.0) ? e_abs / std::fabs(ref.energy) : e_abs;

    double max_ferr = 0.0, max_ferr_rel = 0.0, ferr2 = 0.0;
    int n_bad = 0;
    double fmax_cuda = 0.0;
    auto check = [&](double got, double want) {
        const double d = std::fabs(got - want);
        const double allow = ftol_rel * std::fabs(want) + ftol_abs;
        max_ferr = std::fmax(max_ferr, d);
        if (std::fabs(want) > ftol_abs)
            max_ferr_rel = std::fmax(max_ferr_rel, d / std::fabs(want));
        ferr2 += d * d;
        fmax_cuda = std::fmax(fmax_cuda, std::fabs(got));
        if (d > allow)
            ++n_bad;
    };
    for (int i = 0; i < n; ++i) {
        check(frame.atoms.fx[i], ref.fx[i]);
        check(frame.atoms.fy[i], ref.fy[i]);
        check(frame.atoms.fz[i], ref.fz[i]);
    }
    const double ferr_rms = std::sqrt(ferr2 / (3.0 * n));

    // Guard: an all-zero CUDA force field means the driver did not write
    // forces back to the host frame at this phase point — a vacuous pass.
    if (fmax_cuda == 0.0) {
        std::fprintf(stderr,
                     "FAIL: CUDA force field read back as identically zero — "
                     "cuda::md::MD may not write forces D2H at steps=%d "
                     "(energy write-back is unaffected). Verify the transfer "
                     "path or compare at --steps 1.\n",
                     p.steps);
        return 1;
    }

    const bool energy_ok = e_rel <= etol_rel;
    const bool forces_ok = n_bad == 0;

    std::printf("# cpu↔cuda LJ backend consistency (mixed precision)\n");
    std::printf("  atoms          = %d   steps = %d (dt = %.3f fs)\n", n,
                p.steps, p.dt);
    std::printf("  F_rms (ref)    = %.6g eV/Angstrom\n", f_rms);
    std::printf("  energy cpu     = %.10g eV\n", ref.energy);
    std::printf("  energy cuda    = %.10g eV\n", e_cuda);
    std::printf("  |dE|           = %.3e eV   rel %.3e   (tol_rel %.1e)  %s\n",
                e_abs, e_rel, etol_rel, energy_ok ? "OK" : "FAIL");
    std::printf("  max |dF|       = %.3e eV/Angstrom   max rel %.3e\n",
                max_ferr, max_ferr_rel);
    std::printf("  rms |dF|       = %.3e eV/Angstrom\n", ferr_rms);
    std::printf("  force tol      = %.1e·|F| + %.3e   (%d/%d components out)\n",
                ftol_rel, ftol_abs, n_bad, 3 * n);

    const bool pass = energy_ok && forces_ok;
    std::printf("%s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}

#else
#include <cstdio>
int main() {
    std::printf("bench_lj_consistency_cuda requires -DATV_BUILD_CUDA=ON "
                "(GPU node only)\n");
    return 1;
}
#endif
