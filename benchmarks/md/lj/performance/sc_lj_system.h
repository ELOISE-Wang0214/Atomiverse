/// @file benchmarks/md/lj/performance/sc_lj_system.h
/// @brief Shared LJ benchmark system + timing/RSS helpers for the CPU MPI
///        scaling drivers (`scaling.cpp` MD, `minimize_scaling.cpp` minimizer).
///
/// One source of truth for the benchmark lattice: simple-cubic, sigma = 2.5,
/// spacing = 1.12*sigma (rho* = 0.712), 0.1-A deterministic sinusoidal jitter,
/// matching the GPU benchmark (bench_multi_gpu_scaling.cu) and the LAMMPS
/// decks (in.lj_kk / in.lj_min_kk) so every panel of the study runs the same
/// system with the same neighbors/atom (~24).
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#else
#include <chrono>
#endif

namespace atv_bench {

inline constexpr double kSigma = 2.5;
inline constexpr double kSpacing = 1.12 * kSigma; // rho* = (1/1.12)^3 = 0.712

/// @brief Cube-rounded edge cell count for a target total atom count.
inline int cells_for(int n_target) {
    return std::max(
        1, static_cast<int>(std::cbrt(static_cast<double>(n_target)) + 0.5));
}

// Build a simple-cubic LJ system with `n_target` total atoms (rounded to the
// nearest cube), IDENTICAL to the GPU benchmark (bench_multi_gpu_scaling):
// sigma = 2.5, spacing = 1.12*sigma -> reduced density rho* = 0.712, so the
// CPU and GPU panels of the 2x2 study run the same LJ system. Mild
// deterministic jitter + thermal velocities (not a frozen lattice).
inline ATV::cpu::Frame make_sc_lj_total(int n_target) {
    using namespace ATV;
    const int cells = cells_for(n_target);
    const int n = cells * cells * cells;
    cpu::Frame frame(n);
    int idx = 0;
    for (int ix = 0; ix < cells; ++ix)
        for (int iy = 0; iy < cells; ++iy)
            for (int iz = 0; iz < cells; ++iz) {
                frame.span(ATV::kAtom_x)[idx] =
                    (ix + 0.5) * kSpacing + 0.1 * std::sin(0.7 * idx);
                frame.span(ATV::kAtom_y)[idx] =
                    (iy + 0.5) * kSpacing + 0.1 * std::sin(1.3 * idx);
                frame.span(ATV::kAtom_z)[idx] =
                    (iz + 0.5) * kSpacing + 0.1 * std::sin(1.9 * idx);
                frame.span(ATV::kAtom_mass)[idx] = 1.0;
                frame.span(ATV::kAtom_type)[idx] = 0;
                frame.span(ATV::kAtom_vx)[idx] = 0.01 * std::sin(0.5 * idx);
                frame.span(ATV::kAtom_vy)[idx] = 0.01 * std::sin(0.9 * idx);
                frame.span(ATV::kAtom_vz)[idx] = 0.01 * std::sin(1.1 * idx);
                ++idx;
            }
    frame.simbox = ATV::cpu::Box::cubic(cells * kSpacing);
    return frame;
}

#ifdef ATV_HAS_MPI
// Build ONLY this rank's owned atoms of the same global simple-cubic LJ lattice
// as make_sc_lj_total — no rank ever holds the global frame (cf. LAMMPS
// create_atoms, which fills each proc's subdomain locally). This is the memory
// fix: the old "rank 0 builds global + scatter" path peaked at O(global-N) on
// rank 0 (5.9 GB at 256 ranks weak, 23 GB at 1024 → the dominant per-rank RSS).
//
// Ownership MUST match DomainDecomposition exactly or the first partition_atoms
// drops boundary atoms: same MPI_Dims_create / ATV_FORCE_DIMS, same rank index
// convention (ix*dy*dz + iy*dz + iz), same rank_along_dim = clamp(int(q*d/L)).
// Global atom index `idx` (lattice order) drives the identical jitter, velocity
// and tag, so the distributed system is bit-identical to the single-rank one.
inline ATV::cpu::Frame make_sc_lj_local(int n_target, int my_rank,
                                        int world_size) {
    using namespace ATV;
    const int cells = cells_for(n_target);
    const double L = cells * kSpacing;

    int dims[3] = {0, 0, 0};
    if (const char* env = std::getenv("ATV_FORCE_DIMS")) {
        int dx, dy, dz;
        if (std::sscanf(env, "%d,%d,%d", &dx, &dy, &dz) == 3
            && dx * dy * dz == world_size) {
            dims[0] = dx;
            dims[1] = dy;
            dims[2] = dz;
        }
    }
    if (dims[0] == 0)
        MPI_Dims_create(world_size, 3, dims);
    const int coords[3] = {my_rank / (dims[1] * dims[2]),
                           (my_rank / dims[2]) % dims[1], my_rank % dims[2]};

    auto owner = [&](double q, int d) {
        q -= L * std::floor(q / L); // wrap to [0, L) like partition_atoms
        if (d <= 1)
            return 0;
        return std::max(0, std::min(d - 1, static_cast<int>(q * d / L)));
    };
    // Local lattice-index window (+/-2 margin covers the 0.1-A jitter).
    auto lo = [&](int c, int d) {
        return std::max(0, static_cast<int>(static_cast<double>(c) * cells / d)
                               - 2);
    };
    auto hi = [&](int c, int d) {
        return std::min(cells,
                        static_cast<int>(static_cast<double>(c + 1) * cells / d)
                            + 2);
    };

    std::vector<double> px, py, pz, vx, vy, vz;
    std::vector<std::int64_t> tg;
    for (int ix = lo(coords[0], dims[0]); ix < hi(coords[0], dims[0]); ++ix)
        for (int iy = lo(coords[1], dims[1]); iy < hi(coords[1], dims[1]); ++iy)
            for (int iz = lo(coords[2], dims[2]); iz < hi(coords[2], dims[2]);
                 ++iz) {
                const std::int64_t idx =
                    (static_cast<std::int64_t>(ix) * cells + iy) * cells + iz;
                const double x =
                    (ix + 0.5) * kSpacing + 0.1 * std::sin(0.7 * idx);
                const double y =
                    (iy + 0.5) * kSpacing + 0.1 * std::sin(1.3 * idx);
                const double z =
                    (iz + 0.5) * kSpacing + 0.1 * std::sin(1.9 * idx);
                if (owner(x, dims[0]) != coords[0]
                    || owner(y, dims[1]) != coords[1]
                    || owner(z, dims[2]) != coords[2])
                    continue;
                px.push_back(x);
                py.push_back(y);
                pz.push_back(z);
                vx.push_back(0.01 * std::sin(0.5 * idx));
                vy.push_back(0.01 * std::sin(0.9 * idx));
                vz.push_back(0.01 * std::sin(1.1 * idx));
                tg.push_back(idx);
            }

    const int n = static_cast<int>(px.size());
    cpu::Frame frame(n);
    for (int i = 0; i < n; ++i) {
        frame.span(ATV::kAtom_x)[i] = px[i];
        frame.span(ATV::kAtom_y)[i] = py[i];
        frame.span(ATV::kAtom_z)[i] = pz[i];
        frame.span(ATV::kAtom_vx)[i] = vx[i];
        frame.span(ATV::kAtom_vy)[i] = vy[i];
        frame.span(ATV::kAtom_vz)[i] = vz[i];
        frame.span(ATV::kAtom_mass)[i] = 1.0;
        frame.span(ATV::kAtom_type)[i] = 0;
        frame.span(ATV::kAtom_tag)[i] =
            tg[i]; // GLOBAL tag (ctor default is local row i)
    }
    frame.set_nlocal(n);
    frame.clear_tag_index();
    frame.rebuild_tag_index(0, n);
    frame.simbox = ATV::cpu::Box::cubic(L);
    return frame;
}
#endif

inline double now_seconds() {
#ifdef ATV_HAS_MPI
    return MPI_Wtime();
#else
    using clk = std::chrono::steady_clock;
    static const auto t0 = clk::now();
    const auto t = clk::now();
    return std::chrono::duration<double>(t - t0).count();
#endif
}

// Peak resident set size of this process so far [kB] (Linux /proc/self/status
// VmHWM). Summed across ranks this is the benchmark's total memory footprint.
inline long read_vmhwm_kb() {
    FILE* f = std::fopen("/proc/self/status", "r");
    if (!f)
        return 0;
    char line[256];
    long kb = 0;
    while (std::fgets(line, sizeof(line), f))
        if (std::sscanf(line, "VmHWM: %ld kB", &kb) == 1)
            break;
    std::fclose(f);
    return kb;
}

} // namespace atv_bench
