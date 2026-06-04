/// @file cpu/md/force_scratch.h
/// @brief Shared scaffolding for bonded/pair force kernels: per-thread force
///        accumulators and a Voigt-order virial accumulator.
///
/// Every bonded and pair kernel previously copy-pasted: (1) a triple of
/// per-thread `std::vector<double>` force buffers allocated in build(),
/// zeroed at the top of compute(), scattered into in the hot loop (so two
/// threads never write the same row), and reduced back into `frame` forces;
/// and (2) six `v_xx..v_yz` virial scalars + an unpack into `stress[0..5]`.
/// This header factors both out.
///
/// CRITICAL: the zero/reduce passes iterate `n_force = n_local + n_ghost`
/// (from `ResolveView::n_force()`), NOT `atoms.n_total()`. In multi-rank
/// `atoms.nlocal` is the global atom count, so the old `n_total()`-bounded
/// passes did O(global N) work on every rank — the measured ~100x bonded
/// straggler. Buffers are still *allocated* to `nmax` (capacity), but only
/// the live `[0, n_force)` rows are touched per step.
///
/// Non-virtual by design: `Bond`/`Angle`/`Dihedral` forbid virtual calls in
/// the per-tuple loop, so this is a plain owned member, not a base mixin.
#pragma once

#include <algorithm>
#include <vector>

#include "atomiverse/cpu/frame.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

/// @brief Raw per-thread force pointers for one OMP thread's hot loop.
struct ThreadForce {
    double* __restrict fx;
    double* __restrict fy;
    double* __restrict fz;
};

/// @brief Owns per-thread force accumulators; zero/reduce sized to n_force.
class BondedForceScratch {
public:
    /// @brief Allocate `n_threads` buffers of `nmax` rows each (call in build()).
    void build(int n_threads, int nmax) {
        n_threads_ = n_threads < 1 ? 1 : n_threads;
        fx_.assign(static_cast<std::size_t>(n_threads_),
                   std::vector<double>(static_cast<std::size_t>(nmax), 0.0));
        fy_ = fx_;
        fz_ = fx_;
    }

    int n_threads() const noexcept { return n_threads_; }

    /// @brief Zero thread `tid`'s `[0, n_force)` rows and return its pointers.
    ///        Call once per thread at the top of the parallel region.
    ThreadForce begin_thread(int tid, int n_force) {
        auto& fx = fx_[static_cast<std::size_t>(tid)];
        auto& fy = fy_[static_cast<std::size_t>(tid)];
        auto& fz = fz_[static_cast<std::size_t>(tid)];
        std::fill(fx.begin(), fx.begin() + n_force, 0.0);
        std::fill(fy.begin(), fy.begin() + n_force, 0.0);
        std::fill(fz.begin(), fz.begin() + n_force, 0.0);
        return {fx.data(), fy.data(), fz.data()};
    }

    /// @brief Reduce all threads' `[0, n_force)` contributions into frame
    ///        forces. Call after the parallel region. The owned + ghost rows
    ///        are both reduced (ghost forces ship home via HaloExchangeReverse).
    void reduce_into(cpu::Frame& f, int n_force) const {
        double* __restrict fx = f.atoms.fx.data();
        double* __restrict fy = f.atoms.fy.data();
        double* __restrict fz = f.atoms.fz.data();
        // Parallelize over rows (each row independent — no cross-thread write
        // hazard), gathering every thread's contribution for row i.
#pragma omp parallel for schedule(static)
        for (int i = 0; i < n_force; ++i) {
            double sx = 0.0, sy = 0.0, sz = 0.0;
            for (int t = 0; t < n_threads_; ++t) {
                sx += fx_[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)];
                sy += fy_[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)];
                sz += fz_[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)];
            }
            fx[i] += sx;
            fy[i] += sy;
            fz[i] += sz;
        }
    }

    void clear() {
        fx_.clear();
        fy_.clear();
        fz_.clear();
        n_threads_ = 0;
    }

private:
    int n_threads_ = 0;
    std::vector<std::vector<double>> fx_, fy_, fz_;
};

/// @brief Voigt-order (xx yy zz xy xz yz) virial accumulator. Supports
///        `#pragma omp declare reduction(+ : Voigt6 : omp_out += omp_in)`
///        so kernels reduce one object instead of six scalars.
struct Voigt6 {
    double v[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    /// @brief Add a pair contribution Δ⊗f in Voigt order.
    void add_pair(double dx, double dy, double dz, double fx, double fy,
                  double fz) noexcept {
        v[0] += dx * fx;
        v[1] += dy * fy;
        v[2] += dz * fz;
        v[3] += dx * fy;
        v[4] += dx * fz;
        v[5] += dy * fz;
    }

    Voigt6& operator+=(const Voigt6& o) noexcept {
        for (int k = 0; k < 6; ++k)
            v[k] += o.v[k];
        return *this;
    }

    /// @brief Accumulate into `frame.scalars.stress[0..5]`.
    void flush_to(SystemScalars& s) const noexcept {
        for (int k = 0; k < 6; ++k)
            s.stress[k] += v[k];
    }
};

} // namespace md
} // namespace cpu
} // namespace ATV
