/// @file tests/cpu/mc/mc_test_util.h
/// @brief Test-scope helpers for the MC physics-validation gates
///        (mc-sampling-03-validate). Header-local, no ATV:: public surface:
///        1-D quadrature, shifted / pure-cut LJ energies, a truncated
///        full-frame energy re-evaluation, and block standard-error.
#pragma once

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"

namespace mctest {

/// @brief FCC lattice of 4·nc³ atoms at reduced density ρ* in a cubic box
///        (σ = 1 units). @param nc cells/side. @param rho ρ*. @param n_out
///        receives the atom count. @return the frame (box set, mass 1).
inline ATV::cpu::Frame fcc(int nc, double rho, int& n_out) {
    const int n = 4 * nc * nc * nc;
    const double a = std::cbrt(4.0 / rho); // lattice constant [σ]
    const double L = nc * a;
    ATV::cpu::Frame f(n);
    const double bx[4] = {0.0, 0.5, 0.5, 0.0};
    const double by[4] = {0.0, 0.5, 0.0, 0.5};
    const double bz[4] = {0.0, 0.0, 0.5, 0.5};
    int idx = 0;
    for (int ix = 0; ix < nc; ++ix)
        for (int iy = 0; iy < nc; ++iy)
            for (int iz = 0; iz < nc; ++iz)
                for (int b = 0; b < 4; ++b) {
                    f.span(ATV::kAtom_x)[idx] = (ix + bx[b]) * a;
                    f.span(ATV::kAtom_y)[idx] = (iy + by[b]) * a;
                    f.span(ATV::kAtom_z)[idx] = (iz + bz[b]) * a;
                    f.span(ATV::kAtom_mass)[idx] = 1.0;
                    ++idx;
                }
    f.simbox = ATV::cpu::Box::cubic(L);
    f.meta.get(ATV::kEnergy) = 0.0;
    n_out = n;
    return f;
}

/// @brief Shifted 12-6 LJ pair energy [eV] — matches cpu::PairLJCut
///        (V(r) = 4ε[(σ/r)¹²−(σ/r)⁶] − V(rc) for r<rc, else 0).
inline double lj_shifted(double r, double eps, double sig, double rc) {
    auto raw = [&](double x) {
        const double s6 = std::pow(sig / x, 6.0);
        return 4.0 * eps * (s6 * s6 - s6);
    };
    return (r < rc) ? raw(r) - raw(rc) : 0.0;
}

/// @brief Pure-cut (unshifted) truncated 12-6 LJ pair energy [eV] — the form
///        the NIST SRSW reference is defined against.
inline double lj_cut(double r, double eps, double sig, double rc) {
    if (r >= rc)
        return 0.0;
    const double s6 = std::pow(sig / r, 6.0);
    return 4.0 * eps * (s6 * s6 - s6);
}

/// @brief Trapezoid ∫_a^b f(r) dr over `n` intervals.
template <class F>
double integrate(F f, double a, double b, int n) {
    const double h = (b - a) / n;
    double s = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; ++i)
        s += f(a + i * h);
    return s * h;
}

/// @brief Conditional ⟨V⟩ [eV] over separations in [a,b] with the radial
///        weight r²·exp(−βV(r)) — the analytic reference for the sampled ⟨V⟩.
inline double mean_V_quadrature(double beta, double eps, double sig, double rc,
                                double a, double b, int n = 40000) {
    auto w = [&](double r) {
        return r * r * std::exp(-beta * lj_shifted(r, eps, sig, rc));
    };
    const double num = integrate(
        [&](double r) { return lj_shifted(r, eps, sig, rc) * w(r); }, a, b, n);
    const double den = integrate(w, a, b, n);
    return num / den;
}

/// @brief Unnormalized radial density r²·exp(−βV_shifted(r)) at `r`.
inline double radial_density(double r, double beta, double eps, double sig,
                             double rc) {
    return r * r * std::exp(-beta * lj_shifted(r, eps, sig, rc));
}

/// @brief Block-averaged mean and standard error of a series.
/// @param x Samples. @param nblocks Number of contiguous blocks.
/// @return {mean, standard error of the mean via block averaging}.
inline std::pair<double, double> block_mean_se(const std::vector<double>& x,
                                               int nblocks) {
    const std::size_t n = x.size();
    const std::size_t bs = n / static_cast<std::size_t>(nblocks);
    std::vector<double> bm;
    bm.reserve(nblocks);
    for (int b = 0; b < nblocks; ++b) {
        double s = 0.0;
        for (std::size_t i = b * bs; i < (b + 1) * bs; ++i)
            s += x[i];
        bm.push_back(s / bs);
    }
    double mean = 0.0;
    for (double v : bm)
        mean += v;
    mean /= nblocks;
    double var = 0.0;
    for (double v : bm)
        var += (v - mean) * (v - mean);
    var /= (nblocks - 1);
    return {mean, std::sqrt(var / nblocks)};
}

/// @brief Total pure-cut (unshifted) truncated LJ configurational energy [eV]
///        of a snapshot, summed over all minimum-image pairs with r < rc.
///        O(N²) — for the small NIST cell, evaluated only on stored snapshots.
inline double truncated_lj_energy(const ATV::cpu::Frame& f, double eps,
                                  double sig, double rc) {
    const int n = f.nlocal();
    double e = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dx = static_cast<double>(f.span(ATV::kAtom_x)[i]
                                            - f.span(ATV::kAtom_x)[j]);
            double dy = static_cast<double>(f.span(ATV::kAtom_y)[i]
                                            - f.span(ATV::kAtom_y)[j]);
            double dz = static_cast<double>(f.span(ATV::kAtom_z)[i]
                                            - f.span(ATV::kAtom_z)[j]);
            f.simbox.min_image(dx, dy, dz);
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            e += lj_cut(r, eps, sig, rc);
        }
    }
    return e;
}

} // namespace mctest
