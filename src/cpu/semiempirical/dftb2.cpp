/// @file src/cpu/semiempirical/dftb2.cpp
/// Minimal H-only SCC-DFTB2 for framework wiring.
#include "atomiverse/cpu/semiempirical/dftb2.h"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "atomiverse/core/frame_schema.h"

namespace ATV {
namespace cpu {
namespace semiempirical {
namespace {

// Toy H–H SK-like tables (Å, eV) — order-of-magnitude textbook, not mio.
constexpr double kH0 = -6.5;    // eV
constexpr double kR0 = 0.74;    // Å
constexpr double kD = 0.5;      // Å
constexpr double kSWidth = 0.8; // Å
constexpr double kU_H = 13.6;   // Hubbard U [eV]
constexpr double kRepA = 0.5;   // eV
constexpr double kRepB = 2.0;   // 1/Å

double h0_ss(double r) { return kH0 * std::exp(-(r - kR0) / kD); }
double s_ss(double r) {
    return std::exp(-0.5 * (r / kSWidth) * (r / kSWidth));
}
double e_rep(double r) { return kRepA * std::exp(-kRepB * r); }

double gamma_ab(double r) {
    const double tau = std::sqrt(kU_H);
    return kU_H * std::exp(-tau * r / 27.2114);
}

} // namespace

void DFTB2::build(const BuildContext& /*ctx*/) {}

void DFTB2::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = f.nlocal();
    if (n < 1)
        throw std::invalid_argument("DFTB2: empty frame");
    if (n > 2)
        throw std::runtime_error(
            "DFTB2 skeleton: only N=1 or N=2 H atoms (textbook)");
    if (n % 2 != 0)
        throw std::runtime_error("DFTB2 skeleton: even number of H only");

    // Resolve spans once — Frame resolution budget forbids per-lookup loops.
    const auto type = f.span(kAtom_type);
    const StorageReal* x = f.span(kAtom_x);
    const StorageReal* y = f.span(kAtom_y);
    const StorageReal* z = f.span(kAtom_z);
    AccReal* fx = f.span(kAtom_fx);
    AccReal* fy = f.span(kAtom_fy);
    AccReal* fz = f.span(kAtom_fz);
    double& energy_meta = f.meta.get(kEnergy);

    for (int i = 0; i < n; ++i) {
        if (type[i] != 1)
            throw std::runtime_error(
                "DFTB2 skeleton: only H (Z=1) is supported");
    }

    auto dist = [&](int i, int j) {
        double dx = x[j] - x[i];
        double dy = y[j] - y[i];
        double dz = z[j] - z[i];
        f.simbox.min_image(dx, dy, dz);
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };

    auto idx = [n](int i, int j) {
        return static_cast<std::size_t>(i * n + j);
    };

    std::vector<double> S(static_cast<std::size_t>(n * n), 0.0);
    std::vector<double> H0(static_cast<std::size_t>(n * n), 0.0);
    for (int i = 0; i < n; ++i) {
        S[idx(i, i)] = 1.0;
        H0[idx(i, i)] = -0.5 * kU_H;
        for (int j = i + 1; j < n; ++j) {
            const double r = dist(i, j);
            S[idx(i, j)] = S[idx(j, i)] = s_ss(r);
            H0[idx(i, j)] = H0[idx(j, i)] = h0_ss(r);
        }
    }

    std::vector<double> dq(static_cast<std::size_t>(n), 0.0);
    double e_el = 0;
    const int n_occ = n / 2;

    for (int it = 0; it < max_scc; ++it) {
        std::vector<double> H = H0;
        for (int a = 0; a < n; ++a) {
            double shift = 0;
            for (int b = 0; b < n; ++b) {
                const double r = (a == b) ? 0.0 : dist(a, b);
                const double g = (a == b) ? kU_H : gamma_ab(r);
                shift += g * dq[static_cast<std::size_t>(b)];
            }
            H[idx(a, a)] += shift;
        }

        std::vector<double> eps(static_cast<std::size_t>(n), 0.0);
        std::vector<double> C(static_cast<std::size_t>(n * n), 0.0);
        if (n == 2) {
            const double a = H[idx(0, 0)], b = H[idx(1, 1)], c = H[idx(0, 1)];
            const double s = S[idx(0, 1)];
            const double A = 1.0 - s * s;
            const double B = -(a + b) + 2.0 * c * s;
            const double Cq = a * b - c * c;
            const double disc = std::max(0.0, B * B - 4.0 * A * Cq);
            const double e1 = (-B - std::sqrt(disc)) / (2.0 * A);
            const double e2 = (-B + std::sqrt(disc)) / (2.0 * A);
            eps[0] = e1;
            eps[1] = e2;
            auto fill_vec = [&](double e, int col) {
                double c1 = 1.0;
                double c2 = -(a - e) / (c - e * s + 1e-30);
                const double nrm =
                    std::sqrt(c1 * c1 + c2 * c2 + 2.0 * s * c1 * c2);
                c1 /= nrm;
                c2 /= nrm;
                C[idx(0, col)] = c1;
                C[idx(1, col)] = c2;
            };
            fill_vec(e1, 0);
            fill_vec(e2, 1);
        }

        std::vector<double> qnew(static_cast<std::size_t>(n), 0.0);
        for (int i = 0; i < n; ++i) {
            double qi = 0;
            for (int o = 0; o < n_occ; ++o) {
                const double ci = C[idx(i, o)];
                qi += 2.0 * ci * ci;
            }
            qnew[static_cast<std::size_t>(i)] = qi;
        }
        double max_dq = 0;
        for (int i = 0; i < n; ++i) {
            dq[static_cast<std::size_t>(i)] =
                qnew[static_cast<std::size_t>(i)] - 1.0;
            max_dq =
                std::max(max_dq, std::abs(dq[static_cast<std::size_t>(i)]));
        }
        e_el = 0;
        for (int o = 0; o < n_occ; ++o)
            e_el += 2.0 * eps[static_cast<std::size_t>(o)];
        double e_scc = 0;
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b) {
                const double r = (a == b) ? 0.0 : dist(a, b);
                const double g = (a == b) ? kU_H : gamma_ab(r);
                e_scc += 0.5 * g * dq[static_cast<std::size_t>(a)]
                         * dq[static_cast<std::size_t>(b)];
            }
        e_el += e_scc;
        if (max_dq < scc_tol)
            break;
    }

    double e_repulsive = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            e_repulsive += e_rep(dist(i, j));

    energy_ = e_el + e_repulsive;
    energy_meta += energy_;

    if (with_forces) {
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double dx = x[j] - x[i];
                double dy = y[j] - y[i];
                double dz = z[j] - z[i];
                f.simbox.min_image(dx, dy, dz);
                const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (r < 1e-12)
                    continue;
                const double dUdr = -kRepB * e_rep(r);
                const double fpair = -dUdr / r;
                fx[i] -= fpair * dx;
                fy[i] -= fpair * dy;
                fz[i] -= fpair * dz;
                fx[j] += fpair * dx;
                fy[j] += fpair * dy;
                fz[j] += fpair * dz;
            }
        }
    }
}

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
