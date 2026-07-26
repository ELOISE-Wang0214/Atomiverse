/// @file src/cpu/multiscale/qmmm_polarizable.cpp
#include "atomiverse/cpu/multiscale/qmmm_polarizable.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"

namespace ATV {
namespace cpu {
namespace {

/// Thole exponential scale u = a * r / (α_i α_j)^{1/6}.
inline double thole_u(double r, double ai, double aj, double a) {
    if (!(a > 0.0) || !(ai > 0.0) || !(aj > 0.0))
        return 0.0;
    const double aij = std::pow(ai * aj, 1.0 / 6.0);
    return a * r / aij;
}

/// Damped dipole–dipole interaction: apply T·μ into E_out (Gaussian units).
void accumulate_T_mu(int n, const double* pos, const double* alpha,
                     const double* mu, double thole_a, double* e_out) {
    for (int i = 0; i < n; ++i) {
        const double* ri = pos + 3 * i;
        for (int j = 0; j < n; ++j) {
            if (i == j)
                continue;
            const double* rj = pos + 3 * j;
            const double* muj = mu + 3 * j;
            double dx = ri[0] - rj[0], dy = ri[1] - rj[1], dz = ri[2] - rj[2];
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (!(r2 > 0.0))
                continue;
            const double r = std::sqrt(r2);
            const double r3 = r2 * r;
            const double r5 = r3 * r2;
            // Undamped T_ab = (3 n_a n_b − δ_ab)/r³
            double f3 = 1.0 / r3;
            double f5 = 1.0 / r5;
            if (thole_a > 0.0) {
                const double u = thole_u(r, alpha[i], alpha[j], thole_a);
                // Exponential Thole: λ3 = 1 − exp(−u³), λ5 = 1 − (1+u³)exp(−u³)
                const double u3 = u * u * u;
                const double expu = std::exp(-u3);
                f3 *= (1.0 - expu);
                f5 *= (1.0 - (1.0 + u3) * expu);
            }
            const double mdotr = muj[0] * dx + muj[1] * dy + muj[2] * dz;
            e_out[3 * i + 0] += 3.0 * f5 * mdotr * dx - f3 * muj[0];
            e_out[3 * i + 1] += 3.0 * f5 * mdotr * dy - f3 * muj[1];
            e_out[3 * i + 2] += 3.0 * f5 * mdotr * dz - f3 * muj[2];
        }
    }
}

} // namespace

double solve_induced_dipoles(int n, const double* pos, const double* alpha,
                             const double* e_perm, double thole_a, double tol,
                             int max_iter, std::vector<double>& mu,
                             int* iters_used) {
    if (n < 0)
        throw std::invalid_argument("solve_induced_dipoles: n < 0");
    if (n == 0) {
        mu.clear();
        if (iters_used)
            *iters_used = 0;
        return 0.0;
    }
    if (!pos || !alpha || !e_perm)
        throw std::invalid_argument("solve_induced_dipoles: null pointer");
    if (!(tol > 0.0) || max_iter < 1)
        throw std::invalid_argument("solve_induced_dipoles: bad tol/max_iter");

    mu.assign(static_cast<std::size_t>(3 * n), 0.0);
    std::vector<double> e_tot(static_cast<std::size_t>(3 * n), 0.0);
    std::vector<double> mu_new(static_cast<std::size_t>(3 * n), 0.0);

    int it = 0;
    bool ok = false;
    for (; it < max_iter; ++it) {
        std::fill(e_tot.begin(), e_tot.end(), 0.0);
        for (int k = 0; k < 3 * n; ++k)
            e_tot[k] = e_perm[k];
        accumulate_T_mu(n, pos, alpha, mu.data(), thole_a, e_tot.data());

        double dmax = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int c = 0; c < 3; ++c) {
                const double m = alpha[i] * e_tot[3 * i + c];
                mu_new[3 * i + c] = m;
                dmax = std::max(dmax, std::fabs(m - mu[3 * i + c]));
            }
        }
        mu.swap(mu_new);
        if (dmax < tol) {
            ok = true;
            ++it; // report 1-based count of iterations performed
            break;
        }
    }
    if (iters_used)
        *iters_used = it;
    if (!ok)
        throw std::runtime_error(
            "solve_induced_dipoles: SCC did not converge within max_iter");

    // E_pol = −½ Σ μ · E_perm  (convert Gaussian → eV via coulomb)
    double e = 0.0;
    for (int i = 0; i < n; ++i) {
        e -= 0.5
             * (mu[3 * i + 0] * e_perm[3 * i + 0]
                + mu[3 * i + 1] * e_perm[3 * i + 1]
                + mu[3 * i + 2] * e_perm[3 * i + 2]);
    }
    return units::coulomb * e;
}

QmmmPolarizable::QmmmPolarizable(std::vector<TagInt> qm_tags, double mm_charge,
                                 double alpha, double thole_a, double epsilon,
                                 double sigma, double cutoff, double pol_tol,
                                 int pol_max_iter)
    : mm_q_(mm_charge), alpha_(alpha), thole_a_(thole_a), eps_(epsilon),
      sig_(sigma), rcut_(cutoff), pol_tol_(pol_tol),
      pol_max_iter_(pol_max_iter) {
    if (qm_tags.empty())
        throw std::invalid_argument("QmmmPolarizable: empty QM tags");
    if (!(alpha_ >= 0.0))
        throw std::invalid_argument("QmmmPolarizable: alpha must be >= 0");
    if (!(thole_a_ >= 0.0))
        throw std::invalid_argument("QmmmPolarizable: thole_a must be >= 0");
    if (!(eps_ > 0.0) || !(sig_ > 0.0) || !(cutoff > 0.0))
        throw std::invalid_argument("QmmmPolarizable: invalid LJ");
    if (!(pol_tol_ > 0.0) || pol_max_iter_ < 1)
        throw std::invalid_argument("QmmmPolarizable: invalid pol SCC params");
    for (TagInt t : qm_tags)
        qm_.insert(t);
}

void QmmmPolarizable::build(const BuildContext&) {}

void QmmmPolarizable::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = f.nlocal();
    const auto tags = f.span(kAtom_tag);
    const auto type = f.span(kAtom_type);
    const StorageReal* x = f.span(kAtom_x);
    const StorageReal* y = f.span(kAtom_y);
    const StorageReal* z = f.span(kAtom_z);
    double& E = f.meta.get(kEnergy);

    std::vector<int> qm_rows, mm_rows;
    for (int i = 0; i < n; ++i)
        (qm_.count(tags[i]) ? qm_rows : mm_rows).push_back(i);
    if (qm_rows.empty())
        throw std::runtime_error("QmmmPolarizable: no QM atoms");

    // QM sub-frame + permanent MM point-charge embedding
    cpu::Frame fqm(static_cast<int>(qm_rows.size()));
    for (std::size_t k = 0; k < qm_rows.size(); ++k) {
        const int i = qm_rows[k];
        fqm.span(kAtom_x)[k] = x[i];
        fqm.span(kAtom_y)[k] = y[i];
        fqm.span(kAtom_z)[k] = z[i];
        fqm.span(kAtom_type)[k] = type[i];
        fqm.span(kAtom_tag)[k] = static_cast<TagInt>(k + 1);
        fqm.span(kAtom_mass)[k] = 1.0;
        fqm.span(kAtom_fx)[k] = fqm.span(kAtom_fy)[k] = fqm.span(kAtom_fz)[k] =
            0;
    }
    fqm.box = f.box;

    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-8;
    rhf.max_iter = 80;
    rhf.compute_force = false;
    rhf.point_charges.clear();
    for (int i : mm_rows) {
        cpu::scf::RHF::PointCharge pc;
        pc.x = x[i];
        pc.y = y[i];
        pc.z = z[i];
        pc.q = mm_q_;
        rhf.point_charges.push_back(pc);
    }
    rhf.run(fqm);
    e_qm_ = fqm.meta.get(kEnergy);
    E += e_qm_;

    // Permanent field at MM sites [e/Å², Gaussian]: QM nuclei (Z) + other MM q.
    const int n_mm = static_cast<int>(mm_rows.size());
    std::vector<double> pos(static_cast<std::size_t>(3 * n_mm), 0.0);
    std::vector<double> alpha(static_cast<std::size_t>(n_mm), alpha_);
    std::vector<double> e_perm(static_cast<std::size_t>(3 * n_mm), 0.0);
    for (int a = 0; a < n_mm; ++a) {
        const int i = mm_rows[a];
        pos[3 * a + 0] = x[i];
        pos[3 * a + 1] = y[i];
        pos[3 * a + 2] = z[i];
    }
    for (int a = 0; a < n_mm; ++a) {
        double ex = 0, ey = 0, ez = 0;
        const double xi = pos[3 * a + 0], yi = pos[3 * a + 1],
                     zi = pos[3 * a + 2];
        // Field from QM nuclei treated as point charges Z (Gaussian units).
        // This is a textbook classical embedding approximation for E_perm;
        // electrons are not re-expanded here (v1).
        for (int iq : qm_rows) {
            const double q = static_cast<double>(type[iq]);
            double dx = xi - x[iq], dy = yi - y[iq], dz = zi - z[iq];
            f.box.min_image(dx, dy, dz);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (!(r2 > 0.0))
                continue;
            const double r3 = r2 * std::sqrt(r2);
            const double fac = q / r3;
            ex += fac * dx;
            ey += fac * dy;
            ez += fac * dz;
        }
        // Field from other permanent MM charges
        for (int b = 0; b < n_mm; ++b) {
            if (a == b)
                continue;
            double dx = xi - pos[3 * b + 0], dy = yi - pos[3 * b + 1],
                   dz = zi - pos[3 * b + 2];
            f.box.min_image(dx, dy, dz);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (!(r2 > 0.0))
                continue;
            const double r3 = r2 * std::sqrt(r2);
            const double fac = mm_q_ / r3;
            ex += fac * dx;
            ey += fac * dy;
            ez += fac * dz;
        }
        e_perm[3 * a + 0] = ex;
        e_perm[3 * a + 1] = ey;
        e_perm[3 * a + 2] = ez;
    }

    if (alpha_ > 0.0 && n_mm > 0) {
        e_pol_ = solve_induced_dipoles(
            n_mm, pos.data(), alpha.data(), e_perm.data(), thole_a_, pol_tol_,
            pol_max_iter_, last_mu_, &last_pol_iters_);
        last_pol_converged_ = true;
        E += e_pol_;
    } else {
        e_pol_ = 0;
        last_mu_.assign(static_cast<std::size_t>(3 * n_mm), 0.0);
        last_pol_iters_ = 0;
        last_pol_converged_ = true; // trivially (no SCC)
    }

    // LJ MM + cross
    e_mm_ = 0;
    const double rc2 = rcut_ * rcut_;
    const double sig2 = sig_ * sig_;
    const double s2c = sig2 / rc2;
    const double s6c = s2c * s2c * s2c;
    const double vshift = 4.0 * eps_ * (s6c * s6c - s6c);
    for (int i = 0; i < n; ++i) {
        const bool iqm = qm_.count(tags[i]) > 0;
        for (int j = i + 1; j < n; ++j) {
            const bool jqm = qm_.count(tags[j]) > 0;
            if (iqm && jqm)
                continue;
            double dx = x[j] - x[i], dy = y[j] - y[i], dz = z[j] - z[i];
            f.box.min_image(dx, dy, dz);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= rc2 || r2 == 0.0)
                continue;
            const double r2inv = 1.0 / r2;
            const double s2 = sig2 * r2inv;
            const double s6 = s2 * s2 * s2;
            e_mm_ += 4.0 * eps_ * (s6 * s6 - s6) - vshift;
        }
    }
    E += e_mm_;
}

} // namespace cpu
} // namespace ATV
