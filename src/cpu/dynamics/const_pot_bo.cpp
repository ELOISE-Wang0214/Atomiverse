/// @file cpu/dynamics/const_pot_bo.cpp
/// @brief Dense BO constant-potential electrode charge solve.

#include "atomiverse/cpu/dynamics/const_pot_bo.h"

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/box.h"
#include "atomiverse/cpu/frame/frame.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace ATV {
namespace cpu {
namespace md {

namespace {

constexpr double kSqrt2 = 1.41421356237309504880;
constexpr double kSqrt2OverPi = 0.79788456080286535588; // sqrt(2/π)

/// Dense symmetric positive-definite solve via Cholesky (row-major).
void cholesky_solve(int n, std::vector<double>& A, std::vector<double>& x) {
    // A is overwritten with L (lower); x holds RHS then solution.
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double s = A[static_cast<std::size_t>(i) * n + j];
            for (int k = 0; k < j; ++k)
                s -= A[static_cast<std::size_t>(i) * n + k]
                     * A[static_cast<std::size_t>(j) * n + k];
            if (i == j) {
                if (s <= 0.0)
                    throw std::runtime_error(
                        "ConstPotBO: elastance matrix not SPD (Cholesky)");
                A[static_cast<std::size_t>(i) * n + j] = std::sqrt(s);
            } else {
                A[static_cast<std::size_t>(i) * n + j] =
                    s / A[static_cast<std::size_t>(j) * n + j];
            }
        }
    }
    // Forward L y = b
    for (int i = 0; i < n; ++i) {
        double s = x[static_cast<std::size_t>(i)];
        for (int k = 0; k < i; ++k)
            s -= A[static_cast<std::size_t>(i) * n + k]
                 * x[static_cast<std::size_t>(k)];
        x[static_cast<std::size_t>(i)] =
            s / A[static_cast<std::size_t>(i) * n + i];
    }
    // Back L^T x = y
    for (int i = n - 1; i >= 0; --i) {
        double s = x[static_cast<std::size_t>(i)];
        for (int k = i + 1; k < n; ++k)
            s -= A[static_cast<std::size_t>(k) * n + i]
                 * x[static_cast<std::size_t>(k)];
        x[static_cast<std::size_t>(i)] =
            s / A[static_cast<std::size_t>(i) * n + i];
    }
}

/// Invert SPD A (row-major) into inv (row-major) via n solves.
void invert_spd(int n, const std::vector<double>& A_in,
                std::vector<double>& inv) {
    inv.assign(static_cast<std::size_t>(n) * n, 0.0);
    std::vector<double> Awork(A_in);
    std::vector<double> col(static_cast<std::size_t>(n));
    for (int j = 0; j < n; ++j) {
        Awork = A_in;
        std::fill(col.begin(), col.end(), 0.0);
        col[static_cast<std::size_t>(j)] = 1.0;
        cholesky_solve(n, Awork, col);
        for (int i = 0; i < n; ++i)
            inv[static_cast<std::size_t>(i) * n + j] = col[static_cast<std::size_t>(i)];
    }
}

double min_image_r(const cpu::Box& box, double dx, double dy, double dz) {
    double rx = dx, ry = dy, rz = dz;
    box.min_image(rx, ry, rz);
    return std::sqrt(rx * rx + ry * ry + rz * rz);
}

} // namespace

ConstPotBO::ConstPotBO(double eta_in) : eta(eta_in) {
    if (!(eta > 0.0))
        throw std::invalid_argument("ConstPotBO: eta must be > 0 [1/Å]");
}

double ConstPotBO::capacitance() const {
    const double dpsi = std::abs(delta_psi());
    if (dpsi < 1e-30 || groups_.size() < 2) return 0.0;
    // |Q| on the last electrode group (opposite of first under neutrality).
    return std::abs(total_charge(groups_.size() - 1)) / dpsi;
}

ConstPotBO& ConstPotBO::add_electrode(int type, double psi) & {
    groups_.push_back(ElectrodeGroup{type, psi});
    return *this;
}

ConstPotBO ConstPotBO::add_electrode(int type, double psi) && {
    groups_.push_back(ElectrodeGroup{type, psi});
    return std::move(*this);
}

void ConstPotBO::build(const BuildContext& ctx) {
    if (groups_.empty())
        throw std::runtime_error(
            "ConstPotBO: add at least one electrode via add_electrode(type, Psi)");
    if (l_tf < 0.0)
        throw std::invalid_argument("ConstPotBO: l_tf must be >= 0 [Å]");
    if (l_tf > 0.0 && !(site_density > 0.0))
        throw std::invalid_argument(
            "ConstPotBO: site_density [1/Å³] must be > 0 when l_tf > 0");
    collect_electrodes_(ctx.frame);
    if (elec_atoms_.empty())
        throw std::runtime_error(
            "ConstPotBO: no atoms match registered electrode types");
    build_matrix_(ctx.frame);
    matrix_ready_ = true;
    total_Q_.assign(groups_.size(), 0.0);
}

void ConstPotBO::launch(const LaunchContext& ctx) {
    if (!matrix_ready_ || rebuild_matrix_every_step) {
        collect_electrodes_(ctx.frame);
        build_matrix_(ctx.frame);
        matrix_ready_ = true;
    }
    build_rhs_(ctx.frame);
    solve_();
    write_charges_(ctx.frame);
}

void ConstPotBO::destroy(const DestroyContext&) {
    matrix_ready_ = false;
    elec_atoms_.clear();
    A_.clear();
    Q_.clear();
}

double ConstPotBO::total_charge(std::size_t ie) const {
    if (ie >= total_Q_.size())
        throw std::out_of_range("ConstPotBO::total_charge");
    return total_Q_[ie];
}

double ConstPotBO::delta_psi() const {
    if (groups_.size() < 2) return 0.0;
    return groups_.back().psi - groups_.front().psi;
}

double ConstPotBO::neutrality_residual() const {
    double s = 0.0;
    for (double q : Q_) s += q;
    return std::abs(s);
}

void ConstPotBO::collect_electrodes_(const cpu::Frame& frame) {
    elec_atoms_.clear();
    elec_group_.clear();
    psi_vec_.clear();
    const int n = frame.nlocal();
    auto type = frame.span(kAtom_type);
    for (int i = 0; i < n; ++i) {
        for (std::size_t g = 0; g < groups_.size(); ++g) {
            if (static_cast<int>(type[i]) == groups_[g].type) {
                elec_atoms_.push_back(i);
                elec_group_.push_back(static_cast<int>(g));
                psi_vec_.push_back(groups_[g].psi);
                break;
            }
        }
    }
    const int m = static_cast<int>(elec_atoms_.size());
    Q_.assign(static_cast<std::size_t>(m), 0.0);
    b_.assign(static_cast<std::size_t>(m), 0.0);
    A_.assign(static_cast<std::size_t>(m) * m, 0.0);
}

double ConstPotBO::pair_A_(double r, double eta_i, double eta_j) const {
    if (r < 1e-14) {
        // Should not be called for i==j self; self uses eta * sqrt(2/π).
        return units::coulomb * eta_i * kSqrt2OverPi;
    }
    const double etaij = eta_i * eta_j / std::sqrt(eta_i * eta_i + eta_j * eta_j);
    return units::coulomb * std::erf(etaij * r) / r;
}

double ConstPotBO::pair_b_(double r, double q, double eta_i) const {
    if (r < 1e-14) return 0.0;
    return units::coulomb * q * std::erf(eta_i * r) / r;
}

void ConstPotBO::build_matrix_(const cpu::Frame& frame) {
    const int m = static_cast<int>(elec_atoms_.size());
    if (m == 0) return;
    auto x = frame.span(kAtom_x);
    auto y = frame.span(kAtom_y);
    auto z = frame.span(kAtom_z);
    const cpu::Box& box = frame.simbox;
    const double Lx = box.length(0);
    const double Ly = box.length(1);
    const int R = xy_image_range < 0 ? 0 : xy_image_range;

    std::fill(A_.begin(), A_.end(), 0.0);
    for (int a = 0; a < m; ++a) {
        const int ia = elec_atoms_[static_cast<std::size_t>(a)];
        // Self: k * η * sqrt(2/π)
        A_[static_cast<std::size_t>(a) * m + a] =
            units::coulomb * eta * kSqrt2OverPi;

        for (int b = a + 1; b < m; ++b) {
            const int ib = elec_atoms_[static_cast<std::size_t>(b)];
            double aij = 0.0;
            for (int nx = -R; nx <= R; ++nx) {
                for (int ny = -R; ny <= R; ++ny) {
                    // Primary cell: include once; for a==b images only nx,ny≠0
                    // handled in self loop below for a==b.
                    double dx = x[ia] - x[ib] + nx * Lx;
                    double dy = y[ia] - y[ib] + ny * Ly;
                    double dz = z[ia] - z[ib];
                    // When R==0 use min-image for general orthorhombic cells.
                    double r;
                    if (R == 0) {
                        r = min_image_r(box, x[ia] - x[ib], y[ia] - y[ib],
                                        z[ia] - z[ib]);
                    } else {
                        r = std::sqrt(dx * dx + dy * dy + dz * dz);
                    }
                    if (r < 1e-14) continue;
                    aij += pair_A_(r, eta, eta);
                    if (R == 0) break; // single min-image contrib
                }
                if (R == 0) break;
            }
            A_[static_cast<std::size_t>(a) * m + b] = aij;
            A_[static_cast<std::size_t>(b) * m + a] = aij;
        }

        // Self images (nx,ny) ≠ (0,0) for 2D lattice
        if (R > 0) {
            double self_img = 0.0;
            for (int nx = -R; nx <= R; ++nx) {
                for (int ny = -R; ny <= R; ++ny) {
                    if (nx == 0 && ny == 0) continue;
                    const double r =
                        std::sqrt((nx * Lx) * (nx * Lx) + (ny * Ly) * (ny * Ly));
                    // Same-site images use η/√2 (plate_cap / Gaussian product)
                    self_img +=
                        units::coulomb * std::erf((eta / kSqrt2) * r) / r;
                }
            }
            A_[static_cast<std::size_t>(a) * m + a] += self_img;
        }
    }

    // Thomas–Fermi metallicity (Scalfi et al. JCP 153, 174704 (2020)):
    // U_TF^(2) = (l_TF² d / (2 ε0)) Σ Q_i²  ⇒  A_ii += l_TF² d / ε0
    // with 1/ε0 = 4π k_e and k_e = units::coulomb [eV·Å/e²].
    if (l_tf > 0.0) {
        constexpr double kPi = 3.14159265358979323846;
        const double hardness =
            4.0 * kPi * units::coulomb * l_tf * l_tf * site_density;
        for (int a = 0; a < m; ++a)
            A_[static_cast<std::size_t>(a) * m + a] += hardness;
    }
}

void ConstPotBO::build_rhs_(const cpu::Frame& frame) {
    const int m = static_cast<int>(elec_atoms_.size());
    const int n = frame.nlocal();
    auto x = frame.span(kAtom_x);
    auto y = frame.span(kAtom_y);
    auto z = frame.span(kAtom_z);
    auto q = frame.span(kAtom_charge);
    auto type = frame.span(kAtom_type);
    const cpu::Box& box = frame.simbox;
    const double Lx = box.length(0);
    const double Ly = box.length(1);
    const int R = xy_image_range < 0 ? 0 : xy_image_range;

    // Mark electrode rows for skip
    std::vector<char> is_elec(static_cast<std::size_t>(n), 0);
    for (int ia : elec_atoms_) is_elec[static_cast<std::size_t>(ia)] = 1;

    std::fill(b_.begin(), b_.end(), 0.0);
    for (int a = 0; a < m; ++a) {
        const int ia = elec_atoms_[static_cast<std::size_t>(a)];
        double bi = 0.0;
        for (int j = 0; j < n; ++j) {
            if (is_elec[static_cast<std::size_t>(j)]) continue;
            (void)type; // electrolyte = non-electrode
            const double qj = static_cast<double>(q[j]);
            if (qj == 0.0) continue;
            for (int nx = -R; nx <= R; ++nx) {
                for (int ny = -R; ny <= R; ++ny) {
                    double r;
                    if (R == 0) {
                        r = min_image_r(box, x[ia] - x[j], y[ia] - y[j],
                                        z[ia] - z[j]);
                    } else {
                        const double dx = x[ia] - x[j] + nx * Lx;
                        const double dy = y[ia] - y[j] + ny * Ly;
                        const double dz = z[ia] - z[j];
                        r = std::sqrt(dx * dx + dy * dy + dz * dz);
                    }
                    bi += pair_b_(r, qj, eta);
                    if (R == 0) break;
                }
                if (R == 0) break;
            }
        }
        b_[static_cast<std::size_t>(a)] = bi;
    }
}

void ConstPotBO::solve_() {
    const int m = static_cast<int>(elec_atoms_.size());
    if (m == 0) return;

    // RHS: Ψ - b
    std::vector<double> rhs(static_cast<std::size_t>(m));
    for (int i = 0; i < m; ++i)
        rhs[static_cast<std::size_t>(i)] =
            psi_vec_[static_cast<std::size_t>(i)] - b_[static_cast<std::size_t>(i)];

    if (electroneutrality && m >= 1) {
        // Q = (A^{-1} - A^{-1} e e^T A^{-1} / (e^T A^{-1} e)) (Ψ - b)
        std::vector<double> inv;
        invert_spd(m, A_, inv);
        // u = A^{-1} e
        std::vector<double> u(static_cast<std::size_t>(m), 0.0);
        for (int i = 0; i < m; ++i) {
            double s = 0.0;
            for (int j = 0; j < m; ++j)
                s += inv[static_cast<std::size_t>(i) * m + j];
            u[static_cast<std::size_t>(i)] = s;
        }
        double eT_u = 0.0;
        for (double v : u) eT_u += v;
        if (std::abs(eT_u) < 1e-30)
            throw std::runtime_error(
                "ConstPotBO: singular neutrality projector");
        // v = A^{-1} rhs
        std::vector<double> v(static_cast<std::size_t>(m), 0.0);
        for (int i = 0; i < m; ++i) {
            double s = 0.0;
            for (int j = 0; j < m; ++j)
                s += inv[static_cast<std::size_t>(i) * m + j]
                     * rhs[static_cast<std::size_t>(j)];
            v[static_cast<std::size_t>(i)] = s;
        }
        double eT_v = 0.0;
        for (double t : v) eT_v += t;
        const double scale = eT_v / eT_u;
        for (int i = 0; i < m; ++i)
            Q_[static_cast<std::size_t>(i)] =
                v[static_cast<std::size_t>(i)] - scale * u[static_cast<std::size_t>(i)];
    } else {
        A_work_ = A_;
        Q_ = rhs;
        cholesky_solve(m, A_work_, Q_);
    }

    // Residual of unconstrained equation AQ -(Ψ-b) (will be ||λ|| for symm)
    residual_ = 0.0;
    for (int i = 0; i < m; ++i) {
        double Ax = 0.0;
        for (int j = 0; j < m; ++j)
            Ax += A_[static_cast<std::size_t>(i) * m + j]
                  * Q_[static_cast<std::size_t>(j)];
        residual_ = std::max(residual_,
                             std::abs(Ax - rhs[static_cast<std::size_t>(i)]));
    }

    // Per-group totals
    total_Q_.assign(groups_.size(), 0.0);
    for (int a = 0; a < m; ++a) {
        const int g = elec_group_[static_cast<std::size_t>(a)];
        total_Q_[static_cast<std::size_t>(g)] += Q_[static_cast<std::size_t>(a)];
    }
}

void ConstPotBO::write_charges_(cpu::Frame& frame) const {
    auto q = frame.span(kAtom_charge);
    // Zero electrode charges first (only electrode sites we own).
    for (std::size_t a = 0; a < elec_atoms_.size(); ++a)
        q[elec_atoms_[a]] = static_cast<StorageReal>(Q_[a]);
}

} // namespace md
} // namespace cpu
} // namespace ATV
