/// @file src/cpu/md/fix_rattle.cpp
/// @brief cpu::md::FixRattle — iterative holonomic distance constraints.

#include "atomiverse/cpu/md/fix_rattle.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace md {

FixRattle::FixRattle(double dt, double tolerance, int max_iter)
    : dt_(dt), tol_(tolerance), max_iter_(max_iter) {
    if (dt <= 0.0) throw std::invalid_argument("FixRattle: dt must be positive");
    if (tolerance <= 0.0)
        throw std::invalid_argument("FixRattle: tolerance must be positive");
    if (max_iter < 1)
        throw std::invalid_argument("FixRattle: max_iter must be >= 1");
}

void FixRattle::constrain_bond(int bond_type, double length) {
    if (length <= 0.0)
        throw std::invalid_argument("FixRattle: constrained length must be positive");
    bond_rules_.push_back({bond_type, length});
}

void FixRattle::constrain_angle(int angle_type, double distance_13) {
    if (distance_13 <= 0.0)
        throw std::invalid_argument("FixRattle: constrained 1-3 distance must be "
                                    "positive");
    angle_rules_.push_back({angle_type, distance_13});
}

std::vector<FixRattle::Constraint>
FixRattle::discover_(const cpu::Frame& frame) const {
    std::vector<Constraint> out;
    const auto& a = frame.atoms;
    const int n = a.size();

    // Tag -> local row. Single-rank only (build() throws under MPI), so the frame's
    // own map is the whole story.
    const auto resolve = [&](std::int64_t tag) -> int {
        const auto it = a.map.find(static_cast<TagInt>(tag));
        return (it == a.map.end()) ? -1 : it->second;
    };

    const auto add = [&](int i, int j, double d) {
        Constraint c;
        c.i = i;
        c.j = j;
        c.d2 = d * d;
        c.inv_mi = 1.0 / a.mass[i];
        c.inv_mj = 1.0 / a.mass[j];
        out.push_back(c);
    };

    for (int i = 0; i < n; ++i) {
        for (int m = 0; m < a.num_bond[i]; ++m) {
            for (const Rule& r : bond_rules_) {
                if (a.bond_type(i, m) != r.type) continue;
                const int j = resolve(a.bond_atom(i, m));
                // A constraint whose partner cannot be resolved is a corrupted
                // molecule, not a small missing energy term. Bonded FORCE styles
                // silently drop such tuples and tally them; a constraint must not.
                if (j < 0)
                    throw std::runtime_error(
                        "FixRattle: bond partner tag "
                        + std::to_string(static_cast<long long>(a.bond_atom(i, m)))
                        + " of atom " + std::to_string(i)
                        + " does not resolve — refusing to leave a bond "
                          "unconstrained");
                add(i, j, r.distance);
            }
        }
        for (int m = 0; m < a.num_angle[i]; ++m) {
            for (const Rule& r : angle_rules_) {
                if (a.angle_type(i, m) != r.type) continue;
                // Atom i is the CENTRE of the stored angle; the two outer atoms are
                // the 1 and 3 of the triplet. Constraining their separation is what
                // freezes the angle.
                const int p = resolve(a.angle_atom(i, 2 * m + 0));
                const int q = resolve(a.angle_atom(i, 2 * m + 1));
                if (p < 0 || q < 0)
                    throw std::runtime_error(
                        "FixRattle: angle partner of atom " + std::to_string(i)
                        + " does not resolve — refusing to leave an angle "
                          "unconstrained");
                add(p, q, r.distance);
            }
        }
    }
    return out;
}

int FixRattle::n_constraints(const cpu::Frame& frame) const {
    return static_cast<int>(discover_(frame).size());
}

void FixRattle::build(const cpu::BuildContext& ctx) {
    if (ctx.runtime.mpi_size > 1)
        throw std::runtime_error(
            "FixRattle: multi-rank is not supported yet. The velocity stage needs "
            "the velocities of every atom in a cluster, and the forward halo carries "
            "POSITIONS ONLY — ghost velocities do not exist. Refusing to run with "
            "unconstrained molecules at rank boundaries.");

    constraints_ = discover_(ctx.frame);
    if (constraints_.empty())
        throw std::runtime_error(
            "FixRattle: no constraints matched the topology. Did you call "
            "constrain_bond()/constrain_angle() with the right type indices?");

    const std::size_t n = static_cast<std::size_t>(ctx.frame.atoms.size());
    ref_x_.assign(n, 0.0);
    ref_y_.assign(n, 0.0);
    ref_z_.assign(n, 0.0);
}

void FixRattle::launch(const cpu::LaunchContext& ctx) {
    cpu::Frame& f = ctx.frame;
    auto& a = f.atoms;

    // ── initial_integrate @ -500: snapshot the reference geometry ────────
    // Taken BEFORE any integrator moves anything. SHAKE/RATTLE corrects along the
    // reference bond vector, and VelocityVerlet fuses kick+drift and then wraps, so
    // by post_integrate the pre-move positions are unrecoverable from the frame.
    if (ctx.stage == static_cast<int>(ATV::md::Phase::initial_integrate)) {
        const int n = a.size();
        for (int i = 0; i < n; ++i) {
            ref_x_[static_cast<std::size_t>(i)] = a.x[i];
            ref_y_[static_cast<std::size_t>(i)] = a.y[i];
            ref_z_[static_cast<std::size_t>(i)] = a.z[i];
        }
        return;
    }

    // ── post_integrate: position stage ───────────────────────────────────
    // Project the drifted coordinates back onto sigma = |r|^2 - d^2 = 0, and feed
    // the correction back into the velocities (the atom did not travel as far as the
    // integrator thought it did).
    if (ctx.stage == static_cast<int>(ATV::md::Phase::post_integrate)) {
        const double inv_dt = 1.0 / dt_;
        int iter = 0;
        for (; iter < max_iter_; ++iter) {
            double worst = 0.0;
            for (const Constraint& c : constraints_) {
                double rx = a.x[c.i] - a.x[c.j];
                double ry = a.y[c.i] - a.y[c.j];
                double rz = a.z[c.i] - a.z[c.j];
                f.box.min_image(rx, ry, rz);
                const double r2 = rx * rx + ry * ry + rz * rz;
                const double sigma = r2 - c.d2;

                const double err = std::fabs(sigma) / c.d2;
                if (err > worst) worst = err;
                if (err <= tol_) continue;

                // Reference bond vector — the direction the Lagrange multiplier acts
                // along.
                double sx = ref_x_[static_cast<std::size_t>(c.i)]
                            - ref_x_[static_cast<std::size_t>(c.j)];
                double sy = ref_y_[static_cast<std::size_t>(c.i)]
                            - ref_y_[static_cast<std::size_t>(c.j)];
                double sz = ref_z_[static_cast<std::size_t>(c.i)]
                            - ref_z_[static_cast<std::size_t>(c.j)];
                f.box.min_image(sx, sy, sz);

                const double rs = rx * sx + ry * sy + rz * sz;
                const double denom = 2.0 * (c.inv_mi + c.inv_mj) * rs;
                if (std::fabs(denom) < 1e-12) continue; // degenerate; next sweep
                const double g = -sigma / denom;

                a.x[c.i] += g * c.inv_mi * sx;
                a.y[c.i] += g * c.inv_mi * sy;
                a.z[c.i] += g * c.inv_mi * sz;
                a.x[c.j] -= g * c.inv_mj * sx;
                a.y[c.j] -= g * c.inv_mj * sy;
                a.z[c.j] -= g * c.inv_mj * sz;

                // The constraint force also changed the momentum.
                a.vx[c.i] += g * c.inv_mi * sx * inv_dt;
                a.vy[c.i] += g * c.inv_mi * sy * inv_dt;
                a.vz[c.i] += g * c.inv_mi * sz * inv_dt;
                a.vx[c.j] -= g * c.inv_mj * sx * inv_dt;
                a.vy[c.j] -= g * c.inv_mj * sy * inv_dt;
                a.vz[c.j] -= g * c.inv_mj * sz * inv_dt;
            }
            if (worst <= tol_) break;
        }
        if (iter >= max_iter_)
            throw std::runtime_error(
                "FixRattle: the position stage did not converge in "
                + std::to_string(max_iter_)
                + " iterations. An unconverged constraint is a corrupted molecule, "
                  "not an approximation.");
        last_pos_iters_ = iter;
        return;
    }

    // ── final_integrate @ +200: velocity stage ───────────────────────────
    // A constraint on the coordinates implies one on the velocities:
    //   d(sigma)/dt = 2 (r_i - r_j).(v_i - v_j) = 0
    // i.e. the relative velocity must have NO component along the bond. Without this
    // the bond lengths stay right while the molecule quietly absorbs energy along
    // its own bonds.
    if (ctx.stage == static_cast<int>(ATV::md::Phase::final_integrate)) {
        int iter = 0;
        for (; iter < max_iter_; ++iter) {
            double worst = 0.0;
            for (const Constraint& c : constraints_) {
                double rx = a.x[c.i] - a.x[c.j];
                double ry = a.y[c.i] - a.y[c.j];
                double rz = a.z[c.i] - a.z[c.j];
                f.box.min_image(rx, ry, rz);

                const double vx = a.vx[c.i] - a.vx[c.j];
                const double vy = a.vy[c.i] - a.vy[c.j];
                const double vz = a.vz[c.i] - a.vz[c.j];

                const double rv = rx * vx + ry * vy + rz * vz;
                const double err = std::fabs(rv) / c.d2;
                if (err > worst) worst = err;
                if (err <= tol_) continue;

                const double k = -rv / ((c.inv_mi + c.inv_mj) * c.d2);

                a.vx[c.i] += k * c.inv_mi * rx;
                a.vy[c.i] += k * c.inv_mi * ry;
                a.vz[c.i] += k * c.inv_mi * rz;
                a.vx[c.j] -= k * c.inv_mj * rx;
                a.vy[c.j] -= k * c.inv_mj * ry;
                a.vz[c.j] -= k * c.inv_mj * rz;
            }
            if (worst <= tol_) break;
        }
        if (iter >= max_iter_)
            throw std::runtime_error(
                "FixRattle: the velocity stage did not converge in "
                + std::to_string(max_iter_) + " iterations.");
        last_vel_iters_ = iter;
        return;
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
