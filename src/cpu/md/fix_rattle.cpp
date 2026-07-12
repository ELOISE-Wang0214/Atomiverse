/// @file src/cpu/md/fix_rattle.cpp
/// @brief cpu::md::FixRattle — iterative holonomic distance constraints.

#include "atomiverse/cpu/md/fix_rattle.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "atomiverse/cpu/comm/domain.h"
#include "atomiverse/cpu/frame/frame.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

namespace {

/// Owner rows of this rank: `dd->n_local` under decomposition, else the whole frame.
///
/// `frame.atoms.size()` is NOT the owner count under MPI — `partition_atoms` compacts
/// the owners to the front and leaves the tail as ghost slack, so the array still
/// spans the global atom count. Iterating it would read stale rows.
int owner_count(const cpu::Context& ctx, const cpu::Frame& frame) {
    const auto* dd = ctx.md.dd;
    if (dd != nullptr && ctx.mpi_size > 1) return dd->n_local;
    return frame.atoms.size();
}

} // namespace

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

std::vector<FixRattle::Spec> FixRattle::discover_local_(const cpu::Context& ctx,
                                                        const cpu::Frame& frame) const {
    std::vector<Spec> out;
    const auto& a = frame.atoms;
    const int n_own = owner_count(ctx, frame);

    // Only OWNED rows are scanned, and that is exactly right: LAMMPS bond ownership
    // gives every bond to its smaller-tag atom, and every angle to its centre atom, so
    // summing the per-rank discoveries over all ranks yields each constraint exactly
    // once — no double counting, nothing missed.
    for (int i = 0; i < n_own; ++i) {
        for (int m = 0; m < a.num_bond[i]; ++m)
            for (const Rule& r : bond_rules_)
                if (a.bond_type(i, m) == r.type)
                    out.push_back({a.tag[i], a.bond_atom(i, m), r.distance});

        for (int m = 0; m < a.num_angle[i]; ++m)
            for (const Rule& r : angle_rules_)
                if (a.angle_type(i, m) == r.type)
                    // Atom i is the CENTRE; the outer two atoms of the triplet are the
                    // 1 and 3. Freezing their separation freezes the angle.
                    out.push_back({a.angle_atom(i, 2 * m + 0),
                                   a.angle_atom(i, 2 * m + 1), r.distance});
    }
    return out;
}

int FixRattle::n_constraints(const cpu::Context& ctx, const cpu::Frame& frame) const {
    int n = static_cast<int>(discover_local_(ctx, frame).size());
#ifdef ATV_HAS_MPI
    // Ndof is GLOBAL, so the count must be too. Each constraint is discovered on
    // exactly one rank (the one owning its bond/angle), so a plain sum is correct.
    const auto* dd = ctx.md.dd;
    if (ctx.mpi_size > 1 && dd != nullptr && dd->cart_comm != MPI_COMM_NULL) {
        int local = n;
        MPI_Allreduce(&local, &n, 1, MPI_INT, MPI_SUM, dd->cart_comm);
    }
#endif
    return n;
}

void FixRattle::build(const cpu::BuildContext& ctx) {
    const cpu::Context& rt = ctx.runtime;
    const cpu::Frame& frame = ctx.frame;
    const auto& a = frame.atoms;
    const int n_own = owner_count(rt, frame);

    // ── The constraint table must be GLOBAL, not this rank's share ──────
    //
    // A water's two bonds and its angle all live on the oxygen (smallest tag / centre
    // atom). So the rank owning the O discovers the cluster — and a rank that owns
    // only an H discovers NOTHING, yet it is the rank responsible for correcting that
    // H. Each rank must therefore know every cluster it holds any atom of, which means
    // the table has to be shared. It is tiny (three doubles per constraint) and static,
    // so it is gathered once, here.
    std::vector<Spec> local = discover_local_(rt, frame);

    // Masses, likewise. Ghost rows carry x/y/z/tag/charge only — a ghost's `mass` is
    // never communicated, so the inverse masses have to be baked into the table while
    // we can still see every atom on its owner.
    std::unordered_map<TagInt, double> mass_of;
    for (int i = 0; i < n_own; ++i) mass_of[a.tag[i]] = a.mass[i];

    std::vector<Spec> global = local;

#ifdef ATV_HAS_MPI
    const auto* dd = rt.md.dd;
    if (rt.mpi_size > 1 && dd != nullptr && dd->cart_comm != MPI_COMM_NULL) {
        const int nranks = rt.mpi_size;

        // Gather the constraint specs (tag, tag, distance) as 3 doubles each.
        std::vector<double> sbuf;
        sbuf.reserve(local.size() * 3);
        for (const Spec& s : local) {
            sbuf.push_back(static_cast<double>(s.ti));
            sbuf.push_back(static_cast<double>(s.tj));
            sbuf.push_back(s.distance);
        }
        int my_n = static_cast<int>(sbuf.size());
        std::vector<int> counts(static_cast<std::size_t>(nranks), 0);
        MPI_Allgather(&my_n, 1, MPI_INT, counts.data(), 1, MPI_INT, dd->cart_comm);
        std::vector<int> displs(static_cast<std::size_t>(nranks), 0);
        int total = 0;
        for (int r = 0; r < nranks; ++r) {
            displs[static_cast<std::size_t>(r)] = total;
            total += counts[static_cast<std::size_t>(r)];
        }
        std::vector<double> rbuf(static_cast<std::size_t>(total), 0.0);
        MPI_Allgatherv(sbuf.data(), my_n, MPI_DOUBLE, rbuf.data(), counts.data(),
                       displs.data(), MPI_DOUBLE, dd->cart_comm);

        global.clear();
        for (int k = 0; k + 2 < total; k += 3)
            global.push_back({static_cast<TagInt>(rbuf[static_cast<std::size_t>(k)]),
                              static_cast<TagInt>(rbuf[static_cast<std::size_t>(k + 1)]),
                              rbuf[static_cast<std::size_t>(k + 2)]});

        // Same for (tag, mass).
        std::vector<double> msend;
        msend.reserve(static_cast<std::size_t>(n_own) * 2);
        for (int i = 0; i < n_own; ++i) {
            msend.push_back(static_cast<double>(a.tag[i]));
            msend.push_back(a.mass[i]);
        }
        int my_m = static_cast<int>(msend.size());
        MPI_Allgather(&my_m, 1, MPI_INT, counts.data(), 1, MPI_INT, dd->cart_comm);
        total = 0;
        for (int r = 0; r < nranks; ++r) {
            displs[static_cast<std::size_t>(r)] = total;
            total += counts[static_cast<std::size_t>(r)];
        }
        std::vector<double> mrecv(static_cast<std::size_t>(total), 0.0);
        MPI_Allgatherv(msend.data(), my_m, MPI_DOUBLE, mrecv.data(), counts.data(),
                       displs.data(), MPI_DOUBLE, dd->cart_comm);
        mass_of.clear();
        for (int k = 0; k + 1 < total; k += 2)
            mass_of[static_cast<TagInt>(mrecv[static_cast<std::size_t>(k)])] =
                mrecv[static_cast<std::size_t>(k + 1)];
    }
#endif

    if (global.empty())
        throw std::runtime_error(
            "FixRattle: no constraints matched the topology. Did you call "
            "constrain_bond()/constrain_angle() with the right type indices?");

    constraints_.clear();
    constraints_.reserve(global.size());
    for (const Spec& s : global) {
        const auto mi = mass_of.find(s.ti);
        const auto mj = mass_of.find(s.tj);
        if (mi == mass_of.end() || mj == mass_of.end())
            throw std::runtime_error(
                "FixRattle: a constrained atom has no mass in the global table — the "
                "topology references a tag that no rank owns");
        Constraint c;
        c.ti = s.ti;
        c.tj = s.tj;
        c.d2 = s.distance * s.distance;
        c.inv_mi = 1.0 / mi->second;
        c.inv_mj = 1.0 / mj->second;
        constraints_.push_back(c);
    }

    const std::size_t cap = static_cast<std::size_t>(frame.atoms.nmax);
    ref_x_.assign(cap, 0.0);
    ref_y_.assign(cap, 0.0);
    ref_z_.assign(cap, 0.0);
    rows_.assign(constraints_.size() * 2, -1);
}

void FixRattle::resolve_rows_(const cpu::LaunchContext& ctx) {
    const auto& a = ctx.frame.atoms;
    const auto* dd = ctx.runtime.md.dd;
    const bool multi = (dd != nullptr && ctx.runtime.mpi_size > 1);

    for (std::size_t k = 0; k < constraints_.size(); ++k) {
        const Constraint& c = constraints_[k];
        int i = -1, j = -1;
        if (multi) {
            const auto it = dd->tag_to_local.find(static_cast<std::int64_t>(c.ti));
            const auto jt = dd->tag_to_local.find(static_cast<std::int64_t>(c.tj));
            if (it != dd->tag_to_local.end()) i = it->second;
            if (jt != dd->tag_to_local.end()) j = jt->second;
        } else {
            const auto it = a.map.find(c.ti);
            const auto jt = a.map.find(c.tj);
            if (it != a.map.end()) i = it->second;
            if (jt != a.map.end()) j = jt->second;
        }

        // Neither atom here: the cluster lives entirely on other ranks. Skip it.
        if (i < 0 && j < 0) {
            rows_[2 * k] = -1;
            rows_[2 * k + 1] = -1;
            continue;
        }
        // One atom here and the other nowhere: the halo is too thin to hold the whole
        // cluster, so this rank cannot solve it — and neither can anyone else, because
        // the solve must be redundant across every rank holding a piece. A dropped
        // constraint is a molecule that is not rigid: throw rather than approximate.
        if (i < 0 || j < 0)
            throw std::runtime_error(
                "FixRattle: a constrained pair spans ranks but its partner is not even "
                "a ghost here. The halo must cover the whole constrained cluster — "
                "widen the cutoff/skin, use fewer ranks, or use a bigger box.");

        rows_[2 * k] = i;
        rows_[2 * k + 1] = j;
    }
}

void FixRattle::launch(const cpu::LaunchContext& ctx) {
    cpu::Frame& f = ctx.frame;
    auto& a = f.atoms;
    const auto* dd = ctx.runtime.md.dd;
    const bool multi = (dd != nullptr && ctx.runtime.mpi_size > 1);
    const int n_rows = multi ? (dd->n_local + dd->n_ghost) : a.size();

    // ── initial_integrate @ -500: snapshot the reference geometry ────────
    // Taken BEFORE any integrator moves anything. SHAKE/RATTLE corrects along the
    // reference bond vector, and VelocityVerlet fuses kick+drift and then wraps, so by
    // post_integrate the pre-move positions are unrecoverable from the frame.
    if (ctx.stage == static_cast<int>(ATV::md::Phase::initial_integrate)) {
        resolve_rows_(ctx);
        if (static_cast<int>(ref_x_.size()) < n_rows) {
            ref_x_.assign(static_cast<std::size_t>(n_rows), 0.0);
            ref_y_.assign(static_cast<std::size_t>(n_rows), 0.0);
            ref_z_.assign(static_cast<std::size_t>(n_rows), 0.0);
        }
        for (int i = 0; i < n_rows; ++i) {
            ref_x_[static_cast<std::size_t>(i)] = a.x[i];
            ref_y_[static_cast<std::size_t>(i)] = a.y[i];
            ref_z_[static_cast<std::size_t>(i)] = a.z[i];
        }
        return;
    }

    // Every rank holding any atom of a cluster solves that cluster REDUNDANTLY, from
    // identical inputs (a HaloRefresh has just re-sent the ghosts' x and v), and so
    // arrives at identical corrections. It applies them to all its rows; the ghost-row
    // writes are local scratch, discarded at the next refresh, while each owner ends up
    // with exactly the correction its owning rank computed. No extra reverse exchange
    // is needed — this is LAMMPS's replicated-cluster route.
    if (ctx.stage == static_cast<int>(ATV::md::Phase::post_integrate)) {
        const double inv_dt = 1.0 / dt_;
        int iter = 0;
        for (; iter < max_iter_; ++iter) {
            double worst = 0.0;
            for (std::size_t k = 0; k < constraints_.size(); ++k) {
                const int ci = rows_[2 * k], cj = rows_[2 * k + 1];
                if (ci < 0) continue;
                const Constraint& c = constraints_[k];

                double rx = a.x[ci] - a.x[cj];
                double ry = a.y[ci] - a.y[cj];
                double rz = a.z[ci] - a.z[cj];
                f.box.min_image(rx, ry, rz);
                const double sigma = rx * rx + ry * ry + rz * rz - c.d2;

                const double err = std::fabs(sigma) / c.d2;
                if (err > worst) worst = err;
                if (err <= tol_) continue;

                double sx = ref_x_[static_cast<std::size_t>(ci)]
                            - ref_x_[static_cast<std::size_t>(cj)];
                double sy = ref_y_[static_cast<std::size_t>(ci)]
                            - ref_y_[static_cast<std::size_t>(cj)];
                double sz = ref_z_[static_cast<std::size_t>(ci)]
                            - ref_z_[static_cast<std::size_t>(cj)];
                f.box.min_image(sx, sy, sz);

                const double rs = rx * sx + ry * sy + rz * sz;
                const double denom = 2.0 * (c.inv_mi + c.inv_mj) * rs;
                if (std::fabs(denom) < 1e-12) continue;
                const double g = -sigma / denom;

                a.x[ci] += g * c.inv_mi * sx;
                a.y[ci] += g * c.inv_mi * sy;
                a.z[ci] += g * c.inv_mi * sz;
                a.x[cj] -= g * c.inv_mj * sx;
                a.y[cj] -= g * c.inv_mj * sy;
                a.z[cj] -= g * c.inv_mj * sz;

                // The constraint force also changed the momentum.
                a.vx[ci] += g * c.inv_mi * sx * inv_dt;
                a.vy[ci] += g * c.inv_mi * sy * inv_dt;
                a.vz[ci] += g * c.inv_mi * sz * inv_dt;
                a.vx[cj] -= g * c.inv_mj * sx * inv_dt;
                a.vy[cj] -= g * c.inv_mj * sy * inv_dt;
                a.vz[cj] -= g * c.inv_mj * sz * inv_dt;
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
    //   d(sigma)/dt = 2 (r_i - r_j).(v_i - v_j) = 0
    // The relative velocity must have NO component along the bond. Skip this and the
    // bond lengths stay right while the molecule quietly absorbs energy along its own
    // bonds.
    if (ctx.stage == static_cast<int>(ATV::md::Phase::final_integrate)) {
        // Re-resolve. `HaloSend` runs at pre_neighbor (phase 3), between the position
        // stage (1) and this one (8), and a full rebuild RESHUFFLES the ghost slots —
        // so the rows resolved at initial_integrate are stale by now, and the velocity
        // stage was silently constraining the wrong atoms. It showed up as a
        // non-converging solve, which is the good failure mode; had the shuffle been
        // milder it would simply have produced wrong dynamics.
        resolve_rows_(ctx);
        int iter = 0;
        for (; iter < max_iter_; ++iter) {
            double worst = 0.0;
            for (std::size_t k = 0; k < constraints_.size(); ++k) {
                const int ci = rows_[2 * k], cj = rows_[2 * k + 1];
                if (ci < 0) continue;
                const Constraint& c = constraints_[k];

                double rx = a.x[ci] - a.x[cj];
                double ry = a.y[ci] - a.y[cj];
                double rz = a.z[ci] - a.z[cj];
                f.box.min_image(rx, ry, rz);

                const double vx = a.vx[ci] - a.vx[cj];
                const double vy = a.vy[ci] - a.vy[cj];
                const double vz = a.vz[ci] - a.vz[cj];

                const double rv = rx * vx + ry * vy + rz * vz;
                const double err = std::fabs(rv) / c.d2;
                if (err > worst) worst = err;
                if (err <= tol_) continue;

                const double k_lag = -rv / ((c.inv_mi + c.inv_mj) * c.d2);

                a.vx[ci] += k_lag * c.inv_mi * rx;
                a.vy[ci] += k_lag * c.inv_mi * ry;
                a.vz[ci] += k_lag * c.inv_mi * rz;
                a.vx[cj] -= k_lag * c.inv_mj * rx;
                a.vy[cj] -= k_lag * c.inv_mj * ry;
                a.vz[cj] -= k_lag * c.inv_mj * rz;
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
