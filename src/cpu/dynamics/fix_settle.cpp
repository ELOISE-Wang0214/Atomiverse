/// @file src/cpu/dynamics/fix_settle.cpp
/// @brief Analytic SETTLE for 3-site rigid water (Miyamoto–Kollman 1992).
///
/// Implementation follows the standard triangular reconstruction: conserve COM,
/// restore OH/HH distances in the plane of the three atoms after the
/// unconstrained move (GROMACS/OpenMM-style simplified SETTLE for v1).

#include "atomiverse/cpu/dynamics/fix_settle.h"

#include <cmath>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace md {

FixSettle::FixSettle(double dt, double d_oh, double d_hh)
    : dt_(dt), d_oh_(d_oh), d_hh_(d_hh) {
    if (!(dt > 0.0))
        throw std::invalid_argument("FixSettle: dt must be > 0");
    if (!(d_oh > 0.0) || !(d_hh > 0.0))
        throw std::invalid_argument("FixSettle: distances must be > 0");
    if (d_hh >= 2.0 * d_oh)
        throw std::invalid_argument("FixSettle: triangle inequality violated");
}

void FixSettle::add_water(TagInt tag_O, TagInt tag_H1, TagInt tag_H2) {
    if (tag_O == tag_H1 || tag_O == tag_H2 || tag_H1 == tag_H2)
        throw std::invalid_argument("FixSettle: water tags must be distinct");
    waters_.push_back({tag_O, tag_H1, tag_H2, -1, -1, -1});
}

void FixSettle::resolve_rows_(const Frame& f) {
    const int n = f.nlocal();
    const auto tags = f.span(kAtom_tag);
    for (auto& w : waters_) {
        w.io = w.ih1 = w.ih2 = -1;
        for (int i = 0; i < n; ++i) {
            if (tags[i] == w.o)
                w.io = i;
            if (tags[i] == w.h1)
                w.ih1 = i;
            if (tags[i] == w.h2)
                w.ih2 = i;
        }
        if (w.io < 0 || w.ih1 < 0 || w.ih2 < 0)
            throw std::runtime_error("FixSettle: missing water tags in frame");
    }
}

void FixSettle::build(const BuildContext& ctx) {
    if (waters_.empty())
        throw std::runtime_error("FixSettle: no waters registered");
    resolve_rows_(ctx.frame);
    n_atoms_ = ctx.frame.nlocal();
    x0_.assign(static_cast<size_t>(n_atoms_), 0.0);
    y0_.assign(static_cast<size_t>(n_atoms_), 0.0);
    z0_.assign(static_cast<size_t>(n_atoms_), 0.0);
}

void FixSettle::snapshot_(const Frame& f) {
    const int n = f.nlocal();
    if (n != n_atoms_) {
        n_atoms_ = n;
        x0_.resize(static_cast<size_t>(n));
        y0_.resize(static_cast<size_t>(n));
        z0_.resize(static_cast<size_t>(n));
    }
    for (int i = 0; i < n; ++i) {
        x0_[static_cast<size_t>(i)] = f.span(kAtom_x)[i];
        y0_[static_cast<size_t>(i)] = f.span(kAtom_y)[i];
        z0_[static_cast<size_t>(i)] = f.span(kAtom_z)[i];
    }
}

void FixSettle::settle_positions_(Frame& f) const {
    // Geometry constants for isosceles water triangle.
    // Height from O to HH midpoint: h = sqrt(d_oh² − (d_hh/2)²)
    const double half_hh = 0.5 * d_hh_;
    const double h = std::sqrt(d_oh_ * d_oh_ - half_hh * half_hh);
    const auto mass = f.span(kAtom_mass);

    for (const auto& w : waters_) {
        const int o = w.io, a = w.ih1, b = w.ih2;
        const double mo = mass[o], ma = mass[a], mb = mass[b];
        const double M = mo + ma + mb;
        // Unconstrained positions
        double xo = f.span(kAtom_x)[o], yo = f.span(kAtom_y)[o],
               zo = f.span(kAtom_z)[o];
        double xa = f.span(kAtom_x)[a], ya = f.span(kAtom_y)[a],
               za = f.span(kAtom_z)[a];
        double xb = f.span(kAtom_x)[b], yb = f.span(kAtom_y)[b],
               zb = f.span(kAtom_z)[b];
        // COM
        const double cx = (mo * xo + ma * xa + mb * xb) / M;
        const double cy = (mo * yo + ma * ya + mb * yb) / M;
        const double cz = (mo * zo + ma * za + mb * zb) / M;
        // Plane basis from unconstrained triangle
        double e1x = xa - xo, e1y = ya - yo, e1z = za - zo;
        double e2x = xb - xo, e2y = yb - yo, e2z = zb - zo;
        // Normal
        double nx = e1y * e2z - e1z * e2y;
        double ny = e1z * e2x - e1x * e2z;
        double nz = e1x * e2y - e1y * e2x;
        double nn = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (nn < 1e-18) {
            // Collinear fallback: use previous snapshot plane z-up.
            nx = 0;
            ny = 0;
            nz = 1;
            nn = 1;
        }
        nx /= nn;
        ny /= nn;
        nz /= nn;
        // In-plane axis: bisector of H-O-H from unconstrained
        double bx = 0.5 * (e1x + e2x), by = 0.5 * (e1y + e2y),
               bz = 0.5 * (e1z + e2z);
        double bn = std::sqrt(bx * bx + by * by + bz * bz);
        if (bn < 1e-18) {
            bx = 1;
            by = 0;
            bz = 0;
            bn = 1;
        }
        bx /= bn;
        by /= bn;
        bz /= bn;
        // Orthogonal in-plane
        double px = by * nz - bz * ny;
        double py = bz * nx - bx * nz;
        double pz = bx * ny - by * nx;
        // Place O, H1, H2 relative to COM with TIP3P masses.
        // Standard SETTLE mass-weighted positions along symmetry:
        // O along −bisector, H's at ± half_hh perpendicular.
        // Distance O–COM along bisector for masses mO, mH, mH:
        // r_O = −2 mH * h / M  (if height measured O to HH line)
        // Actually: place HH at y=0, O at y=h in body frame then shift COM to 0.
        const double y_O = h;
        const double y_H = 0.0;
        // Body coords
        const double body_Ox = 0, body_Oy = y_O;
        const double body_Ax = -half_hh, body_Ay = y_H;
        const double body_Bx = half_hh, body_By = y_H;
        const double com_x =
            (mo * body_Ox + ma * body_Ax + mb * body_Bx) / M;
        const double com_y =
            (mo * body_Oy + ma * body_Ay + mb * body_By) / M;
        auto place = [&](double bx_, double by_, int idx) {
            const double lx = (bx_ - com_x);
            const double ly = (by_ - com_y);
            // x = cx + lx * p + ly * b  (p ⊥ b in plane)
            f.span(kAtom_x)[idx] = cx + lx * px + ly * bx;
            f.span(kAtom_y)[idx] = cy + lx * py + ly * by;
            f.span(kAtom_z)[idx] = cz + lx * pz + ly * bz;
        };
        place(body_Ox, body_Oy, o);
        place(body_Ax, body_Ay, a);
        place(body_Bx, body_By, b);
    }
}

void FixSettle::settle_velocities_(Frame& f) const {
    // Coupled 3-bond velocity constraints: iterate RATTLE-style projections.
    const auto mass = f.span(kAtom_mass);
    for (int iter = 0; iter < 20; ++iter) {
        double max_rv = 0;
        for (const auto& w : waters_) {
            const int idxs[3][2] = {{w.io, w.ih1}, {w.io, w.ih2}, {w.ih1, w.ih2}};
            for (int c = 0; c < 3; ++c) {
                const int i = idxs[c][0], j = idxs[c][1];
                double dx = f.span(kAtom_x)[i] - f.span(kAtom_x)[j];
                double dy = f.span(kAtom_y)[i] - f.span(kAtom_y)[j];
                double dz = f.span(kAtom_z)[i] - f.span(kAtom_z)[j];
                double dvx = f.span(kAtom_vx)[i] - f.span(kAtom_vx)[j];
                double dvy = f.span(kAtom_vy)[i] - f.span(kAtom_vy)[j];
                double dvz = f.span(kAtom_vz)[i] - f.span(kAtom_vz)[j];
                const double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 < 1e-30)
                    continue;
                const double rv = dx * dvx + dy * dvy + dz * dvz;
                max_rv = std::max(max_rv, std::fabs(rv));
                const double mi = mass[i], mj = mass[j];
                const double inv_mu = 1.0 / mi + 1.0 / mj;
                const double lam = -rv / (r2 * inv_mu);
                f.span(kAtom_vx)[i] += lam * dx / mi;
                f.span(kAtom_vy)[i] += lam * dy / mi;
                f.span(kAtom_vz)[i] += lam * dz / mi;
                f.span(kAtom_vx)[j] -= lam * dx / mj;
                f.span(kAtom_vy)[j] -= lam * dy / mj;
                f.span(kAtom_vz)[j] -= lam * dz / mj;
            }
        }
        if (max_rv < 1e-12)
            break;
    }
}

void FixSettle::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    if (ctx.at(ATV::md::Phase::initial_integrate)) {
        snapshot_(f);
        return;
    }
    if (ctx.at(ATV::md::Phase::post_integrate)) {
        settle_positions_(f);
        return;
    }
    if (ctx.at(ATV::md::Phase::final_integrate)) {
        settle_velocities_(f);
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
