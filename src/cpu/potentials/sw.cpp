/// @file src/cpu/potentials/sw.cpp
/// @brief cpu::StillingerWeber — three-body force style.

#include "atomiverse/cpu/potentials/sw.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/neighbor/pair_list.h"
#include "atomiverse/cpu/potentials/sw_data/sw_tables.h"

namespace ATV {
namespace cpu {

SWPotential::SWPotential(Element element, std::string_view name) {
    const int z = static_cast<int>(element);
    for (unsigned long k = 0; k < sw_data::kStillingerWeberCount; ++k) {
        const sw_data::SWRow& row = sw_data::kStillingerWeber[k];
        if (row.z == z && name == row.name) {
            eps_ = row.eps;
            sigma_ = row.sigma;
            A_ = row.A;
            B_ = row.B;
            p_ = row.p;
            q_ = row.q;
            a_ = row.a;
            lambda_ = row.lambda;
            gamma_ = row.gamma;
            cos_theta0_ = row.cos_theta0;
            return;
        }
    }
    throw std::out_of_range("SWPotential: no parameterization '" + std::string(name)
                            + "' for Z=" + std::to_string(z));
}

StillingerWeber::StillingerWeber(SWPotential potential) : pot_(potential) {}

void StillingerWeber::build(const cpu::BuildContext& ctx) {
    if (ctx.runtime.mpi_size > 1) {
        throw std::runtime_error(
            "StillingerWeber: multi-rank is not supported yet. A triplet with its "
            "vertex on a ghost needs that ghost's own complete neighbour set, which "
            "requires the halo to reach 2*cutoff. Refusing to compute "
            "silently-wrong forces at rank boundaries.");
    }
    // The angular term enumerates every (j,k) pair of i's neighbours. A half list
    // shows i only half of them, so most triplets at i would be invisible.
    ctx.runtime.md.full_list_requested = true;
}

void StillingerWeber::launch(const cpu::LaunchContext& ctx) {
    cpu::Frame& f = ctx.frame;
    const cpu::md::PairList* pl = ctx.runtime.md.pair_list;
    if (pl == nullptr)
        throw std::runtime_error("StillingerWeber: no neighbor list published");

    const int n = f.atoms.size();
    const double rc = pot_.cutoff();
    const double rc2 = rc * rc;
    const double cos0 = pot_.cos_theta0();
    const double le = pot_.lambda_eps();

    double e_total = 0.0;
    double v[6] = {0, 0, 0, 0, 0, 0};

    for (int i = 0; i < n; ++i) {
        const int off = pl->firstneigh[i];
        const int nn = pl->numneigh[i];

        // ── Two-body ────────────────────────────────────────────────────
        // The full list visits each pair twice, so the energy and virial are
        // halved and the force is not: this visit writes only atom i.
        for (int jj = 0; jj < nn; ++jj) {
            const int j = pl->neighbors[off + jj];
            double dx = f.atoms.x[i] - f.atoms.x[j];
            double dy = f.atoms.y[i] - f.atoms.y[j];
            double dz = f.atoms.z[i] - f.atoms.z[j];
            f.box.min_image(dx, dy, dz);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= rc2) continue;

            const double r = std::sqrt(r2);
            e_total += 0.5 * pot_.phi2(r);
            const double fpair = -pot_.dphi2(r) / r;

            f.atoms.fx[i] += fpair * dx;
            f.atoms.fy[i] += fpair * dy;
            f.atoms.fz[i] += fpair * dz;

            v[0] += 0.5 * dx * fpair * dx;
            v[1] += 0.5 * dy * fpair * dy;
            v[2] += 0.5 * dz * fpair * dz;
            v[3] += 0.5 * dx * fpair * dy;
            v[4] += 0.5 * dx * fpair * dz;
            v[5] += 0.5 * dy * fpair * dz;
        }

        // ── Three-body ──────────────────────────────────────────────────
        // Vertex at i, legs to j and k. Each unordered {j,k} is enumerated once
        // (kk > jj), and each triplet is counted once per vertex — so no halving
        // here: unlike a pair, a triplet has a distinguished atom.
        //
        //   phi3 = lambda*eps * (cos(theta) - cos0)^2 * g(r1) * g(r2)
        //
        // The angle is what makes this irreducibly three-body: no function of r1
        // and r2 alone can express it.
        for (int jj = 0; jj < nn; ++jj) {
            const int j = pl->neighbors[off + jj];
            double ax = f.atoms.x[i] - f.atoms.x[j];
            double ay = f.atoms.y[i] - f.atoms.y[j];
            double az = f.atoms.z[i] - f.atoms.z[j];
            f.box.min_image(ax, ay, az);
            const double r1sq = ax * ax + ay * ay + az * az;
            if (r1sq >= rc2) continue;
            const double r1 = std::sqrt(r1sq);

            for (int kk = jj + 1; kk < nn; ++kk) {
                const int k = pl->neighbors[off + kk];
                double bx = f.atoms.x[i] - f.atoms.x[k];
                double by = f.atoms.y[i] - f.atoms.y[k];
                double bz = f.atoms.z[i] - f.atoms.z[k];
                f.box.min_image(bx, by, bz);
                const double r2sq = bx * bx + by * by + bz * bz;
                if (r2sq >= rc2) continue;
                const double r2 = std::sqrt(r2sq);

                const double g1 = pot_.g(r1), g2 = pot_.g(r2);
                const double dg1 = pot_.dg(r1), dg2 = pot_.dg(r2);

                const double cs = (ax * bx + ay * by + az * bz) / (r1 * r2);
                const double delta = cs - cos0;

                e_total += le * delta * delta * g1 * g2;

                // d(phi3)/d(cos) = 2*lambda*eps*delta*g1*g2
                const double dcos = 2.0 * le * delta * g1 * g2;
                // d(phi3)/dr1 and /dr2 (through the exponential legs)
                const double dE_dr1 = le * delta * delta * dg1 * g2;
                const double dE_dr2 = le * delta * delta * g1 * dg2;

                // d(cos)/d(a) and d(cos)/d(b), a = r_i - r_j, b = r_i - r_k
                const double inv_r1r2 = 1.0 / (r1 * r2);
                const double dcos_ax = bx * inv_r1r2 - cs * ax / r1sq;
                const double dcos_ay = by * inv_r1r2 - cs * ay / r1sq;
                const double dcos_az = bz * inv_r1r2 - cs * az / r1sq;
                const double dcos_bx = ax * inv_r1r2 - cs * bx / r2sq;
                const double dcos_by = ay * inv_r1r2 - cs * by / r2sq;
                const double dcos_bz = az * inv_r1r2 - cs * bz / r2sq;

                // Force on j: -d(phi3)/d(r_j). Since a = r_i - r_j,
                // d(a)/d(r_j) = -I and d(r1)/d(r_j) = -a/r1, so the two signs
                // cancel into a plain +.
                const double fjx = dE_dr1 * (ax / r1) + dcos * dcos_ax;
                const double fjy = dE_dr1 * (ay / r1) + dcos * dcos_ay;
                const double fjz = dE_dr1 * (az / r1) + dcos * dcos_az;

                const double fkx = dE_dr2 * (bx / r2) + dcos * dcos_bx;
                const double fky = dE_dr2 * (by / r2) + dcos * dcos_by;
                const double fkz = dE_dr2 * (bz / r2) + dcos * dcos_bz;

                // Translation invariance fixes the vertex force — it is not an
                // independent derivative, and computing it as one would be a
                // second chance to get a sign wrong.
                f.atoms.fx[j] += fjx;
                f.atoms.fy[j] += fjy;
                f.atoms.fz[j] += fjz;
                f.atoms.fx[k] += fkx;
                f.atoms.fy[k] += fky;
                f.atoms.fz[k] += fkz;
                f.atoms.fx[i] -= fjx + fkx;
                f.atoms.fy[i] -= fjy + fky;
                f.atoms.fz[i] -= fjz + fkz;

                // Virial of a three-body term, decomposed onto its two legs
                // (LAMMPS `ThreeBody::threebody` convention): with
                // delr1 = r_j - r_i = -a and delr2 = r_k - r_i = -b,
                //   v += delr1 (x) f_j + delr2 (x) f_k.
                v[0] += (-ax) * fjx + (-bx) * fkx;
                v[1] += (-ay) * fjy + (-by) * fky;
                v[2] += (-az) * fjz + (-bz) * fkz;
                v[3] += (-ax) * fjy + (-bx) * fky;
                v[4] += (-ax) * fjz + (-bx) * fkz;
                v[5] += (-ay) * fjz + (-by) * fkz;
            }
        }
    }

    f.scalars.energy += e_total;
    f.scalars.energy_pair += e_total;
    for (int c = 0; c < 6; ++c) f.scalars.stress[c] += v[c];
}

} // namespace cpu
} // namespace ATV
