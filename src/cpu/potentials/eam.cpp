/// @file src/cpu/potentials/eam.cpp
/// @brief cpu::EAM + EAMPotential — analytic embedded-atom many-body style.

#include "atomiverse/cpu/potentials/eam.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/neighbor/pair_list.h"
#include "atomiverse/cpu/potentials/eam_data/eam_tables.h"

namespace ATV {
namespace cpu {

EAMPotential::EAMPotential(Element element, std::string_view name) {
    const int z = static_cast<int>(element);
    for (unsigned long k = 0; k < eam_data::kSuttonChenCount; ++k) {
        const eam_data::SuttonChenRow& row = eam_data::kSuttonChen[k];
        if (row.z == z && name == row.name) {
            n_ = row.n;
            m_ = row.m;
            eps_ = row.eps;
            c_ = row.c;
            a_ = row.a;
            return;
        }
    }
    // No fallback: an unsupported element or name is a programming error, not a
    // reason to silently substitute something plausible.
    throw std::out_of_range("EAMPotential: no parameterization '" + std::string(name)
                            + "' for Z=" + std::to_string(z));
}

EAM::EAM(EAMPotential potential, double cutoff)
    : pot_(potential), cutoff_(cutoff) {
    if (cutoff <= 0.0)
        throw std::invalid_argument("EAM: cutoff must be positive");
    // Sutton-Chen is a power law with no intrinsic range. Shift the pair term so
    // it reaches zero at the cutoff; without this a neighbour crossing the cutoff
    // steps the energy discontinuously and NVE stops conserving.
    phi_shift_ = pot_.phi(cutoff);
}

void EAM::build(const cpu::BuildContext& ctx) {
    if (ctx.runtime.mpi_size > 1) {
        throw std::runtime_error(
            "EAM: multi-rank is not supported yet. Pass 2 needs F'(rho_j) for "
            "ghost j, and a ghost's rho cannot be completed locally unless the "
            "halo reaches 2*cutoff; the engine exchanges exactly one forward halo "
            "(pre_neighbor) and one reverse halo (post_force, payload fx/fy/fz), "
            "so there is no mid-force slot for a derived per-atom quantity. "
            "Refusing to compute silently-wrong forces at rank boundaries.");
    }

    // rho_i is a sum over ALL neighbours of i. A half list gives each atom only
    // the half of them that landed in its own row, so ask for the full list.
    ctx.runtime.md.full_list_requested = true;

    const std::size_t n = static_cast<std::size_t>(ctx.frame.atoms.size());
    rho_.assign(n, 0.0);
    dF_.assign(n, 0.0);
}

void EAM::launch(const cpu::LaunchContext& ctx) {
    cpu::Frame& f = ctx.frame;
    const cpu::md::PairList* pl = ctx.runtime.md.pair_list;
    if (pl == nullptr)
        throw std::runtime_error("EAM: no neighbor list published");

    const int n = f.atoms.size();
    if (static_cast<int>(rho_.size()) < n) { // the frame grew under us
        rho_.assign(static_cast<std::size_t>(n), 0.0);
        dF_.assign(static_cast<std::size_t>(n), 0.0);
    }

    const double rc2 = cutoff_ * cutoff_;

    // ── Pass 1: rho_i = sum_{j != i} f(r_ij) ────────────────────────────
    // Every atom is an i here (not just pl->ilist), because pass 2 needs
    // F'(rho_j) for j on the other end of each pair. Single-rank, so every row
    // of the full list is complete.
    for (int i = 0; i < n; ++i) {
        double rho = 0.0;
        const int off = pl->firstneigh[i];
        const int nn = pl->numneigh[i];
        for (int jj = 0; jj < nn; ++jj) {
            const int j = pl->neighbors[off + jj];
            double dx = f.atoms.x[i] - f.atoms.x[j];
            double dy = f.atoms.y[i] - f.atoms.y[j];
            double dz = f.atoms.z[i] - f.atoms.z[j];
            f.box.min_image(dx, dy, dz);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= rc2) continue; // the list is a Verlet list; screen here
            rho += pot_.f(std::sqrt(r2));
        }
        rho_[static_cast<std::size_t>(i)] = rho;
        dF_[static_cast<std::size_t>(i)] = pot_.dF(rho);
    }

    // ── Pass 2: energy and forces ───────────────────────────────────────
    // E = sum_i F(rho_i) + 1/2 sum_i sum_{j!=i} phi(r_ij)
    //
    // dE/dr_ij collects the pair term once and BOTH embedding derivatives:
    //   dE/dr = phi'(r) + (F'(rho_i) + F'(rho_j)) * f'(r)
    // The second term is the many-body part — the force between i and j depends
    // on how densely embedded BOTH of them are.
    double e_total = 0.0;
    double v[6] = {0, 0, 0, 0, 0, 0};

    for (int i = 0; i < n; ++i) {
        e_total += pot_.F(rho_[static_cast<std::size_t>(i)]);

        const int off = pl->firstneigh[i];
        const int nn = pl->numneigh[i];
        for (int jj = 0; jj < nn; ++jj) {
            const int j = pl->neighbors[off + jj];
            double dx = f.atoms.x[i] - f.atoms.x[j];
            double dy = f.atoms.y[i] - f.atoms.y[j];
            double dz = f.atoms.z[i] - f.atoms.z[j];
            f.box.min_image(dx, dy, dz);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= rc2) continue;

            const double r = std::sqrt(r2);
            // The full list carries each pair twice (i->j and j->i), so every
            // per-pair quantity is halved here and the two visits sum to one.
            e_total += 0.5 * (pot_.phi(r) - phi_shift_);

            const double dE_dr =
                pot_.dphi(r)
                + (dF_[static_cast<std::size_t>(i)] + dF_[static_cast<std::size_t>(j)])
                      * pot_.df(r);
            // fpair = -(1/r) dE/dr, at FULL strength — not halved. The energy and
            // the virial are halved above/below because the full list visits each
            // pair twice and both visits accumulate into the same total. The force
            // is different: each visit writes only atom i's force, so the i-visit
            // must deposit the whole force that pair (i,j) exerts on i. (Halving
            // here was the first draft's bug; the FD oracle caught it immediately.)
            const double fpair = -dE_dr / r;

            f.atoms.fx[i] += fpair * dx;
            f.atoms.fy[i] += fpair * dy;
            f.atoms.fz[i] += fpair * dz;

            // Virial: each visit contributes its half of dx (x) f.
            v[0] += 0.5 * dx * fpair * dx;
            v[1] += 0.5 * dy * fpair * dy;
            v[2] += 0.5 * dz * fpair * dz;
            v[3] += 0.5 * dx * fpair * dy;
            v[4] += 0.5 * dx * fpair * dz;
            v[5] += 0.5 * dy * fpair * dz;
        }
    }

    f.scalars.energy += e_total;
    f.scalars.energy_pair += e_total;
    for (int k = 0; k < 6; ++k) f.scalars.stress[k] += v[k];
}

} // namespace cpu
} // namespace ATV
