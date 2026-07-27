/// @file regressions/md/maze-symplectic-charge-dynamics/reg_maze-symplectic-charge-dynamics.cpp
/// @brief E2E: true mass-zero constrained (SHAKE-like) MaZe charge dynamics —
///        time-reversibility + energy conservation vs the CG+ASPC twin.
///
/// Public API only: `ConstPotMaZe` + `PairGaussianCoul` components driven by
/// an explicit NVE velocity-Verlet loop (the elastance/RHS build is a full
/// min-image sum with no cutoff, so the pair twin must also see every
/// min-image pair — a driver NeighborList cannot cover that; the explicit
/// loop also makes the velocity-flip palindrome exact).
///
/// Literature pin (no MetalWalls binary): Coretti et al., J. Chem. Phys.
/// 152, 194701 (2020) — see REFERENCE.md. Exits non-zero on failure.

#include <cmath>
#include <cstdio>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/dynamics/const_pot_maze.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/potentials/pair_gaussian_coul.h"

using namespace ATV;

namespace {

constexpr double kLxy = 6.0;    // [Å]
constexpr double kLz = 20.0;    // [Å]
constexpr double kEta = 2.0;    // [1/Å]
constexpr double kPsiBot = 0.0; // [V]  type 1, z = 4 Å
constexpr double kPsiTop = 1.0; // [V]  type 2, z = 16 Å
constexpr double kCut = 12.0;   // [Å]  > max min-image distance (10.9 Å)
constexpr double kDt = 0.25;    // [fs]
constexpr double kMion = 40.0;  // [amu]
constexpr int kSteps = 400;     // forward steps (then 400 back)
constexpr int kIon0 = 8;        // first ion row (walls are rows 0–7)

/// 2×2 wall sites per plate (types 1/2, static) + two point ions (type 3).
cpu::Frame make_capacitor() {
    const int N = 10;
    cpu::Frame f(N);
    int a = 0;
    for (int wall = 0; wall < 2; ++wall) {
        const double z = (wall == 0) ? 4.0 : 16.0;
        const int typ = (wall == 0) ? 1 : 2;
        for (int iy = 0; iy < 2; ++iy)
            for (int ix = 0; ix < 2; ++ix) {
                f.span(kAtom_x)[a] = (ix + 0.5) * (kLxy / 2);
                f.span(kAtom_y)[a] = (iy + 0.5) * (kLxy / 2);
                f.span(kAtom_z)[a] = z;
                f.span(kAtom_type)[a] = typ;
                f.span(kAtom_charge)[a] = 0.0;
                f.span(kAtom_mass)[a] = 196.97;
                f.span(kAtom_tag)[a] = a + 1;
                ++a;
            }
    }
    // Point ions ±0.5 e, laterally separated (min-image xy distance 2.97 Å;
    // offsets chosen so no ion–wall pair sits near the L/2 min-image kink).
    const double ion_x[2] = {2.1, 4.2};
    const double ion_y[2] = {2.1, 4.2};
    const double ion_z[2] = {8.5, 11.5};
    const double ion_q[2] = {0.5, -0.5};
    for (int i = 0; i < 2; ++i) {
        f.span(kAtom_x)[a] = ion_x[i];
        f.span(kAtom_y)[a] = ion_y[i];
        f.span(kAtom_z)[a] = ion_z[i];
        f.span(kAtom_type)[a] = 3;
        f.span(kAtom_charge)[a] = ion_q[i];
        f.span(kAtom_mass)[a] = kMion;
        f.span(kAtom_tag)[a] = a + 1;
        ++a;
    }
    f.box = cpu::Box::orthorhombic(kLxy, kLxy, kLz);
    return f;
}

struct RunResult {
    double drift = 0.0;       ///< max |E(t) − E_ref| after settling [eV]
    double ret_pos = 0.0;     ///< max ion |Δr| vs start after palindrome [Å]
    double ret_q = 0.0;       ///< max electrode |ΔQ| vs start [e]
    double max_resid = 0.0;   ///< max constraint residual over run [V]
    double max_neut = 0.0;    ///< max |ΣQ| over run [e]
    double min_sep = 1e300;   ///< min ion–ion min-image distance [Å]
    double q_top_first = 0.0; ///< total Q on the +1 V plate, first solve [e]
};

/// Forward kSteps, velocity flip, backward kSteps. NVE, static walls.
RunResult run_case(bool symplectic, double foothold_cg_tol) {
    cpu::Frame f = make_capacitor();

    cpu::md::ConstPotMaZe maze(kEta);
    maze.add_electrode(1, kPsiBot);
    maze.add_electrode(2, kPsiTop);
    maze.symplectic = symplectic;
    if (!symplectic) {
        maze.solver().cg_tol = foothold_cg_tol; // loose adiabatic twin
        maze.solver().predictor_order = 2;      // Kolafa ASPC warm start
    }

    // Per-type Gaussian widths for the pair energy/forces: walls Gaussian
    // (η = kEta), ions point charges — identical kernels to the elastance.
    std::vector<double> etas = {0.0, kEta, kEta, 0.0};
    cpu::PairGaussianCoul pg(kCut, etas, /*include_self=*/true);

    cpu::Context ctx;
    maze.build({ctx, f});
    pg.build({ctx, f});

    const double x0[2] = {f.span(kAtom_x)[kIon0], f.span(kAtom_x)[kIon0 + 1]};
    const double y0[2] = {f.span(kAtom_y)[kIon0], f.span(kAtom_y)[kIon0 + 1]};
    const double z0[2] = {f.span(kAtom_z)[kIon0], f.span(kAtom_z)[kIon0 + 1]};
    double vel[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    RunResult res;
    auto q_plate = [&](int typ) {
        double s = 0.0;
        for (int i = 0; i < f.nlocal(); ++i)
            if (f.span(kAtom_type)[i] == typ)
                s += f.span(kAtom_charge)[i];
        return s;
    };
    // Charge solve (pre_force) then pair energy/forces at current positions.
    int step = 0;
    auto eval = [&]() {
        for (int i = 0; i < f.nlocal(); ++i)
            f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] =
                0.0;
        f.meta.get(kEnergy) = 0.0;
        maze.launch({ctx, f, 0, step});
        pg.launch({ctx, f, 0, step});
        ++step;
    };
    auto conserved = [&]() {
        // E = KE + U_pair − Σ_i Ψ_i Q_i  (constant-potential grand energy)
        double mv2 = 0.0;
        for (int k = 0; k < 2; ++k)
            mv2 += kMion
                   * (vel[k][0] * vel[k][0] + vel[k][1] * vel[k][1]
                      + vel[k][2] * vel[k][2]);
        return 0.5 * units::mass_v2_to_energy * mv2 + f.meta.get(kEnergy)
               - (kPsiBot * q_plate(1) + kPsiTop * q_plate(2));
    };
    auto observe = [&]() {
        if (symplectic)
            res.max_resid = std::max(res.max_resid, maze.constraint_residual());
        res.max_neut =
            std::max(res.max_neut, maze.solver().neutrality_residual());
        double dx = f.span(kAtom_x)[kIon0] - f.span(kAtom_x)[kIon0 + 1];
        double dy = f.span(kAtom_y)[kIon0] - f.span(kAtom_y)[kIon0 + 1];
        double dz = f.span(kAtom_z)[kIon0] - f.span(kAtom_z)[kIon0 + 1];
        f.box.min_image(dx, dy, dz);
        res.min_sep = std::min(res.min_sep,
                               std::sqrt(dx * dx + dy * dy + dz * dz));
    };
    auto half_kick = [&]() {
        for (int k = 0; k < 2; ++k) {
            const int i = kIon0 + k;
            vel[k][0] += 0.5 * kDt * f.span(kAtom_fx)[i]
                         * units::force_to_accel / kMion;
            vel[k][1] += 0.5 * kDt * f.span(kAtom_fy)[i]
                         * units::force_to_accel / kMion;
            vel[k][2] += 0.5 * kDt * f.span(kAtom_fz)[i]
                         * units::force_to_accel / kMion;
        }
    };

    eval(); // initial charges + forces
    observe();
    res.q_top_first = q_plate(2);
    std::vector<double> q_start(maze.solver().charges());

    const int settle = 5; // skip bootstrap / ASPC-fill transient
    double e_ref = 0.0;
    int sampled = 0;
    auto do_leg = [&]() {
        for (int s = 0; s < kSteps; ++s) {
            half_kick();
            for (int k = 0; k < 2; ++k) {
                const int i = kIon0 + k;
                f.span(kAtom_x)[i] += kDt * vel[k][0];
                f.span(kAtom_y)[i] += kDt * vel[k][1];
                f.span(kAtom_z)[i] += kDt * vel[k][2];
            }
            eval();
            half_kick();
            observe();
            if (sampled == settle)
                e_ref = conserved();
            if (sampled > settle)
                res.drift =
                    std::max(res.drift, std::abs(conserved() - e_ref));
            ++sampled;
        }
    };

    do_leg(); // forward
    for (int k = 0; k < 2; ++k)
        for (int d = 0; d < 3; ++d)
            vel[k][d] = -vel[k][d]; // velocity flip
    do_leg();                       // backward (palindrome)

    for (int k = 0; k < 2; ++k) {
        const int i = kIon0 + k;
        res.ret_pos =
            std::max(res.ret_pos, std::abs(f.span(kAtom_x)[i] - x0[k]));
        res.ret_pos =
            std::max(res.ret_pos, std::abs(f.span(kAtom_y)[i] - y0[k]));
        res.ret_pos =
            std::max(res.ret_pos, std::abs(f.span(kAtom_z)[i] - z0[k]));
    }
    const auto& q_end = maze.solver().charges();
    for (std::size_t i = 0; i < q_end.size() && i < q_start.size(); ++i)
        res.ret_q = std::max(res.ret_q, std::abs(q_end[i] - q_start[i]));
    return res;
}

} // namespace

int main() {
    int fail = 0;

    const RunResult maze = run_case(/*symplectic=*/true, 0.0);
    const RunResult aspc = run_case(/*symplectic=*/false, /*cg_tol=*/1e-3);

    std::printf("maze-symplectic: drift=%.3e eV  ret_pos=%.3e Å  ret_q=%.3e e"
                "  max_resid=%.3e V  neut=%.3e e  min_sep=%.3f Å  Qtop=%.6f\n",
                maze.drift, maze.ret_pos, maze.ret_q, maze.max_resid,
                maze.max_neut, maze.min_sep, maze.q_top_first);
    std::printf("aspc-loose-twin: drift=%.3e eV  ret_pos=%.3e Å  ret_q=%.3e e"
                "  neut=%.3e e  min_sep=%.3f Å  Qtop=%.6f\n",
                aspc.drift, aspc.ret_pos, aspc.ret_q, aspc.max_neut,
                aspc.min_sep, aspc.q_top_first);

    // ── Validity of the scenario itself ─────────────────────────────────
    // Ions must stay well separated (bound smooth motion, no collapse).
    if (!(maze.min_sep > 1.8)) {
        std::printf("FAIL: ion collapse (min_sep=%.3f Å <= 1.8)\n",
                    maze.min_sep);
        ++fail;
    }
    // +1 V plate polarizes positive vs the 0 V plate (physical polarity).
    if (!(maze.q_top_first > 0.0)) {
        std::printf("FAIL: Q_top=%.3e <= 0 at ΔΨ=+1 V\n", maze.q_top_first);
        ++fail;
    }

    // ── MaZe constraint pins (Coretti 2020: constraint at solver precision,
    //    electroneutrality exact) ─────────────────────────────────────────
    if (!(maze.max_resid <= 1e-10)) {
        std::printf("FAIL: constraint residual %.3e V > 1e-10\n",
                    maze.max_resid);
        ++fail;
    }
    if (!(maze.max_neut <= 1e-9)) {
        std::printf("FAIL: electroneutrality %.3e e > 1e-9\n", maze.max_neut);
        ++fail;
    }

    // ── Time reversibility (Coretti 2020 §III: mass-zero constrained
    //    dynamics is time-reversible; dissipative predictor twins are not) ─
    if (!(maze.ret_pos <= 1e-6)) {
        std::printf("FAIL: MaZe palindrome position error %.3e Å > 1e-6\n",
                    maze.ret_pos);
        ++fail;
    }
    if (!(maze.ret_q <= 1e-6)) {
        std::printf("FAIL: MaZe palindrome charge error %.3e e > 1e-6\n",
                    maze.ret_q);
        ++fail;
    }
    if (!(maze.ret_pos < 0.1 * aspc.ret_pos)) {
        std::printf("FAIL: MaZe not >10x more reversible than loose ASPC twin "
                    "(maze=%.3e aspc=%.3e Å)\n",
                    maze.ret_pos, aspc.ret_pos);
        ++fail;
    }

    // ── Energy conservation (Coretti 2020 fig. 2-style: no secular drift
    //    for MaZe; loose-tolerance adiabatic twin drifts) ────────────────
    if (!(maze.drift <= 1e-6)) {
        std::printf("FAIL: MaZe energy drift %.3e eV > 1e-6\n", maze.drift);
        ++fail;
    }
    if (!(maze.drift < aspc.drift)) {
        std::printf("FAIL: MaZe drift %.3e not below loose ASPC twin %.3e\n",
                    maze.drift, aspc.drift);
        ++fail;
    }

    if (fail == 0)
        std::printf("PASS: maze-symplectic-charge-dynamics — constraint at "
                    "1e-10 V, reversible palindrome, drift below the loose "
                    "CG+ASPC twin\n");
    return fail ? 1 : 0;
}
