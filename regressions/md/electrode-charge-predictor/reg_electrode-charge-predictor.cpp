/// @file regressions/md/electrode-charge-predictor/reg_electrode-charge-predictor.cpp
/// @brief E2E: ASPC electrode charge predictor reduces CG iters vs cold start.
///
/// Public API only. Offline pin: regressions/md/electrode-charge-predictor/REFERENCE.md
/// (oracle=metalwalls, Kolafa ASPC + MetalWalls CG warm-start algorithm).
/// Exits non-zero on failure.

#include <cmath>
#include <cstdio>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/dynamics/const_pot_bo.h"
#include "atomiverse/cpu/electrodes/charge_aspc.h"
#include "atomiverse/cpu/frame/frame.h"

using namespace ATV;

namespace {

/// Small madelung for CG-vs-dense pin (2 electrode sites).
cpu::Frame make_madelung() {
    cpu::Frame f(4);
    const double z[4] = {8.0, 12.0, 9.0, 11.0};
    const double q0[4] = {0.0, 0.0, 0.5, -0.5};
    const int typ[4] = {1, 2, 3, 3};
    for (int i = 0; i < 4; ++i) {
        f.span(kAtom_x)[i] = 0.0;
        f.span(kAtom_y)[i] = 0.0;
        f.span(kAtom_z)[i] = z[i];
        f.span(kAtom_charge)[i] = q0[i];
        f.span(kAtom_type)[i] = typ[i];
        f.span(kAtom_mass)[i] = (typ[i] < 3) ? 196.97 : 1.0;
        f.span(kAtom_tag)[i] = i + 1;
    }
    f.simbox = cpu::Box::orthorhombic(1.0, 1.0, 20.0);
    return f;
}

/// Multi-site parallel plate — enough electrodes for CG iters to matter.
cpu::Frame make_multi(int n_side = 3) {
    const int n_per = n_side * n_side;
    const int n_ions = 4;
    const int N = 2 * n_per + n_ions;
    cpu::Frame f(N);
    const double Lxy = 6.0;
    int idx = 0;
    for (int wall = 0; wall < 2; ++wall) {
        const double z = (wall == 0) ? 4.0 : 16.0;
        const int typ = (wall == 0) ? 1 : 2;
        for (int ix = 0; ix < n_side; ++ix)
            for (int iy = 0; iy < n_side; ++iy) {
                f.span(kAtom_x)[idx] = (ix + 0.5) * (Lxy / n_side);
                f.span(kAtom_y)[idx] = (iy + 0.5) * (Lxy / n_side);
                f.span(kAtom_z)[idx] = z;
                f.span(kAtom_type)[idx] = typ;
                f.span(kAtom_charge)[idx] = 0.0;
                f.span(kAtom_mass)[idx] = 12.0;
                f.span(kAtom_tag)[idx] = idx + 1;
                ++idx;
            }
    }
    for (int i = 0; i < n_ions; ++i) {
        const bool pos = (i % 2 == 0);
        f.span(kAtom_x)[idx] = (0.5 + (i % 2) * 0.4) * Lxy;
        f.span(kAtom_y)[idx] = (0.3 + (i / 2) * 0.4) * Lxy;
        f.span(kAtom_z)[idx] = pos ? 8.0 : 12.0;
        f.span(kAtom_type)[idx] = 3;
        f.span(kAtom_charge)[idx] = pos ? 0.5 : -0.5;
        f.span(kAtom_mass)[idx] = 1.0;
        f.span(kAtom_tag)[idx] = idx + 1;
        ++idx;
    }
    f.simbox = cpu::Box::orthorhombic(Lxy, Lxy, 20.0);
    return f;
}

cpu::md::ConstPotBO make_solver(int pred_order) {
    cpu::md::ConstPotBO c(2.0);
    c.add_electrode(1, 0.0);
    c.add_electrode(2, 1.0);
    c.electroneutrality = true;
    c.xy_image_range = 0;
    c.use_cg = true;
    c.cg_tol = 1e-10;
    c.max_cg_iters = 500;
    c.predictor_order = pred_order;
    return c;
}

} // namespace

int main() {
    // ── ac-002 / ac-003: ASPC K=2 coefficients + linear predict ──────
    {
        cpu::electrodes::ChargeASPC aspc(2);
        if (!(aspc.coefficient(0) == 1.0 && aspc.coefficient(1) == -3.0
              && aspc.coefficient(2) == 3.0)) {
            std::printf(
                "FAIL: ASPC K=2 coeffs got {%g,%g,%g} expect {1,-3,3}\n",
                aspc.coefficient(0), aspc.coefficient(1), aspc.coefficient(2));
            return 1;
        }
        for (int n = 0; n <= 2; ++n)
            aspc.record(std::vector<double>{static_cast<double>(n)});
        std::vector<double> pred;
        if (!aspc.predict(pred) || pred.size() != 1
            || std::fabs(pred[0] - 3.0) > 1e-12) {
            std::printf("FAIL: linear ASPC predict (got %g)\n",
                        pred.empty() ? -1.0 : pred[0]);
            return 1;
        }
        std::printf("ASPC K=2: omega={1,-3,3}, linear predict OK\n");
    }

    // ── ac-004: CG matches dense on madelung ─────────────────────────
    {
        auto f_d = make_madelung();
        auto f_c = make_madelung();
        cpu::md::ConstPotBO dense(2.0);
        dense.add_electrode(1, 0.0).add_electrode(2, 1.0);
        dense.use_cg = false;
        dense.xy_image_range = 0;
        cpu::md::ConstPotBO cg = make_solver(/*pred_order=*/0);
        cpu::Context ctx;
        dense.build({ctx, f_d});
        dense.launch({ctx, f_d, 0, 0});
        cg.build({ctx, f_c});
        cg.launch({ctx, f_c, 0, 0});
        const double dQ = std::fabs(cg.total_charge(1) - dense.total_charge(1));
        if (dQ > 1e-8) {
            std::printf("FAIL: CG vs dense |ΔQ_top|=%.3e > 1e-8\n", dQ);
            return 1;
        }
        std::printf("CG vs dense: Q_top dense=%.8e cg=%.8e Δ=%.2e\n",
                    dense.total_charge(1), cg.total_charge(1), dQ);
    }

    // ── ac-005 / ac-006: predictor reduces mean CG iters ─────────────
    {
        auto cold = make_solver(0);
        auto warm = make_solver(2);
        auto f_cold = make_multi(3);
        auto f_warm = make_multi(3);
        cpu::Context ctx;
        cold.build({ctx, f_cold});
        warm.build({ctx, f_warm});

        const int n = f_cold.nlocal();
        const int ion0 = n - 4;
        std::printf("Multi-site: n_atoms=%d n_electrode=%d\n", n,
                    cold.n_electrode());

        long sum_cold = 0, sum_warm = 0;
        int n_post = 0;
        constexpr int n_steps = 16;
        for (int s = 0; s < n_steps; ++s) {
            const double dz = 0.03 * (s + 1);
            for (auto* f : {&f_cold, &f_warm}) {
                for (int k = 0; k < 4; ++k) {
                    const double z0 = (k % 2 == 0) ? 8.0 : 12.0;
                    const double sign = (k % 2 == 0) ? 1.0 : -0.5;
                    f->span(kAtom_z)[ion0 + k] = z0 + sign * dz;
                }
            }
            cold.launch({ctx, f_cold, 0, s});
            warm.launch({ctx, f_warm, 0, s});
            if (s >= 3) {
                sum_cold += cold.last_iters();
                sum_warm += warm.last_iters();
                ++n_post;
                if (!warm.last_used_predictor()) {
                    std::printf("FAIL: expected ASPC warm-start at step %d\n",
                                s);
                    return 1;
                }
            }
            if (!(warm.neutrality_residual() <= 1e-8)) {
                std::printf("FAIL: neutrality %.3e at step %d\n",
                            warm.neutrality_residual(), s);
                return 1;
            }
        }
        if (n_post <= 0) {
            std::printf("FAIL: no post-fill steps\n");
            return 1;
        }
        const double mean_cold = static_cast<double>(sum_cold) / n_post;
        const double mean_warm = static_cast<double>(sum_warm) / n_post;
        std::printf("Mean CG iters (post-fill n=%d): cold=%.3f warm=%.3f\n",
                    n_post, mean_cold, mean_warm);
        if (!(mean_warm < mean_cold)) {
            std::printf("FAIL: ASPC did not reduce mean CG iters "
                        "(warm=%.3f cold=%.3f)\n",
                        mean_warm, mean_cold);
            return 1;
        }
        if (!(warm.total_charge(1) > 0.0 && warm.total_charge(0) < 0.0)) {
            std::printf("FAIL: charge polarity lost Q_bot=%.6e Q_top=%.6e\n",
                        warm.total_charge(0), warm.total_charge(1));
            return 1;
        }
    }

    std::printf("PASS: electrode-charge-predictor — ASPC coeffs, CG=dense, "
                "fewer CG iters with ASPC\n");
    return 0;
}
