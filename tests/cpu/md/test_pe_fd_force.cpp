/// @file tests/cpu/md/test_pe_fd_force.cpp
/// @brief Finite-difference consistency check for each bonded force kernel.
///
/// For a perturbed C4 backbone (so every kernel sees a non-trivial geometry),
/// compare the analytic force returned by HarmonicBond, HarmonicAngle and
/// MultiHarmonicDihedral against a central-difference gradient of the same
/// kernel's energy. Mismatch pinpoints which kernel's force formula is bad.

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/context.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/potentials/harmonic_angle.h"
#include "atomiverse/cpu/potentials/harmonic_bond.h"
#include "atomiverse/cpu/potentials/multi_harmonic_dihedral.h"

using namespace ATV;

namespace {

cpu::Frame make_perturbed_butane() {
    cpu::Frame frame(4);
    // Off-trans geometry so cos(phi) != -1, sin(phi) != 0.
    frame.atoms.x[0] = 0.0;
    frame.atoms.y[0] = 0.0;
    frame.atoms.z[0] = 0.0;
    frame.atoms.x[1] = 1.5;
    frame.atoms.y[1] = 0.0;
    frame.atoms.z[1] = 0.0;
    frame.atoms.x[2] = 2.1;
    frame.atoms.y[2] = 1.3;
    frame.atoms.z[2] = 0.0;
    frame.atoms.x[3] = 3.5;
    frame.atoms.y[3] = 1.1;
    frame.atoms.z[3] = 0.4;
    for (int i = 0; i < 4; ++i) {
        frame.atoms.mass[i] = 14.0;
        frame.atoms.type[i] = 0;
    }
    frame.box = Mat3::cubic(50.0);

    frame.atoms.map_clear();
    frame.atoms.map_set(0, frame.atoms.nlocal);
    for (int k = 0; k < 3; ++k) {
        frame.atoms.add_bond(k, k + 1, 0);
    }
    for (int k = 0; k < 2; ++k) {
        frame.atoms.add_angle(k, k + 1, k + 2, 0);
    }
    frame.atoms.add_dihedral(0, 1, 2, 3, 0);
    frame.atoms.build_special_list();
    return frame;
}

// Return the analytic force on atom `atom_idx`, axis `d` (0=x,1=y,2=z),
// after running `launch` on a freshly-zeroed-force copy of `frame`.
double analytic_force(cpu::Frame frame, std::function<void(cpu::Frame&)> launch,
                      int atom_idx, int axis) {
    for (int i = 0; i < frame.atoms.size(); ++i) {
        frame.atoms.fx[i] = 0;
        frame.atoms.fy[i] = 0;
        frame.atoms.fz[i] = 0;
    }
    frame.scalars.energy = 0;
    launch(frame);
    if (axis == 0)
        return frame.atoms.fx[atom_idx];
    if (axis == 1)
        return frame.atoms.fy[atom_idx];
    return frame.atoms.fz[atom_idx];
}

// Central-difference estimate of dE/dx[atom, axis].
double fd_gradient(cpu::Frame frame, std::function<void(cpu::Frame&)> launch,
                   int atom_idx, int axis, double h) {
    auto eval = [&](cpu::Frame& f) {
        for (int i = 0; i < f.atoms.size(); ++i) {
            f.atoms.fx[i] = 0;
            f.atoms.fy[i] = 0;
            f.atoms.fz[i] = 0;
        }
        f.scalars.energy = 0;
        launch(f);
        return f.scalars.energy;
    };
    auto& pos = (axis == 0)   ? frame.atoms.x
                : (axis == 1) ? frame.atoms.y
                              : frame.atoms.z;
    const double orig = pos[atom_idx];
    pos[atom_idx] = orig + h;
    const double e_plus = eval(frame);
    pos[atom_idx] = orig - h;
    const double e_minus = eval(frame);
    pos[atom_idx] = orig;
    return (e_plus - e_minus) / (2.0 * h);
}

void check_kernel(const std::string& label,
                  std::function<void(cpu::Frame&)> launch, double tol) {
    auto frame = make_perturbed_butane();
    for (int a = 0; a < frame.atoms.size(); ++a) {
        for (int d = 0; d < 3; ++d) {
            const double f_analytic = analytic_force(frame, launch, a, d);
            const double grad = fd_gradient(frame, launch, a, d, 1e-4);
            // Analytic force = -grad; compare.
            EXPECT_NEAR(f_analytic, -grad, tol)
                << label << " atom " << a << " axis " << d
                << "  analytic=" << f_analytic << "  -grad=" << -grad;
        }
    }
}

} // namespace

TEST(PeFdForce, HarmonicBond) {
    check_kernel(
        "HarmonicBond",
        [](cpu::Frame& f) {
            cpu::Context ctx;
            cpu::BuildContext bctx{ctx, f};
            cpu::LaunchContext lctx{ctx, f, 0, 0};
            HarmonicBond c({1.0}, {1.5});
            c.build(bctx);
            c.launch(lctx);
            c.destroy({ctx, f});
        },
        1e-6);
}

TEST(PeFdForce, HarmonicAngle) {
    check_kernel(
        "HarmonicAngle",
        [](cpu::Frame& f) {
            cpu::Context ctx;
            cpu::BuildContext bctx{ctx, f};
            cpu::LaunchContext lctx{ctx, f, 0, 0};
            HarmonicAngle c({1.0}, {2.0}); // K=1, theta0=2 rad
            c.build(bctx);
            c.launch(lctx);
            c.destroy({ctx, f});
        },
        1e-6);
}

TEST(PeFdForce, MultiHarmonicDihedral) {
    check_kernel(
        "MultiHarmonicDihedral",
        [](cpu::Frame& f) {
            cpu::Context ctx;
            cpu::BuildContext bctx{ctx, f};
            cpu::LaunchContext lctx{ctx, f, 0, 0};
            MultiHarmonicDihedral c({{0.1, -0.2, 0.05, 0.3, -0.1}});
            c.build(bctx);
            c.launch(lctx);
            c.destroy({ctx, f});
        },
        1e-6);
}
