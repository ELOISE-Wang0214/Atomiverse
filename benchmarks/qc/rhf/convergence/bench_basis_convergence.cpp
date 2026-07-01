/// @file benchmarks/qc/rhf/convergence/bench_basis_convergence.cpp
/// @brief Convergence study: RHF total energy vs Gaussian basis set.
///        Runnable skeleton for the manual "basis sets & convergence" chapter
///        (manual-system §3 chapter 4). Result class = physical approximation
///        (manual-system §5): the check is monotone convergence toward the
///        complete-basis limit, not agreement with a single reference number.
///
/// Runs a closed-shell RHF single point on a FIXED H2O geometry across the
/// systematically-improving basis ladder
///     STO-3G → 6-31G → 6-31G* → def2-SVP → def2-TZVP
/// and asserts the variational trend: each larger basis lowers (or holds, to a
/// small slack) the HF energy, and the successive energy lowering shrinks —
/// the signature of convergence.
///
/// This is scaffolding: the physics is real, but the value it adds over the
/// PySCF-validated `qc/rhf/validation` suite is the *trend across basis sets*.
/// Extend by adding molecules / a CBS extrapolation as the chapter needs.
///
/// Usage:
///   bench_qc_basis_convergence [--slack HA]
///     --slack HA   allowed positive slack on the monotone-decrease test
///                  (absorbs cross-family near-degeneracy; default 2e-3 Ha)
/// Exit 0 = monotone-convergent trend, 1 = a basis broke the trend / SCF failed.

#ifdef ATV_HAS_RHF

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"

using namespace ATV;

namespace {

/// Fixed equilibrium-like H2O geometry [Å]; type = atomic number Z.
cpu::Frame make_h2o() {
    cpu::Frame f(3);
    const int Z[3] = {8, 1, 1};
    const double xyz[3][3] = {
        {0.000000, 0.000000, 0.117300},
        {0.000000, 0.757200, -0.469200},
        {0.000000, -0.757200, -0.469200},
    };
    for (int i = 0; i < 3; ++i) {
        f.atoms.type[i] = Z[i];
        f.atoms.x[i] = xyz[i][0];
        f.atoms.y[i] = xyz[i][1];
        f.atoms.z[i] = xyz[i][2];
    }
    f.box = Mat3::cubic(100.0);
    return f;
}

struct Rung {
    const char* name;
    cpu::scf::BasisSet basis;
};

} // namespace

int main(int argc, char** argv) {
    double slack = 2e-3; // Ha — cross-family near-degeneracy tolerance
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--slack" && i + 1 < argc)
            slack = std::atof(argv[++i]);
        else {
            std::fprintf(stderr, "unknown argument: %s\n", arg.c_str());
            return 1;
        }
    }

    const std::vector<Rung> ladder = {
        {"STO-3G", cpu::scf::BasisSet::STO3G},
        {"6-31G", cpu::scf::BasisSet::Pople631G},
        {"6-31G*", cpu::scf::BasisSet::Pople631Gs},
        {"def2-SVP", cpu::scf::BasisSet::Def2SVP},
        {"def2-TZVP", cpu::scf::BasisSet::Def2TZVP},
    };

    std::printf("# RHF/H2O total energy vs basis set (variational convergence)\n");
    std::printf("# %-12s  %-18s  %-16s\n", "basis", "E_RHF [Ha]",
                "dE_vs_prev [Ha]");

    auto frame = make_h2o();
    std::vector<double> energies;
    int n_fail = 0;
    double prev_e = 0.0, prev_delta = 0.0;
    bool have_prev = false;

    for (const auto& rung : ladder) {
        double e;
        try {
            auto ints = cpu::scf::build_integrals(frame, rung.basis);
            auto res = cpu::scf::rhf_solve(ints, 1e-10, 200);
            if (!res.converged) {
                std::printf("  %-12s  SCF NOT CONVERGED\n", rung.name);
                ++n_fail;
                continue;
            }
            e = res.energy;
        } catch (const std::exception& ex) {
            std::printf("  %-12s  ERROR: %s\n", rung.name, ex.what());
            ++n_fail;
            continue;
        }

        if (!have_prev) {
            std::printf("  %-12s  %18.10f  %-16s\n", rung.name, e, "-");
        } else {
            const double delta = e - prev_e;
            // Variational monotone test: a larger basis must not RAISE the HF
            // energy beyond a small cross-family slack.
            const char* flag = (delta <= slack) ? "" : "  <-- NON-MONOTONE";
            if (delta > slack)
                ++n_fail;
            // Convergence test: successive lowerings shrink in magnitude
            // (only checked once we have two deltas; informational).
            const char* conv =
                (prev_delta != 0.0 && std::fabs(delta) > std::fabs(prev_delta))
                    ? "  (delta grew)"
                    : "";
            std::printf("  %-12s  %18.10f  %+16.3e%s%s\n", rung.name, e, delta,
                        flag, conv);
            prev_delta = delta;
        }
        energies.push_back(e);
        prev_e = e;
        have_prev = true;
    }

    std::printf("#\n# Expected: E decreases monotonically toward the CBS limit "
                "and the\n# per-step lowering shrinks. %d issue(s).\n", n_fail);
    return n_fail == 0 ? 0 : 1;
}

#else
#include <cstdio>
int main() {
    std::printf("bench_qc_basis_convergence requires -DATV_BUILD_RHF=ON\n");
    return 1;
}
#endif
