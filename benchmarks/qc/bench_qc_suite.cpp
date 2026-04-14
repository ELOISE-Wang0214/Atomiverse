/// @file benchmarks/qc/bench_qc_suite.cpp
/// @brief Phase 1.2 QC benchmark suite — LBFGS + SCFEvaluator at RHF
///        / 6-31G, against PySCF geomeTRIC-optimized references.
///
/// RHF / 6-31G is the largest basis the 2e-gradient path currently
/// produces physically correct forces at. def2-SVP energy works (SCF
/// parity to 1e-10 Ha) but the d-shell gradient path still disagrees
/// with FD; tracked in phase1-productization.md §1.2.
///
/// Build / run:
///     cmake --build build --target bench_qc_suite
///     ./build/benchmarks/bench_qc_suite
///
/// Regenerate references with:
///     python benchmarks/qc/reference/generate_references.py \
///       --method rhf --basis 6-31g
#ifdef ATV_HAS_RHF

#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;

// ── Molecule table ──────────────────────────────────────────────────
//
// Each entry holds the starting geometry (Angstrom, atoms listed Z then
// xyz) plus the reference-optimized energy in Hartree from PySCF
// PBE/def2-SVP (cart=True). Keep in sync with the xyz and reference
// JSON under benchmarks/qc/.

struct Atom { int Z; double x, y, z; };

struct MolSpec {
    const char* name;
    std::vector<Atom> atoms;
    double e_ref_hartree_optimized;  // from reference/<name>.json
};

// Reference energies: PySCF RHF / 6-31G / cart=True on each molecule's
// geomeTRIC-optimized geometry. Regenerate via
// `python benchmarks/qc/reference/generate_references.py --method rhf
// --basis 6-31g` and sync both files if geometries change.
static std::vector<MolSpec> build_molecule_suite() {
    return {
        {
            "h2o",
            {
                {8,  0.000000,  0.000000,  0.117790},
                {1,  0.000000,  0.755453, -0.471161},
                {1,  0.000000, -0.755453, -0.471161},
            },
            -75.9853591763,
        },
        {
            "nh3",
            {
                {7,  0.000000,  0.000000,  0.113880},
                {1,  0.000000,  0.938035, -0.265825},
                {1,  0.812328, -0.469017, -0.265825},
                {1, -0.812328, -0.469017, -0.265825},
            },
            -56.1655213295,
        },
        {
            "ch4",
            {
                {6,  0.000000,  0.000000,  0.000000},
                {1,  0.627640,  0.627640,  0.627640},
                {1, -0.627640, -0.627640,  0.627640},
                {1,  0.627640, -0.627640, -0.627640},
                {1, -0.627640,  0.627640, -0.627640},
            },
            -40.1805541783,
        },
    };
}

static cpu::Frame to_frame(const MolSpec& m) {
    cpu::Frame f(static_cast<int>(m.atoms.size()));
    for (int i = 0; i < static_cast<int>(m.atoms.size()); ++i) {
        f.atoms.type[i] = m.atoms[i].Z;
        f.atoms.x[i]    = m.atoms[i].x;
        f.atoms.y[i]    = m.atoms[i].y;
        f.atoms.z[i]    = m.atoms[i].z;
    }
    f.box = Mat3::cubic(100.0);
    return f;
}

// ── Per-molecule run ────────────────────────────────────────────────

struct RunResult {
    bool   converged;
    int    opt_steps;
    double wall_ms;
    double final_energy_eV;
    double final_energy_hartree;
    double delta_hartree_vs_ref;
};

static RunResult run_molecule(const MolSpec& m) {
    auto frame = to_frame(m);

    auto integrals = cpu::scf::build_integrals(
        frame, cpu::scf::BasisSet::Pople631G);

    auto cfg = cpu::scf::rhf_config(integrals, /*tol=*/1e-8);
    cfg.max_iter = 150;

    // Tight tolerances — we compare against the PySCF-geometric
    // optimum, which uses ~1e-4 Ha/Bohr on forces. ~1e-3 eV/Å matches.
    cpu::minimizer::LBFGS opt;
    opt.add(cpu::minimizer::LBFGSStep{
        /*memory=*/5,
        /*alpha0=*/0.01,
        /*ftol  =*/1e-3,
        /*max_ls=*/10,
        /*rms_ftol=*/5e-4,
        /*max_dtol=*/1e-2,
        /*rms_dtol=*/5e-3,
        /*etol=*/1e-6,
    });
    opt.on(ATV::minimizer::Phase::force).add(
        cpu::scf::SCFEvaluator{cfg, cpu::scf::BasisSet::Pople631G});
    opt.build();

    const auto t0 = std::chrono::steady_clock::now();
    const int max_steps = 100;
    opt.run(frame, max_steps);
    const double wall_ms =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t0).count();

    // Re-evaluate final energy in Hartree on the converged geometry for
    // direct comparison with the reference (which is in Hartree).
    auto ints_final = cpu::scf::build_integrals(
        frame, cpu::scf::BasisSet::Pople631G);
    auto cfg_final = cpu::scf::rhf_config(ints_final, /*tol=*/1e-10);
    auto res = cpu::scf::scf_solve(ints_final, cfg_final);

    RunResult out;
    out.converged             = res.converged;
    out.opt_steps             = max_steps;  // LBFGSStep doesn't expose a public step-count hook
    out.wall_ms               = wall_ms;
    out.final_energy_eV       = frame.energy;
    out.final_energy_hartree  = res.energy;
    out.delta_hartree_vs_ref  = res.energy - m.e_ref_hartree_optimized;
    return out;
}

int main() {
    std::printf("# Phase 1.2 QC benchmark suite — LBFGS + SCFEvaluator at RHF / 6-31G\n");
    std::printf("#\n");
    std::printf("# NOTE: LBFGSStep has no line search yet, so fixed-step\n");
    std::printf("# convergence can leave a sub-mHa residue vs the true PySCF\n");
    std::printf("# minimum. Threshold is 1 mHa; tighter agreement needs a\n");
    std::printf("# proper line search (tracked in phase1-productization.md).\n");
    std::printf("#\n");
    std::printf("# %-6s  %-10s  %-10s  %-16s  %-14s  %-10s\n",
                "mol", "conv", "wall_ms", "E_final [Ha]", "dE_vs_ref [Ha]", "status");

    const double threshold_ha = 1e-3;
    int n_fail = 0;

    for (const auto& mol : build_molecule_suite()) {
        const auto r = run_molecule(mol);
        const bool within = std::abs(r.delta_hartree_vs_ref) <= threshold_ha;
        const char* status =
            !r.converged ? "SCF FAIL"
                         : (within ? "PASS" : "DRIFT");

        std::printf("  %-6s  %-10s  %10.1f  %16.10f  %+14.2e  %-10s\n",
                    mol.name,
                    r.converged ? "yes" : "no",
                    r.wall_ms,
                    r.final_energy_hartree,
                    r.delta_hartree_vs_ref,
                    status);

        if (!r.converged || !within) ++n_fail;
    }

    std::printf("#\n# %d / %d molecules within %.0e Ha of PySCF reference.\n",
                static_cast<int>(build_molecule_suite().size()) - n_fail,
                static_cast<int>(build_molecule_suite().size()),
                threshold_ha);

    return (n_fail == 0) ? 0 : 1;
}

#else
#include <cstdio>
int main() {
    std::printf("bench_qc_suite requires -DATV_BUILD_RHF=ON\n");
    return 1;
}
#endif
