/// @file benchmarks/pw/si/convergence/bench_ecut_kmesh_convergence.cpp
/// @brief Convergence study: plane-wave total energy vs cutoff (ecut) and vs
///        k-point mesh, for bulk Si. Runnable skeleton for the manual
///        "plane-wave DFT for solids" chapter (manual-system §3 chapter 8,
///        which calls for a convergence study alongside the EOS validation).
///
/// Result class = physical approximation (manual-system §5): the check is
/// monotone convergence toward the converged limit, not a single reference
/// energy (absolute PW energies are not comparable across codes/pseudos —
/// see the EOS validation bench for the pseudo-insensitive a0/B0 observables).
///
/// Two independent 1-D sweeps on a FIXED 2-atom fcc Si cell (a = 5.43 Å,
/// ONCVPSP LDA pseudopotential shipped under the pw/si validation reference):
///   1. ecut sweep at fixed k-mesh — VARIATIONAL: a larger plane-wave basis
///      lowers (or holds) the total energy; the lowering shrinks and plateaus.
///   2. k-mesh sweep at fixed ecut — Brillouin-zone sampling converges: the
///      energy is NOT variational in k but successive changes shrink.
///
/// Usage:
///   bench_pw_convergence [--backend cpu|cuda] [--ecut-list R,R,...]
///                        [--kmesh-fixed N] [--kmesh-list N,N,...]
///                        [--ecut-fixed R] [--slack RY] [--upf PATH]
/// Exit 0 = both sweeps converge monotonically, 1 = a point broke the trend.
///
/// NOTE: CPU plane-wave SCF is single-threaded and slow (minutes/point even at
/// low ecut). On the GH200 node prefer `--backend cuda`. The defaults are a
/// small, fast smoke sweep — production convergence uses larger ecut/k lists.

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/pw/pw_scf.h"
#ifdef ATV_BUILD_CUDA
#include "atomiverse/cuda/pw/pw_scf.cuh"
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

using namespace ATV;

namespace {

cpu::Frame make_si_fcc(double a) {
    cpu::Frame f(2);
    for (int i = 0; i < 2; ++i) {
        f.atoms.type[i] = 14;
        f.atoms.mass[i] = 28.0855;
    }
    f.atoms.x[0] = 0;
    f.atoms.y[0] = 0;
    f.atoms.z[0] = 0;
    f.atoms.x[1] = a / 4.0;
    f.atoms.y[1] = a / 4.0;
    f.atoms.z[1] = a / 4.0;
    const double h = a / 2.0;
    Mat3 cell;
    cell(0, 0) = 0;   cell(0, 1) = h;   cell(0, 2) = h;
    cell(1, 0) = h;   cell(1, 1) = 0;   cell(1, 2) = h;
    cell(2, 0) = h;   cell(2, 1) = h;   cell(2, 2) = 0;
    f.box = cell;
    return f;
}

std::vector<double> parse_doubles(const std::string& csv) {
    std::vector<double> out;
    std::stringstream ss(csv);
    std::string tok;
    while (std::getline(ss, tok, ','))
        if (!tok.empty())
            out.push_back(std::atof(tok.c_str()));
    return out;
}

std::vector<int> parse_ints(const std::string& csv) {
    std::vector<int> out;
    std::stringstream ss(csv);
    std::string tok;
    while (std::getline(ss, tok, ','))
        if (!tok.empty())
            out.push_back(std::atoi(tok.c_str()));
    return out;
}

double run_pw(const cpu::Frame& frame, const std::string& upf, double ecut_ry,
              int kmesh, const std::string& backend, bool& converged) {
    cpu::pw::PWResult r;
    if (backend == "cuda") {
#ifdef ATV_BUILD_CUDA
        r = cuda::pw::pw_scf(frame, upf, ecut_ry, kmesh, 1e-7, 60);
#else
        std::fprintf(stderr, "built without ATV_BUILD_CUDA\n");
        converged = false;
        return 0.0;
#endif
    } else {
        r = cpu::pw::pw_scf(frame, upf, ecut_ry, kmesh, 1e-7, 60);
    }
    converged = r.converged;
    return r.energy; // [Ry]
}

} // namespace

int main(int argc, char** argv) {
    std::string backend = "cpu";
    std::string upf = std::string(ATV_PW_DATA_DIR) + "/Si_lda_nc.upf";
    std::vector<double> ecut_list = {15.0, 20.0, 25.0, 30.0};
    int kmesh_fixed = 2;
    std::vector<int> kmesh_list = {2, 3, 4};
    double ecut_fixed = 20.0;
    double slack = 1e-4; // Ry — variational slack on the ecut monotone test
    const double a_si = 5.43;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&]() -> std::string {
            return (i + 1 < argc) ? argv[++i] : "";
        };
        if (arg == "--backend")
            backend = next();
        else if (arg == "--upf")
            upf = next();
        else if (arg == "--ecut-list")
            ecut_list = parse_doubles(next());
        else if (arg == "--kmesh-fixed")
            kmesh_fixed = std::atoi(next().c_str());
        else if (arg == "--kmesh-list")
            kmesh_list = parse_ints(next());
        else if (arg == "--ecut-fixed")
            ecut_fixed = std::atof(next().c_str());
        else if (arg == "--slack")
            slack = std::atof(next().c_str());
        else {
            std::fprintf(stderr, "unknown argument: %s\n", arg.c_str());
            return 1;
        }
    }

    auto frame = make_si_fcc(a_si);
    int n_fail = 0;

    // ── Sweep 1: ecut at fixed k-mesh (variational) ───────────────────
    std::printf("# PW Si total energy vs ecut  (k-mesh %dx%dx%d, %s)\n",
                kmesh_fixed, kmesh_fixed, kmesh_fixed, backend.c_str());
    std::printf("# %-10s  %-18s  %-16s\n", "ecut[Ry]", "E[Ry]", "dE_vs_prev[Ry]");
    double prev = 0.0, prev_delta = 0.0;
    bool have_prev = false;
    for (double ecut : ecut_list) {
        bool conv = false;
        const double e = run_pw(frame, upf, ecut, kmesh_fixed, backend, conv);
        if (!conv) {
            std::printf("  %-10.1f  SCF NOT CONVERGED\n", ecut);
            ++n_fail;
            continue;
        }
        if (!have_prev) {
            std::printf("  %-10.1f  %18.8f  %-16s\n", ecut, e, "-");
        } else {
            const double d = e - prev;
            const char* flag = (d <= slack) ? "" : "  <-- NON-MONOTONE";
            if (d > slack)
                ++n_fail;
            const char* conv_note =
                (prev_delta != 0.0 && std::fabs(d) > std::fabs(prev_delta))
                    ? "  (delta grew)"
                    : "";
            std::printf("  %-10.1f  %18.8f  %+16.3e%s%s\n", ecut, e, d, flag,
                        conv_note);
            prev_delta = d;
        }
        prev = e;
        have_prev = true;
    }

    // ── Sweep 2: k-mesh at fixed ecut (BZ sampling) ───────────────────
    std::printf("#\n# PW Si total energy vs k-mesh  (ecut %.1f Ry, %s)\n",
                ecut_fixed, backend.c_str());
    std::printf("# %-10s  %-18s  %-16s\n", "kmesh", "E[Ry]", "dE_vs_prev[Ry]");
    prev = 0.0;
    prev_delta = 0.0;
    have_prev = false;
    for (int k : kmesh_list) {
        bool conv = false;
        const double e = run_pw(frame, upf, ecut_fixed, k, backend, conv);
        if (!conv) {
            std::printf("  %-10d  SCF NOT CONVERGED\n", k);
            ++n_fail;
            continue;
        }
        if (!have_prev) {
            std::printf("  %dx%dx%d      %18.8f  %-16s\n", k, k, k, e, "-");
        } else {
            const double d = e - prev;
            // k-mesh energy is NOT variational — only require the successive
            // change to SHRINK (Brillouin-zone sampling converging).
            const char* flag =
                (prev_delta != 0.0 && std::fabs(d) > std::fabs(prev_delta))
                    ? "  <-- NOT CONVERGING (delta grew)"
                    : "";
            if (prev_delta != 0.0 && std::fabs(d) > std::fabs(prev_delta))
                ++n_fail;
            std::printf("  %dx%dx%d      %18.8f  %+16.3e%s\n", k, k, k, e, d,
                        flag);
            prev_delta = d;
        }
        prev = e;
        have_prev = true;
    }

    std::printf("#\n# Expected: E(ecut) decreases monotonically and plateaus; "
                "E(k-mesh)\n# converges (|dE| shrinks). %d issue(s).\n", n_fail);
    return n_fail == 0 ? 0 : 1;
}
