/// @file benchmarks/pw/validation/bench_pw_si_cp2k_eos.cpp
/// @brief Scientific validation: Si fcc equation of state vs CP2K.
///
/// Runs the plane-wave SCF over a lattice-constant scan of bulk Si
/// (2-atom fcc primitive cell, ONCVPSP LDA pseudopotential), fits the
/// E(V) curve, and compares the equilibrium lattice constant a0 and
/// bulk modulus B0 against a CP2K Quickstep reference (same scan,
/// GTH-PADE-q4 / DZVP-GTH, LDA) fitted with the *same* code path.
///
/// Absolute energies are NOT comparable across codes (different
/// pseudopotentials and discretizations); EOS-derived observables are.
/// If the CP2K reference CSV has not been generated yet
/// (reference/cp2k/eos_reference.csv), falls back to literature LDA
/// values (a0 = 5.40 Å, B0 = 96 GPa) with a notice.
///
/// Usage: bench_pw_si_cp2k_eos [--ecut RY] [--kmesh N] [--backend cpu|cuda]
///                             [--out csv_path]
/// Exit code 0 = PASS, 1 = FAIL.

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/pw/pw_scf.h"
#ifdef ATV_BUILD_CUDA
#include "atomiverse/cuda/pw/pw_scf.cuh"
#endif

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace ATV;

namespace {

constexpr double kRyToEv = 13.605693122994;
constexpr double kHaToEv = 2.0 * kRyToEv;
constexpr double kEvPerA3ToGPa = 160.21766208;

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
    double h = a / 2.0;
    Mat3 cell;
    cell(0, 0) = 0;
    cell(0, 1) = h;
    cell(0, 2) = h;
    cell(1, 0) = h;
    cell(1, 1) = 0;
    cell(1, 2) = h;
    cell(2, 0) = h;
    cell(2, 1) = h;
    cell(2, 2) = 0;
    f.box = cell;
    return f;
}

struct EosFit {
    double v0 = 0;  ///< Equilibrium volume per 2-atom cell [Å³].
    double a0 = 0;  ///< Equilibrium fcc lattice constant [Å].
    double b0 = 0;  ///< Bulk modulus [GPa].
    bool ok = false;
};

/// Least-squares quadratic fit E(V) = c0 + c1·V + c2·V², then
/// V0 = -c1/(2c2), B0 = V0·2c2. Energies in eV, volumes in Å³.
EosFit fit_eos(const std::vector<double>& v, const std::vector<double>& e) {
    EosFit fit;
    size_t n = v.size();
    if (n < 3)
        return fit;
    // Normal equations for [1, V, V²] basis.
    double s[5] = {0, 0, 0, 0, 0}, t[3] = {0, 0, 0};
    for (size_t i = 0; i < n; ++i) {
        double p = 1;
        for (int k = 0; k < 5; ++k) {
            s[k] += p;
            if (k < 3)
                t[k] += e[i] * p;
            p *= v[i];
        }
    }
    double m[9] = {s[0], s[1], s[2], s[1], s[2], s[3], s[2], s[3], s[4]};
    double c[3] = {t[0], t[1], t[2]};
    // Gaussian elimination, 3×3.
    for (int k = 0; k < 3; ++k) {
        int piv = k;
        for (int i = k + 1; i < 3; ++i)
            if (std::abs(m[i * 3 + k]) > std::abs(m[piv * 3 + k]))
                piv = i;
        for (int j = 0; j < 3; ++j)
            std::swap(m[k * 3 + j], m[piv * 3 + j]);
        std::swap(c[k], c[piv]);
        if (std::abs(m[k * 3 + k]) < 1e-30)
            return fit;
        for (int i = k + 1; i < 3; ++i) {
            double f = m[i * 3 + k] / m[k * 3 + k];
            for (int j = k; j < 3; ++j)
                m[i * 3 + j] -= f * m[k * 3 + j];
            c[i] -= f * c[k];
        }
    }
    double coef[3];
    for (int k = 2; k >= 0; --k) {
        double sum = c[k];
        for (int j = k + 1; j < 3; ++j)
            sum -= m[k * 3 + j] * coef[j];
        coef[k] = sum / m[k * 3 + k];
    }
    if (coef[2] <= 0)
        return fit; // not a minimum
    fit.v0 = -coef[1] / (2.0 * coef[2]);
    fit.a0 = std::cbrt(4.0 * fit.v0); // primitive cell V = a³/4
    fit.b0 = fit.v0 * 2.0 * coef[2] * kEvPerA3ToGPa;
    fit.ok = true;
    return fit;
}

/// Read "a_angstrom,energy" CSV (# comments skipped). Energy unit is
/// declared by a "# energy_unit: hartree|rydberg|ev" header line.
bool read_eos_csv(const std::string& path, std::vector<double>& a,
                  std::vector<double>& e_ev) {
    std::ifstream f(path);
    if (!f)
        return false;
    double to_ev = kHaToEv; // default: hartree (CP2K native)
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty())
            continue;
        if (line[0] == '#') {
            if (line.find("energy_unit:") != std::string::npos) {
                if (line.find("rydberg") != std::string::npos)
                    to_ev = kRyToEv;
                else if (line.find("ev") != std::string::npos)
                    to_ev = 1.0;
            }
            continue;
        }
        std::istringstream ss(line);
        std::string col;
        double ai, ei;
        if (!std::getline(ss, col, ','))
            continue;
        ai = std::atof(col.c_str());
        if (!std::getline(ss, col, ','))
            continue;
        ei = std::atof(col.c_str());
        a.push_back(ai);
        e_ev.push_back(ei * to_ev);
    }
    return a.size() >= 3;
}

} // namespace

int main(int argc, char** argv) {
    double ecut_ry = 30.0;
    int kmesh = 2;
    std::string backend = "cpu";
    std::string out_csv = "atomiverse_si_eos.csv";
    std::string upf = std::string(ATV_PW_DATA_DIR) + "/Si_lda_nc.upf";
    std::string ref_csv =
        std::string(ATV_PW_DATA_DIR) + "/cp2k/eos_reference.csv";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto next = [&]() -> std::string {
            return (i + 1 < argc) ? argv[++i] : "";
        };
        if (arg == "--ecut")
            ecut_ry = std::atof(next().c_str());
        else if (arg == "--kmesh")
            kmesh = std::atoi(next().c_str());
        else if (arg == "--backend")
            backend = next();
        else if (arg == "--out")
            out_csv = next();
        else if (arg == "--upf")
            upf = next();
        else if (arg == "--ref")
            ref_csv = next();
    }

    const std::array<double, 7> scan = {5.25, 5.30, 5.35, 5.40,
                                        5.45, 5.50, 5.55};

    std::cout << "=== Si fcc EOS: Atomiverse pw (" << backend
              << ", ecut=" << ecut_ry << " Ry, k=" << kmesh << "x" << kmesh
              << "x" << kmesh << ") vs CP2K ===" << std::endl;

    std::vector<double> vols, energies_ev, lats;
    bool all_converged = true;
    for (double a : scan) {
        auto frame = make_si_fcc(a);
        cpu::pw::PWResult r;
        if (backend == "cuda") {
#ifdef ATV_BUILD_CUDA
            r = cuda::pw::pw_scf(frame, upf, ecut_ry, kmesh, 1e-7, 60);
#else
            std::cerr << "FAIL: built without ATV_BUILD_CUDA" << std::endl;
            return 1;
#endif
        } else {
            r = cpu::pw::pw_scf(frame, upf, ecut_ry, kmesh, 1e-7, 60);
        }
        if (!r.converged) {
            std::cerr << "  WARNING: SCF not converged at a=" << a
                      << std::endl;
            all_converged = false;
        }
        lats.push_back(a);
        vols.push_back(a * a * a / 4.0);
        energies_ev.push_back(r.energy * kRyToEv);
        std::printf("  a=%.2f A   E=%+.6f Ry (%d iters)\n", a, r.energy,
                    r.iterations);
    }

    {
        std::ofstream out(out_csv);
        out << "# Atomiverse pw Si fcc EOS (" << backend
            << "), 2-atom primitive cell\n";
        out << "# energy_unit: rydberg\n# a_angstrom,energy\n";
        for (size_t i = 0; i < lats.size(); ++i) {
            char buf[64];
            std::snprintf(buf, sizeof buf, "%.2f,%.8f\n", lats[i],
                          energies_ev[i] / kRyToEv);
            out << buf;
        }
    }

    auto atv = fit_eos(vols, energies_ev);
    if (!atv.ok || !all_converged) {
        std::cerr << "FAIL: Atomiverse EOS fit failed (no minimum in scan "
                     "range or unconverged SCF)"
                  << std::endl;
        return 1;
    }

    // ── Reference: CP2K scan if generated, literature LDA otherwise ──
    std::vector<double> ref_a, ref_e;
    EosFit ref;
    bool have_cp2k = read_eos_csv(ref_csv, ref_a, ref_e);
    double tol_a0 = 0.08, tol_b0 = 25.0;
    std::string ref_name;
    if (have_cp2k) {
        std::vector<double> ref_v;
        for (double ai : ref_a)
            ref_v.push_back(ai * ai * ai / 4.0);
        ref = fit_eos(ref_v, ref_e);
        ref_name = "CP2K (GTH-PADE-q4, DZVP-GTH, LDA)";
        if (!ref.ok) {
            std::cerr << "FAIL: CP2K reference fit failed (" << ref_csv
                      << ")" << std::endl;
            return 1;
        }
    } else {
        // Literature LDA: a0 ≈ 5.40 Å, B0 ≈ 96 GPa (see SOURCES.md).
        ref.a0 = 5.40;
        ref.b0 = 96.0;
        ref.ok = true;
        ref_name = "literature LDA (CP2K CSV not generated yet — run "
                   "reference/cp2k/run_eos.sh)";
        tol_a0 = 0.10;
        tol_b0 = 30.0;
    }

    std::printf("\n  %-28s a0 = %.4f A   B0 = %6.1f GPa\n", "Atomiverse pw:",
                atv.a0, atv.b0);
    std::printf("  %-28s a0 = %.4f A   B0 = %6.1f GPa\n",
                (ref_name.substr(0, 4) + ":").c_str(), ref.a0, ref.b0);
    std::printf("  reference: %s\n", ref_name.c_str());
    std::printf("  |d a0| = %.4f A (tol %.2f)   |d B0| = %.1f GPa (tol "
                "%.0f)\n",
                std::abs(atv.a0 - ref.a0), tol_a0, std::abs(atv.b0 - ref.b0),
                tol_b0);

    bool pass = std::abs(atv.a0 - ref.a0) <= tol_a0
                && std::abs(atv.b0 - ref.b0) <= tol_b0;
    std::cout << (pass ? "\nPASS" : "\nFAIL") << std::endl;
    return pass ? 0 : 1;
}
