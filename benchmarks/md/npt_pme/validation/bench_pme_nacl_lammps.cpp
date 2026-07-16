/// @file benchmarks/md/npt_pme/validation/bench_pme_nacl_lammps.cpp
/// @brief Validation: Atomiverse CoulLong + KSpacePME vs LAMMPS coul/long +
///        pppm on three NaCl rock-salt cells (s = 1.00, 0.99, 1.02).
///
/// Reference: validation/reference/frames_lammps.txt + SOURCES.md
/// (LAMMPS 22 Jul 2025 Update 4, metal units, pinned α/mesh/order).
///
/// Per frame checks energy [eV], six virial components [eV] (from LAMMPS
/// bar pressure via nktv2p), volume, box edge, and force sum/max (cubic
/// symmetry → forces vanish). Not a ctest gate — run under
/// -DATV_BUILD_BENCHMARKS=ON at acceptance time.

#include "bench_check.h"

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/driver/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/kspace/pme.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/coul_long.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace ATV;

namespace {

constexpr double kD = 2.82;
constexpr int kSide = 4;
constexpr double kCutoff = 1.9 * kD;
constexpr double kAlpha = 3.5 / kCutoff;
constexpr double kSpacing = 0.5;
constexpr int kOrder = 6;
// LAMMPS metal: P[bar] = virial[eV] / V * nktv2p with nktv2p = 1/(eV/Å³ per bar).
constexpr double kBarToEvA3 = 6.2415e-7;
constexpr double kNktv2p = 1.0 / kBarToEvA3;

struct LammpsFrame {
    double s = 1.0;
    double E = 0.0;
    double p[6] = {}; // bars, Voigt
    double V = 0.0;
    double lx = 0.0;
};

cpu::Frame make_nacl(double scale) {
    const int N = kSide * kSide * kSide;
    cpu::Frame frame(N);
    int idx = 0;
    for (int i = 0; i < kSide; ++i)
        for (int j = 0; j < kSide; ++j)
            for (int k = 0; k < kSide; ++k) {
                frame.span(ATV::kAtom_x)[idx] = (i * kD) * scale;
                frame.span(ATV::kAtom_y)[idx] = (j * kD) * scale;
                frame.span(ATV::kAtom_z)[idx] = (k * kD) * scale;
                frame.span(ATV::kAtom_charge)[idx] =
                    ((i + j + k) % 2 == 0) ? 1.0 : -1.0;
                frame.span(ATV::kAtom_mass)[idx] =
                    ((i + j + k) % 2 == 0) ? 22.9898 : 35.453;
                frame.span(ATV::kAtom_type)[idx] =
                    ((i + j + k) % 2 == 0) ? 0 : 1;
                frame.span(ATV::kAtom_tag)[idx] = idx + 1;
                ++idx;
            }
    frame.simbox = cpu::Box::cubic(kSide * kD * scale);
    return frame;
}

void evaluate_electrostatics(cpu::Frame& frame) {
    cpu::Evaluate eval;
    eval.set_neighbor(cpu::md::NeighborList<>{kCutoff, /*skin=*/0.3});
    eval.on(evaluate::Phase::evaluate).add(cpu::CoulLong{kAlpha, kCutoff});
    eval.on(evaluate::Phase::evaluate)
        .add(cpu::KSpacePME{kAlpha, kSpacing, kOrder});
    eval.build();
    eval.run(frame);
}

std::vector<LammpsFrame> load_reference(const std::string& path) {
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("cannot open LAMMPS reference: " + path);
    std::vector<LammpsFrame> frames;
    std::string line;
    while (std::getline(in, line)) {
        if (line.find("E=") == std::string::npos)
            continue;
        LammpsFrame f;
        auto get = [&](const char* key) {
            const auto p = line.find(key);
            if (p == std::string::npos)
                throw std::runtime_error(std::string("missing key ") + key);
            return std::stod(
                line.substr(p + std::char_traits<char>::length(key)));
        };
        // Optional s= field
        if (line.find("s=") != std::string::npos)
            f.s = get("s=");
        f.E = get("E=");
        f.p[0] = get("pxx=");
        f.p[1] = get("pyy=");
        f.p[2] = get("pzz=");
        f.p[3] = get("pxy=");
        f.p[4] = get("pxz=");
        f.p[5] = get("pyz=");
        f.V = get("V=");
        f.lx = get("lx=");
        frames.push_back(f);
    }
    if (frames.size() != 3)
        throw std::runtime_error("expected 3 LAMMPS frames, got "
                                 + std::to_string(frames.size()));
    return frames;
}

// Locate reference file relative to CWD (build/benchmarks or source tree).
std::string find_reference() {
    const char* candidates[] = {
        "md/npt_pme/validation/reference/frames_lammps.txt",
        "../benchmarks/md/npt_pme/validation/reference/frames_lammps.txt",
        "benchmarks/md/npt_pme/validation/reference/frames_lammps.txt",
        "../../benchmarks/md/npt_pme/validation/reference/frames_lammps.txt",
    };
    for (const char* c : candidates) {
        std::ifstream in(c);
        if (in)
            return c;
    }
    return candidates[0];
}

} // namespace

BENCH_CASE(PmeNacl_LammpsStaticFrames_EnergyForceStressVolume) {
    const auto ref = load_reference(find_reference());
    const double scales[3] = {1.00, 0.99, 1.02};

    // Discretization budget: PME grid + real-space cutoff. Energy relative
    // 1e-4; virial absolute ~0.5 eV on ~95 eV diagonals (~0.5 %).
    constexpr double kERel = 1e-4;
    constexpr double kVirAbs = 0.5; // eV
    constexpr double kVolRel = 1e-12;
    constexpr double kForceAbs = 1e-8; // eV/Å — rock-salt forces vanish

    for (int fi = 0; fi < 3; ++fi) {
        auto frame = make_nacl(scales[fi]);
        evaluate_electrostatics(frame);

        const double E = frame.meta.get(ATV::kEnergy);
        const double V = frame.simbox.volume();
        const double lx = frame.simbox.length(0);
        const auto& S = frame.meta.get(ATV::kStress);

        BM_NEAR(E, ref[fi].E, kERel * std::max(1.0, std::fabs(ref[fi].E)))
            << "frame " << fi << " energy";
        BM_NEAR(V, ref[fi].V, kVolRel * ref[fi].V)
            << "frame " << fi << " volume";
        BM_NEAR(lx, ref[fi].lx, kVolRel * ref[fi].lx)
            << "frame " << fi << " lx";

        for (int c = 0; c < 6; ++c) {
            const double vir_ref = ref[fi].p[c] * ref[fi].V / kNktv2p;
            BM_NEAR(S[c], vir_ref, kVirAbs)
                << "frame " << fi << " virial Voigt " << c << " atv=" << S[c]
                << " lmp=" << vir_ref;
        }

        // Scalar pressure from virial only (T=0): P = tr(σ)/(3V) [eV/Å³].
        const double P_atv = (S[0] + S[1] + S[2]) / (3.0 * V);
        const double P_lmp =
            (ref[fi].p[0] + ref[fi].p[1] + ref[fi].p[2]) / 3.0 * kBarToEvA3;
        BM_NEAR(P_atv, P_lmp, 1e-4 * std::max(1e-3, std::fabs(P_lmp)))
            << "frame " << fi << " pressure";

        double fsum[3] = {0, 0, 0};
        double fmax = 0.0;
        for (int i = 0; i < frame.nlocal(); ++i) {
            const double fx = frame.span(ATV::kAtom_fx)[i];
            const double fy = frame.span(ATV::kAtom_fy)[i];
            const double fz = frame.span(ATV::kAtom_fz)[i];
            fsum[0] += fx;
            fsum[1] += fy;
            fsum[2] += fz;
            fmax = std::max(fmax,
                            std::max(std::fabs(fx),
                                     std::max(std::fabs(fy), std::fabs(fz))));
        }
        BM_NEAR(fsum[0], 0.0, kForceAbs) << "frame " << fi << " force sum x";
        BM_NEAR(fsum[1], 0.0, kForceAbs) << "frame " << fi << " force sum y";
        BM_NEAR(fsum[2], 0.0, kForceAbs) << "frame " << fi << " force sum z";
        BM_NEAR(fmax, 0.0, kForceAbs) << "frame " << fi << " max |F|";

        std::printf("  frame %d s=%.2f  E_atv=%.6f E_lmp=%.6f  "
                    "vir_xx_atv=%.4f vir_xx_lmp=%.4f  |F|_max=%.2e\n",
                    fi, scales[fi], E, ref[fi].E, S[0],
                    ref[fi].p[0] * ref[fi].V / kNktv2p, fmax);
    }
}

BENCH_MAIN("MD coul/long+pppm — NaCl PME LAMMPS static-frame parity");
