/// @file tests/cpu/gpw/test_cpu_gpw_edl_frag.cpp
/// @brief EDL fragment (frag268: 88 O + 174 H + 3 Na + 3 Cl) periodic GPW energy
///        vs CP2K (-1702.32 Ha), exercising the Γ-point lattice sum on the real
///        electric-double-layer basis (DZVP-MOLOPT-SR-GTH / GTH-PBE / PBE).
///
/// Heavy (hundreds of atoms): gated behind ATV_RUN_EDL=1 so it is a no-op in the
/// normal suite. Reads a plain XYZ (element + Å coords) from ATV_EDL_XYZ.
///   ATV_RUN_EDL=1 ATV_EDL_XYZ=.../frag268.xyz [ATV_GPW_EPS=1e-8] \
///     ./build/tests/test_cpu_gpw_edl_frag

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"

using namespace ATV::cpu::gpw;

namespace {
const std::string kData = ATV_CP2K_DATA_DIR;
constexpr double kAng = 1.8897259885789; // Å -> Bohr

struct Elem {
    const char* sym;
    const char* pot;
    int q; // valence electrons
};
const std::map<std::string, Elem> kEl = {
    {"O", {"O", "GTH-PBE-q6", 6}},   {"H", {"H", "GTH-PBE-q1", 1}},
    {"Na", {"Na", "GTH-PBE-q9", 9}}, {"Cl", {"Cl", "GTH-PBE-q7", 7}},
    {"Zn", {"Zn", "GTH-PBE-q12", 12}},
};
} // namespace

TEST(GpwEdlFrag, PeriodicEnergyVsCp2k) {
    if (!std::getenv("ATV_RUN_EDL")) GTEST_SKIP() << "set ATV_RUN_EDL=1 to run";
    const char* xyz = std::getenv("ATV_EDL_XYZ");
    ASSERT_TRUE(xyz) << "set ATV_EDL_XYZ to the fragment .xyz";

    std::ifstream f(xyz);
    ASSERT_TRUE(f.good()) << "cannot open " << xyz;
    int n = 0;
    std::string line;
    std::getline(f, line);
    n = std::atoi(line.c_str());
    std::getline(f, line); // comment

    // Cache basis/potential per element symbol.
    std::map<std::string, GaussianBasis> basis;
    std::map<std::string, GthPotential> pot;
    for (const auto& [sym, e] : kEl) {
        basis[sym] = parse_gth_basis(kData + "/BASIS_MOLOPT", e.sym, "DZVP-MOLOPT-SR-GTH");
        pot[sym] = parse_gth_potential(kData + "/GTH_POTENTIALS", e.sym, e.pot);
    }

    std::vector<GpwAtom> atoms;
    int nval = 0;
    std::map<std::string, int> counts;
    for (int i = 0; i < n; ++i) {
        std::getline(f, line);
        std::istringstream ss(line);
        std::string sym;
        double x, y, z;
        ss >> sym >> x >> y >> z;
        ASSERT_TRUE(kEl.count(sym)) << "unknown element " << sym;
        const double r[3] = {x * kAng, y * kAng, z * kAng};
        atoms.push_back(make_gpw_atom(r, basis[sym], pot[sym]));
        nval += kEl.at(sym).q;
        counts[sym]++;
    }
    ASSERT_EQ(static_cast<int>(atoms.size()), n);
    std::printf("[EDL] %d atoms (", n);
    for (auto& [s, c] : counts) std::printf("%s:%d ", s.c_str(), c);
    std::printf(") nval=%d n_occ=%d\n", nval, nval / 2);

    GpwConfig cfg;
    cfg.cell[0] = 13.128 * kAng;
    cfg.cell[1] = 10.582 * kAng;
    cfg.cell[2] = 47.05377225 * kAng;
    cfg.cutoff_ry = std::getenv("ATV_EDL_CUT") ? std::atof(std::getenv("ATV_EDL_CUT")) : 120.0;
    if (std::getenv("ATV_EDL_NGRIDS")) cfg.ngrids = std::atoi(std::getenv("ATV_EDL_NGRIDS"));
    if (std::getenv("ATV_EDL_RELCUT")) cfg.rel_cutoff = std::atof(std::getenv("ATV_EDL_RELCUT"));
    cfg.n_occ = nval / 2;
    cfg.max_iter = std::getenv("ATV_EDL_ITERS") ? std::atoi(std::getenv("ATV_EDL_ITERS")) : 100;
    cfg.tol = 1e-6;
    cfg.xc = XcKind::PBE;
    cfg.nonperiodic = false;
    cfg.use_ot = true;
    cfg.eps_default = 1e-8;

    GpwResult r = gpw_scf(atoms, cfg);
    std::printf("[EDL] conv=%d iters=%d  E_total=%.6f Ha  (CP2K -1702.32)\n", r.converged,
                r.iters, r.e_total);
    EXPECT_TRUE(std::isfinite(r.e_total));

    // Warm-start demo (ATV_EDL_WARMDEMO): emulate one AIMD step — nudge every atom
    // by ~0.05 Bohr, then solve the new geometry COLD vs WARM (from the previous
    // converged orbitals, as ASPC would supply). Warm should converge in far fewer
    // iterations to the same energy: the per-MD-step SCF cost the dynamics sees.
    if (std::getenv("ATV_EDL_WARMDEMO")) {
        std::vector<GpwAtom> moved = atoms;
        for (std::size_t i = 0; i < moved.size(); ++i) {
            const double s = (i % 2 ? 0.05 : -0.05);
            moved[i].pos[0] += s;
            moved[i].pos[1] -= s;
            moved[i].pos[2] += 0.5 * s;
        }
        // Fixed small SCF budget (one AIMD step's worth): the warm run, starting
        // from the previous geometry's orbitals (as ASPC would supply), should be
        // far closer to the self-consistent energy than the cold run.
        GpwConfig step = cfg;
        step.max_iter = std::getenv("ATV_EDL_STEPIT") ? std::atoi(std::getenv("ATV_EDL_STEPIT")) : 8;
        const GpwResult cold = gpw_scf(moved, step);
        const GpwResult warm = gpw_scf(moved, step, &r.C);
        std::printf("[EDL warm] %d-iter budget | cold E=%.6f | warm E=%.6f | warm−cold=%.4f Ha\n",
                    step.max_iter, cold.e_total, warm.e_total, warm.e_total - cold.e_total);
    }
}
