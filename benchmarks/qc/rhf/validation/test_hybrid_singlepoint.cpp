/// @file tests/cpu/scf/test_hybrid_singlepoint.cpp
/// @brief B3LYP / PBE0 global-hybrid single-point parity against PySCF.
///
/// Acceptance criteria for spec `phase3-hybrid-dft-02-hybrid-singlepoint`:
///   - ac-001: B3LYP total energy vs PySCF, <= 1e-7 Ha (5 molecules).
///   - ac-002: PBE0 total energy vs PySCF, <= 1e-7 Ha (5 molecules).
///   - ac-003: HOMO / LUMO orbital energies vs PySCF, <= 1e-5 eV.
///   - ac-004: alpha_hf=0 hybrid path reproduces pure-PBE bit-for-bit.
///   - ac-005: energy decomposition closes to 1e-10 Ha.
///   - ac-006: b3lyp_config / pbe0_config set xc_type + alpha_hf.
///
/// Ground truth is genuine PySCF 2.13.0 `RKS(mol, xc=...)` at 6-31G*,
/// checked into the per-molecule JSON files under `benchmarks/qc/rhf/validation/reference/`
/// (keyed `b3lyp` / `pbe0` blocks, each carrying `pyscf_version`). The
/// geometries below are copied verbatim from `benchmarks/qc/rhf/validation/molecules/*.xyz`,
/// the exact geometries those PySCF runs used.

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/scf/scf_solver.h"
#include "atomiverse/cpu/scf/xc_functional.h"

using namespace ATV;
using ATV::cpu::scf::BasisSet;
using ATV::cpu::scf::XCType;

namespace {

constexpr double kHartreeToEV = 27.211386245988;

// ── Minimal JSON value extraction ───────────────────────────────────────
//
// The reference JSONs are pinned from a one-off PySCF run (generator removed) with
// a stable, flat-per-block shape. Rather than pull in a JSON dependency,
// we extract one numeric value from one keyed block. The block is found by
// `"<block>": {` and the key by `"<key>":` inside it.

std::string read_reference_file(const std::string& stem) {
    const std::string path = std::string(ATV_SOURCE_DIR)
                             + "/benchmarks/qc/rhf/validation/reference/" + stem + ".json";
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("cannot open reference JSON: " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

double json_value(const std::string& doc, const std::string& block,
                  const std::string& key) {
    const std::string blk = "\"" + block + "\"";
    std::size_t bpos = doc.find(blk);
    if (bpos == std::string::npos)
        throw std::runtime_error("block not found: " + block);
    std::size_t brace = doc.find('{', bpos);
    if (brace == std::string::npos)
        throw std::runtime_error("malformed block: " + block);
    std::size_t bend = doc.find('}', brace);
    const std::string region = doc.substr(brace, bend - brace);

    const std::string kq = "\"" + key + "\"";
    std::size_t kpos = region.find(kq);
    if (kpos == std::string::npos)
        throw std::runtime_error("key not found: " + block + "/" + key);
    std::size_t colon = region.find(':', kpos);
    std::size_t vstart = region.find_first_of("-0123456789.eE", colon);
    std::size_t vend = region.find_first_of(",}\n", vstart);
    return std::stod(region.substr(vstart, vend - vstart));
}

// ── Molecule geometries (verbatim from benchmarks/qc/rhf/validation/molecules/*.xyz) ────

cpu::Frame make_h2o() {
    cpu::Frame f(3);
    f.atoms.type = {8, 1, 1};
    f.atoms.x = {0.000000, 0.000000, 0.000000};
    f.atoms.y = {0.000000, 0.755453, -0.755453};
    f.atoms.z = {0.117790, -0.471161, -0.471161};
    f.box = Mat3::cubic(100.0);
    return f;
}

cpu::Frame make_ch4() {
    cpu::Frame f(5);
    f.atoms.type = {6, 1, 1, 1, 1};
    f.atoms.x = {0.000000, 0.627640, -0.627640, 0.627640, -0.627640};
    f.atoms.y = {0.000000, 0.627640, -0.627640, -0.627640, 0.627640};
    f.atoms.z = {0.000000, 0.627640, 0.627640, -0.627640, -0.627640};
    f.box = Mat3::cubic(100.0);
    return f;
}

cpu::Frame make_nh3() {
    cpu::Frame f(4);
    f.atoms.type = {7, 1, 1, 1};
    f.atoms.x = {0.000000, 0.000000, 0.812328, -0.812328};
    f.atoms.y = {0.000000, 0.938035, -0.469017, -0.469017};
    f.atoms.z = {0.113880, -0.265825, -0.265825, -0.265825};
    f.box = Mat3::cubic(100.0);
    return f;
}

cpu::Frame make_benzene() {
    cpu::Frame f(12);
    f.atoms.type = {6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1};
    f.atoms.x = {-1.362315, -0.956919, 0.405396, 1.362315, 0.956919, -0.405396,
                 -2.416240, -1.697218, 0.719021, 2.416240, 1.697218, -0.719021};
    f.atoms.y = {-0.318457, 1.020526, 1.338982, 0.318457, -1.020526, -1.338982,
                 -0.564824, 1.810033, 2.374857, 0.564824, -1.810033, -2.374857};
    f.atoms.z = {0.002346, 0.012443, 0.010097, -0.002345, -0.012442, -0.010097,
                 0.004159, 0.022066, 0.017908, -0.004160, -0.022067, -0.017909};
    f.box = Mat3::cubic(100.0);
    return f;
}

cpu::Frame make_glycine() {
    cpu::Frame f(10);
    f.atoms.type = {7, 6, 6, 8, 8, 1, 1, 1, 1, 1};
    f.atoms.x = {-1.400839, -0.370019, 0.925792,  1.233153,  1.771618,
                 -1.625527, -2.268277, -0.238214, -0.673752, 2.646065};
    f.atoms.y = {0.271160, -0.558813, 0.184484,  0.868576,  0.156585,
                 1.052674, -0.307372, -1.486583, -0.848738, 0.668026};
    f.atoms.z = {-0.425197, 0.187422,  0.223822,  1.236362, -0.880405,
                 0.233297,  -0.505344, -0.410315, 1.218266, -0.877908};
    f.box = Mat3::cubic(100.0);
    return f;
}

cpu::Frame make_molecule(const std::string& stem) {
    if (stem == "h2o")
        return make_h2o();
    if (stem == "ch4")
        return make_ch4();
    if (stem == "nh3")
        return make_nh3();
    if (stem == "benzene")
        return make_benzene();
    if (stem == "glycine")
        return make_glycine();
    throw std::runtime_error("unknown molecule: " + stem);
}

const char* kMolecules[] = {"h2o", "ch4", "nh3", "benzene", "glycine"};

// Run a hybrid single-point and return the converged result.
cpu::scf::SCFResult run_hybrid_uncached(const std::string& stem, XCType xc) {
    auto frame = make_molecule(stem);
    auto integrals = cpu::scf::build_integrals(frame, BasisSet::Pople631Gs);
    cpu::scf::SCFConfig cfg = (xc == XCType::B3LYP)
                                  ? cpu::scf::b3lyp_config(integrals)
                                  : cpu::scf::pbe0_config(integrals);
    cfg.tol = 1e-10;
    cpu::scf::SCFSolver solver{cfg};
    solver.build(integrals);
    return solver.run(integrals);
}

// Hybrid SCF is K-heavy (~2x pure PBE). The energy, HOMO/LUMO and
// decomposition criteria all read the SAME converged result, so each
// (molecule, functional) SCF is run once and cached for the whole test
// binary -- avoiding ~3x redundant benzene/glycine runs.
const cpu::scf::SCFResult& run_hybrid(const std::string& stem, XCType xc) {
    static std::map<std::pair<std::string, int>, cpu::scf::SCFResult> cache;
    const auto key = std::make_pair(stem, static_cast<int>(xc));
    auto it = cache.find(key);
    if (it == cache.end())
        it = cache.emplace(key, run_hybrid_uncached(stem, xc)).first;
    return it->second;
}

} // namespace

// ── ac-006: config helpers ──────────────────────────────────────────────

TEST(HybridSinglePoint, ConfigHelpersSetFunctionalAndAlpha) {
    auto frame = make_h2o();
    auto integrals = cpu::scf::build_integrals(frame, BasisSet::Pople631Gs);

    auto b3 = cpu::scf::b3lyp_config(integrals);
    EXPECT_EQ(b3.method, cpu::scf::MeanFieldType::KS);
    EXPECT_EQ(b3.spin, cpu::scf::SpinType::Restricted);
    EXPECT_EQ(b3.xc_type, XCType::B3LYP);
    EXPECT_DOUBLE_EQ(cpu::scf::hybrid_exchange_fraction(XCType::B3LYP), 0.20);

    auto pbe0 = cpu::scf::pbe0_config(integrals);
    EXPECT_EQ(pbe0.method, cpu::scf::MeanFieldType::KS);
    EXPECT_EQ(pbe0.xc_type, XCType::PBE0);
    EXPECT_DOUBLE_EQ(cpu::scf::hybrid_exchange_fraction(XCType::PBE0), 0.25);

    // Non-hybrid functionals carry zero exact-exchange fraction.
    EXPECT_DOUBLE_EQ(cpu::scf::hybrid_exchange_fraction(XCType::PBE), 0.0);
    EXPECT_DOUBLE_EQ(cpu::scf::hybrid_exchange_fraction(XCType::LSDA), 0.0);

    // is_gga predicate (reused by link 3).
    EXPECT_TRUE(cpu::scf::is_gga(XCType::PBE));
    EXPECT_TRUE(cpu::scf::is_gga(XCType::B3LYP));
    EXPECT_TRUE(cpu::scf::is_gga(XCType::PBE0));
    EXPECT_FALSE(cpu::scf::is_gga(XCType::LSDA));
}

// ── ac-001: B3LYP total energy vs PySCF ─────────────────────────────────

TEST(HybridSinglePoint, B3LYPEnergyMatchesPySCF) {
    for (const char* stem : kMolecules) {
        const std::string doc = read_reference_file(stem);
        const double e_ref = json_value(doc, "b3lyp", "energy_hartree");
        auto result = run_hybrid(stem, XCType::B3LYP);
        ASSERT_TRUE(result.converged) << "B3LYP did not converge: " << stem;
        EXPECT_NEAR(result.energy, e_ref, 1e-7)
            << "B3LYP/6-31G* energy mismatch for " << stem;
    }
}

// ── ac-002: PBE0 total energy vs PySCF ──────────────────────────────────

TEST(HybridSinglePoint, PBE0EnergyMatchesPySCF) {
    for (const char* stem : kMolecules) {
        const std::string doc = read_reference_file(stem);
        const double e_ref = json_value(doc, "pbe0", "energy_hartree");
        auto result = run_hybrid(stem, XCType::PBE0);
        ASSERT_TRUE(result.converged) << "PBE0 did not converge: " << stem;
        EXPECT_NEAR(result.energy, e_ref, 1e-7)
            << "PBE0/6-31G* energy mismatch for " << stem;
    }
}

// ── ac-003: HOMO / LUMO orbital energies vs PySCF ───────────────────────

TEST(HybridSinglePoint, HomoLumoMatchesPySCF) {
    for (const char* stem : kMolecules) {
        const std::string doc = read_reference_file(stem);
        for (const char* xc_name : {"b3lyp", "pbe0"}) {
            const XCType xc = (std::string(xc_name) == "b3lyp") ? XCType::B3LYP
                                                                : XCType::PBE0;
            const double homo_ref = json_value(doc, xc_name, "homo_hartree");
            const double lumo_ref = json_value(doc, xc_name, "lumo_hartree");

            auto result = run_hybrid(stem, xc);
            ASSERT_TRUE(result.converged) << xc_name << " " << stem;

            const int n_occ = result.n_alpha;
            const double homo = result.eigenvalues_alpha(n_occ - 1);
            const double lumo = result.eigenvalues_alpha(n_occ);

            EXPECT_NEAR(homo * kHartreeToEV, homo_ref * kHartreeToEV, 1e-5)
                << "HOMO mismatch " << xc_name << " " << stem;
            EXPECT_NEAR(lumo * kHartreeToEV, lumo_ref * kHartreeToEV, 1e-5)
                << "LUMO mismatch " << xc_name << " " << stem;
        }
    }
}

// ── ac-005: energy decomposition closes ─────────────────────────────────

TEST(HybridSinglePoint, EnergyDecompositionCloses) {
    for (const char* stem : kMolecules) {
        for (XCType xc : {XCType::B3LYP, XCType::PBE0}) {
            auto result = run_hybrid(stem, xc);
            ASSERT_TRUE(result.converged);
            const double recomposed =
                result.one_electron_energy + result.coulomb_energy
                + result.exchange_energy + result.xc_energy
                + result.nuclear_repulsion;
            EXPECT_NEAR(result.energy, recomposed, 1e-10)
                << "decomposition did not close for " << stem;
        }
    }
}

// ── ac-004: alpha_hf=0 reproduces pure-PBE bit-for-bit ──────────────────
//
// Driving the PBE0 hybrid path (KS + Coulomb + K + Vxc^sl) with the
// exact-exchange fraction forced to 0 must collapse exactly onto the
// pre-existing pure-PBE SCF result. The PBE0 semi-local part is
// 0.75*PBE-X + 1.00*PBE-C, so the comparison uses a custom XCType whose
// semi-local mix is pure PBE; here we instead compare a PBE run against a
// PBE0 run with alpha forced to 0 AND the semi-local exchange coefficient
// restored to 1.0 -- which is exactly pure PBE. The cleanest invariant
// the implementation must satisfy: a hybrid Fock build with alpha_hf=0
// equals the non-hybrid Fock build for the same semi-local functional.
//
// We test the load-bearing convention guarantee directly: run pure PBE
// through the standard path, and run the hybrid code path with alpha_hf
// monkey-patched to 0 via an internal hook; both must be bit-identical.

TEST(HybridSinglePoint, AlphaZeroReproducesPurePBE) {
    // Bit-for-bit equality requires a deterministic reduction order: the
    // direct J / K builders use OpenMP thread-local accumulation whose
    // final sum is non-associative and therefore thread-count dependent.
    // Pin to one thread so the two SCF runs below execute the identical
    // floating-point operation sequence -- then alpha_hf=0 collapsing the
    // hybrid path onto pure PBE is a genuine bitwise check.
#ifdef _OPENMP
    const int saved_threads = omp_get_max_threads();
    omp_set_num_threads(1);
#endif

    auto frame = make_h2o();
    auto integrals = cpu::scf::build_integrals(frame, BasisSet::Pople631Gs);

    // Pure PBE through default_config (non-hybrid KS path).
    auto pbe_cfg = cpu::scf::default_config(integrals);
    pbe_cfg.xc_type = XCType::PBE;
    pbe_cfg.tol = 1e-10;
    cpu::scf::SCFSolver pbe_solver{pbe_cfg};
    pbe_solver.build(integrals);
    const double e_pbe = pbe_solver.run(integrals).energy;

    // Hybrid code path with alpha_hf forced to 0. With alpha_hf=0 the
    // K term and exact-exchange energy both vanish; the semi-local
    // functional must remain pure PBE. SCFConfig exposes
    // `hybrid_alpha_override` for exactly this regression.
    auto hyb_cfg = cpu::scf::default_config(integrals);
    hyb_cfg.xc_type = XCType::PBE;
    hyb_cfg.tol = 1e-10;
    hyb_cfg.hybrid_alpha_override = 0.0; // force the K branch with alpha=0
    cpu::scf::SCFSolver hyb_solver{hyb_cfg};
    hyb_solver.build(integrals);
    const double e_hyb = hyb_solver.run(integrals).energy;

#ifdef _OPENMP
    omp_set_num_threads(saved_threads);
#endif

    EXPECT_EQ(e_pbe, e_hyb)
        << "alpha_hf=0 hybrid path must reproduce pure-PBE bit-for-bit";
}

#endif // ATV_HAS_RHF
