/// @file tests/cpu/md/test_msu_reference_inputs.cpp
/// @brief Phase D.1 — well-formedness gate for the MSU amorphous-PE
///        LAMMPS reference deck, the `polyethylene_tension` benchmark's
///        trajectory/RDF CLI surface, and the format of any committed
///        reference CSV artifacts.
///
/// These are static text checks: they verify the inputs that the
/// 24h-class operator run (.claude/specs/msu-pe-phase-d-01-reference-runs)
/// consumes are present and faithful, without running LAMMPS or the
/// 360000-step protocol. The committed-artifact tests skip cleanly when
/// the operator has not yet produced the reference CSVs.
///
/// @see .claude/specs/msu-pe-phase-d-01-reference-runs.md

#include <gtest/gtest.h>

#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

const std::string kRoot   = ATV_SOURCE_DIR;
const std::string kRefDir = kRoot + "/benchmarks/md/reference/polyethylene";

/// @brief Read an entire file into a string ("" if missing).
std::string slurp(const std::string& path) {
    std::ifstream in(path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

/// @brief Split text into lines.
std::vector<std::string> lines_of(const std::string& text) {
    std::vector<std::string> out;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) out.push_back(line);
    return out;
}

/// @brief Strip leading whitespace.
std::string ltrim(const std::string& s) {
    size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    return s.substr(i);
}

/// @brief Whitespace-tokenize a line.
std::vector<std::string> tokens_of(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream in(line);
    std::string tok;
    while (in >> tok) out.push_back(tok);
    return out;
}

} // namespace

// ─────────────────────────────────────────────────────────────────────
// lammps.in — deck well-formedness  (ac-001 / ac-002 / ac-003)
// ─────────────────────────────────────────────────────────────────────

TEST(MsuLammpsDeck, ExistsAndReadsMsuData) {
    const std::string deck = slurp(kRefDir + "/lammps.in");
    ASSERT_FALSE(deck.empty()) << "lammps.in is missing or empty";
    EXPECT_NE(deck.find("read_data"), std::string::npos)
        << "deck has no read_data command";
    EXPECT_NE(deck.find("msu_pe.data"), std::string::npos)
        << "deck does not read msu_pe.data";
}

TEST(MsuLammpsDeck, HasFiveRunsWithCorrectStepCounts) {
    const std::string deck = slurp(kRefDir + "/lammps.in");
    ASSERT_FALSE(deck.empty());
    std::vector<long> runs;
    for (const auto& raw : lines_of(deck)) {
        const auto toks = tokens_of(ltrim(raw));
        if (toks.size() >= 2 && toks[0] == "run") {
            runs.push_back(std::stol(toks[1]));
        }
    }
    const std::vector<long> expect{10000, 50000, 50000, 50000, 200000};
    EXPECT_EQ(runs, expect)
        << "deck must run exactly the five MSU stages in order";
}

TEST(MsuLammpsDeck, PinsDreidingUaCoefficients) {
    const std::string deck = slurp(kRefDir + "/lammps.in");
    ASSERT_FALSE(deck.empty());
    EXPECT_NE(deck.find("bond_coeff"), std::string::npos);
    EXPECT_NE(deck.find("350 1.53"), std::string::npos)
        << "bond harmonic k/r0 must be 350 1.53";
    EXPECT_NE(deck.find("angle_coeff"), std::string::npos);
    EXPECT_NE(deck.find("60 109.5"), std::string::npos)
        << "angle harmonic k/theta0 must be 60 109.5";
    EXPECT_NE(deck.find("dihedral_coeff"), std::string::npos);
    EXPECT_NE(deck.find("1.73 -4.49 0.776 6.99 0.0"), std::string::npos)
        << "dihedral multi/harmonic A1..A5 must match SOURCES.md";
    EXPECT_NE(deck.find("pair_coeff"), std::string::npos);
    EXPECT_NE(deck.find("0.112 4.01"), std::string::npos)
        << "LJ eps/sigma must be 0.112 4.01";
    EXPECT_NE(deck.find("lj/cut 10.5"), std::string::npos)
        << "LJ cutoff must be 10.5";
}

TEST(MsuLammpsDeck, Stage3IsNvtCoolingStages2to5UseHalfFsTimestep) {
    const std::string deck = slurp(kRefDir + "/lammps.in");
    ASSERT_FALSE(deck.empty());
    // stage 3 is NVT with a 500 -> 100 K ramp, not NPT.
    EXPECT_NE(deck.find("nvt temp 500 100"), std::string::npos)
        << "stage 3 must be fix nvt with a temp 500 100 cooling ramp";
    // timestep budget: exactly one 1.0 fs line (stage 1) and at least
    // four 0.5 fs lines (stages 2-5).
    int n_half = 0, n_full = 0;
    for (const auto& raw : lines_of(deck)) {
        const auto toks = tokens_of(ltrim(raw));
        if (toks.size() >= 2 && toks[0] == "timestep") {
            const double dt = std::stod(toks[1]);
            if (dt == 0.5) ++n_half;
            else if (dt == 1.0) ++n_full;
        }
    }
    EXPECT_GE(n_half, 4) << "stages 2-5 must each use timestep 0.5";
    EXPECT_EQ(n_full, 1) << "only stage 1 uses timestep 1.0";
}

// ─────────────────────────────────────────────────────────────────────
// polyethylene_tension benchmark — trajectory/RDF CLI surface  (ac-004)
// ─────────────────────────────────────────────────────────────────────

TEST(PeTensionBenchmark, ExposesTrajEveryAndRdfFlags) {
    const std::string src =
        slurp(kRoot + "/benchmarks/md/polyethylene_tension.cpp");
    ASSERT_FALSE(src.empty()) << "polyethylene_tension.cpp not found";
    EXPECT_NE(src.find("--traj-every"), std::string::npos)
        << "benchmark must accept a --traj-every flag";
    EXPECT_NE(src.find("--rdf"), std::string::npos)
        << "benchmark must accept a --rdf flag";
}

// ─────────────────────────────────────────────────────────────────────
// Committed reference artifacts — column format  (ac-008)
// Skips cleanly when the operator has not yet run the protocol.
// ─────────────────────────────────────────────────────────────────────

namespace {

/// @brief Assert a present CSV's header and that its first column is a
///        strictly-increasing numeric series. No-op when the file is
///        absent (operator has not produced it yet).
void check_reference_csv(const std::string& path,
                         const std::string& want_header) {
    const std::string body = slurp(path);
    if (body.empty()) {
        GTEST_SKIP() << "reference artifact not present: " << path;
    }
    const auto rows = lines_of(body);
    ASSERT_GE(rows.size(), 3u)
        << path << " must have a header and a plausible row count";
    EXPECT_EQ(rows.front(), want_header)
        << path << " has the wrong header column set";
    double prev = -1e300;
    int data_rows = 0;
    for (size_t i = 1; i < rows.size(); ++i) {
        if (rows[i].empty()) continue;
        const double key = std::stod(rows[i].substr(0, rows[i].find(',')));
        EXPECT_GT(key, prev)
            << path << " row " << i << " first column not increasing";
        prev = key;
        ++data_rows;
    }
    EXPECT_GT(data_rows, 0) << path << " has no data rows";
}

} // namespace

TEST(MsuReferenceArtifacts, ThermoCsvsHaveCanonicalColumns) {
    const std::string header = "step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz";
    for (const char* eng : {"atomiverse", "lammps"}) {
        for (const char* rk : {"1r", "8r"}) {
            check_reference_csv(
                kRefDir + "/thermo_" + eng + "_" + rk + ".csv", header);
        }
    }
}

TEST(MsuReferenceArtifacts, StressStrainCsvsHaveCanonicalColumns) {
    const std::string header = "strain,sigma_xx";
    for (const char* eng : {"atomiverse", "lammps"}) {
        for (const char* rk : {"1r", "8r"}) {
            check_reference_csv(
                kRefDir + "/stress_strain_" + eng + "_" + rk + ".csv",
                header);
        }
    }
}
