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

const std::string kRoot = ATV_SOURCE_DIR;
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
    while (std::getline(in, line))
        out.push_back(line);
    return out;
}

/// @brief Strip leading whitespace.
std::string ltrim(const std::string& s) {
    size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i])))
        ++i;
    return s.substr(i);
}

/// @brief Whitespace-tokenize a line.
std::vector<std::string> tokens_of(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream in(line);
    std::string tok;
    while (in >> tok)
        out.push_back(tok);
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
    // Step counts are byte-for-byte the canonical lammps_msu_paper.in.
    const std::vector<long> expect{10000, 50000, 50000, 50000, 171800};
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

TEST(MsuLammpsDeck, MatchesCanonicalProtocolEnsembles) {
    // The deck's physics must be byte-for-byte the canonical upstream
    // protocol (lammps_msu_paper.in): stage 1 is the overlap-relaxation
    // stage (nve/limit + Langevin), stages 2-4 are NPT, stage 3 carries
    // the 500 -> 100 K cooling ramp, stage 5 is NPT(y,z) + deform(x).
    const std::string deck = slurp(kRefDir + "/lammps.in");
    ASSERT_FALSE(deck.empty());
    EXPECT_NE(deck.find("nve/limit"), std::string::npos)
        << "stage 1 must use fix nve/limit (overlap relaxation)";
    EXPECT_NE(deck.find("langevin"), std::string::npos)
        << "stage 1 must use a Langevin thermostat";
    EXPECT_NE(deck.find("npt temp 500 100"), std::string::npos)
        << "stage 3 must be fix npt with a temp 500 100 cooling ramp";
    EXPECT_NE(deck.find("deform 1 x erate"), std::string::npos)
        << "stage 5 must drive the box with fix deform along x";
    // dt budget: stages 1 and 5 use 1.0 fs, stages 2-4 use 0.5 fs.
    int n_half = 0, n_full = 0;
    for (const auto& raw : lines_of(deck)) {
        const auto toks = tokens_of(ltrim(raw));
        if (toks.size() >= 2 && toks[0] == "timestep") {
            const double dt = std::stod(toks[1]);
            if (dt == 0.5)
                ++n_half;
            else if (dt == 1.0)
                ++n_full;
        }
    }
    EXPECT_EQ(n_half, 3) << "stages 2-4 use timestep 0.5";
    EXPECT_EQ(n_full, 2) << "stages 1 and 5 use timestep 1.0";
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
// polyethylene_tension benchmark — per-stage protocol mirrors lammps.in
// (ac-006: each stage's integrator / ensemble / timestep matches the
//  corresponding lammps.in stage).
// ─────────────────────────────────────────────────────────────────────

namespace {

/// @brief Helper section accessor: returns the substring of `src`
///        between the `// ── Stage <id> ` marker and the next stage
///        marker (or end of file). Empty when the marker is absent.
std::string stage_section(const std::string& src, int id) {
    const std::string marker = "// ── Stage " + std::to_string(id) + " ";
    const size_t start = src.find(marker);
    if (start == std::string::npos)
        return "";
    const std::string next = "// ── Stage " + std::to_string(id + 1) + " ";
    const size_t end = src.find(next, start);
    return src.substr(start, end == std::string::npos ? std::string::npos
                                                      : end - start);
}

} // namespace

TEST(PeTensionBenchmark, BuildsAFreshDriverPerStage) {
    // ac-001 — each stage constructs and build()s its own cpu::md::MD.
    const std::string src =
        slurp(kRoot + "/benchmarks/md/polyethylene_tension.cpp");
    ASSERT_FALSE(src.empty());
    // There must be a per-stage construction helper that builds a driver.
    EXPECT_NE(src.find("cpu::md::MD"), std::string::npos)
        << "benchmark must construct cpu::md::MD drivers";
    // The obsolete single-dt-collapse comment and --cool-segments flag
    // must be gone.
    EXPECT_EQ(src.find("--cool-segments"), std::string::npos)
        << "--cool-segments flag must be removed (stage 3 is a single ramp)";
    EXPECT_EQ(src.find("cool_segments"), std::string::npos)
        << "cool_segments micro-segment loop must be removed";
    // The obsolete single-dt-collapse comment ("dt is uniform within a
    // stage block") must be gone — each stage now builds at its own dt.
    EXPECT_EQ(src.find("dt is uniform within a"), std::string::npos)
        << "obsolete 'dt is uniform within a stage block' comment must be gone";
}

TEST(PeTensionBenchmark, Stage1UsesNveLimitPlusLangevinAtDt1) {
    // ac-002 — stage 1 mirrors lammps.in `nve/limit 0.05` + `langevin 500`.
    const std::string src =
        slurp(kRoot + "/benchmarks/md/polyethylene_tension.cpp");
    const std::string s1 = stage_section(src, 1);
    ASSERT_FALSE(s1.empty()) << "stage 1 section marker not found";
    EXPECT_NE(s1.find("NVELimit"), std::string::npos)
        << "stage 1 must use NVELimit (lammps.in `fix nve/limit`)";
    EXPECT_NE(s1.find("FixLangevin"), std::string::npos)
        << "stage 1 must use FixLangevin (lammps.in `fix langevin`)";
    EXPECT_NE(s1.find("InitVelocity"), std::string::npos)
        << "stage 1 must seed velocities (lammps.in `velocity create`)";
    EXPECT_NE(s1.find("dt_stage1"), std::string::npos)
        << "stage 1 must run at dt_stage1 (1.0 fs)";
}

TEST(PeTensionBenchmark, Stages2And4UseIsotropicNoseHooverNpt) {
    // ac-004 — stages 2 and 4 are isotropic NoseHooverNPT (T_stop==T_start)
    // at dt 0.5 fs targeting P=0.
    const std::string src =
        slurp(kRoot + "/benchmarks/md/polyethylene_tension.cpp");
    for (int id : {2, 4}) {
        const std::string s = stage_section(src, id);
        ASSERT_FALSE(s.empty()) << "stage " << id << " section not found";
        EXPECT_NE(s.find("NoseHooverNPT"), std::string::npos)
            << "stage " << id << " must use NoseHooverNPT";
        EXPECT_NE(s.find("dt_stage" + std::to_string(id)), std::string::npos)
            << "stage " << id << " must run at its own dt (0.5 fs)";
        // Isotropic: no couple_axis override in stages 2/4.
        EXPECT_EQ(s.find("couple_axis"), std::string::npos)
            << "stage " << id << " must be isotropic (no couple_axis override)";
    }
}

TEST(PeTensionBenchmark, Stage3IsSingleRampedNoseHooverNpt) {
    // ac-003 — stage 3 is one NoseHooverNPT with T_start=500, T_stop=100
    // and a set_run_span, not a cool-segment loop.
    const std::string src =
        slurp(kRoot + "/benchmarks/md/polyethylene_tension.cpp");
    const std::string s3 = stage_section(src, 3);
    ASSERT_FALSE(s3.empty()) << "stage 3 section not found";
    EXPECT_NE(s3.find("NoseHooverNPT"), std::string::npos)
        << "stage 3 must use NoseHooverNPT";
    EXPECT_NE(s3.find("T_stop"), std::string::npos)
        << "stage 3 must set T_stop for the 500->100 K ramp";
    EXPECT_NE(s3.find("set_run_span"), std::string::npos)
        << "stage 3 must call set_run_span for the temperature ramp";
    EXPECT_EQ(s3.find("for ("), std::string::npos)
        << "stage 3 must run in one sim.run, not a micro-segment loop";
}

TEST(PeTensionBenchmark, Stage5IsAnisotropicNoseHooverNptPlusDeform) {
    // ac-005 — stage 5 is NoseHooverNPT couple_axis={false,true,true}
    // plus Deform on axis x at dt 1.0 fs.
    const std::string src =
        slurp(kRoot + "/benchmarks/md/polyethylene_tension.cpp");
    const std::string s5 = stage_section(src, 5);
    ASSERT_FALSE(s5.empty()) << "stage 5 section not found";
    EXPECT_NE(s5.find("NoseHooverNPT"), std::string::npos)
        << "stage 5 must use NoseHooverNPT";
    EXPECT_NE(s5.find("couple_axis"), std::string::npos)
        << "stage 5 must set couple_axis (x rigid, y/z barostatted)";
    EXPECT_NE(s5.find("{false, true, true}"), std::string::npos)
        << "stage 5 couple_axis must be {false,true,true}";
    EXPECT_NE(s5.find("Deform"), std::string::npos)
        << "stage 5 must use Deform along x";
    EXPECT_NE(s5.find("dt_stage5"), std::string::npos)
        << "stage 5 must run at dt_stage5 (1.0 fs)";
}

TEST(PeTensionBenchmark, PerStageTimestepsMatchLammpsIn) {
    // ac-006 — the per-stage dt defaults are 1.0/0.5/0.5/0.5/1.0 fs,
    // byte-aligned with the lammps.in `timestep` budget.
    const std::string src =
        slurp(kRoot + "/benchmarks/md/polyethylene_tension.cpp");
    ASSERT_FALSE(src.empty());
    EXPECT_NE(src.find("dt_stage1"), std::string::npos);
    EXPECT_NE(src.find("dt_stage2"), std::string::npos);
    EXPECT_NE(src.find("dt_stage3"), std::string::npos);
    EXPECT_NE(src.find("dt_stage4"), std::string::npos);
    EXPECT_NE(src.find("dt_stage5"), std::string::npos);
    // Default values: stages 1 and 5 at 1.0 fs, stages 2-4 at 0.5 fs.
    EXPECT_NE(src.find("dt_stage1 = 1.0"), std::string::npos)
        << "stage 1 default dt must be 1.0 fs";
    EXPECT_NE(src.find("dt_stage2 = 0.5"), std::string::npos)
        << "stage 2 default dt must be 0.5 fs";
    EXPECT_NE(src.find("dt_stage3 = 0.5"), std::string::npos)
        << "stage 3 default dt must be 0.5 fs";
    EXPECT_NE(src.find("dt_stage4 = 0.5"), std::string::npos)
        << "stage 4 default dt must be 0.5 fs";
    EXPECT_NE(src.find("dt_stage5 = 1.0"), std::string::npos)
        << "stage 5 default dt must be 1.0 fs";
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
        if (rows[i].empty())
            continue;
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
            check_reference_csv(kRefDir + "/thermo_" + eng + "_" + rk + ".csv",
                                header);
        }
    }
}

TEST(MsuReferenceArtifacts, StressStrainCsvsHaveCanonicalColumns) {
    const std::string header = "strain,sigma_xx";
    for (const char* eng : {"atomiverse", "lammps"}) {
        for (const char* rk : {"1r", "8r"}) {
            check_reference_csv(
                kRefDir + "/stress_strain_" + eng + "_" + rk + ".csv", header);
        }
    }
}
