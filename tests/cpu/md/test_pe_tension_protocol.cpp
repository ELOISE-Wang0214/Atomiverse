/// @file tests/cpu/md/test_pe_tension_protocol.cpp
/// @brief Protocol-pinning gate for the `polyethylene_tension` PE
///        performance benchmark: each stage's integrator / ensemble /
///        timestep must keep matching the published MSU five-stage
///        protocol, and the trajectory/RDF CLI surface must stay
///        exposed. Static text checks over the benchmark source — no
///        LAMMPS run involved.

#include <gtest/gtest.h>

#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

const std::string kRoot = ATV_SOURCE_DIR;
const std::string kBenchSrc =
    kRoot + "/benchmarks/md/pe/performance/polyethylene_tension.cpp";

/// @brief Read an entire file into a string ("" if missing).
std::string slurp(const std::string& path) {
    std::ifstream in(path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

} // namespace

// ─────────────────────────────────────────────────────────────────────
// polyethylene_tension benchmark — trajectory/RDF CLI surface  (ac-004)
// ─────────────────────────────────────────────────────────────────────

TEST(PeTensionBenchmark, ExposesTrajEveryAndRdfFlags) {
    const std::string src =
        slurp(kBenchSrc);
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
        slurp(kBenchSrc);
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
        slurp(kBenchSrc);
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
        slurp(kBenchSrc);
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
        slurp(kBenchSrc);
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
        slurp(kBenchSrc);
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
        slurp(kBenchSrc);
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
