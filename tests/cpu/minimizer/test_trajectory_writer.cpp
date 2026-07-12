/// @file tests/cpu/minimizer/test_trajectory_writer.cpp
/// @brief Tests for cpu::minimizer::TrajectoryWriter and OptLog + FinalGeometry.
///
/// These three components share the same driver binding style, so they are
/// exercised together: an H2 L-BFGS run writes a multi-frame trajectory, a
/// TSV log, and a single-frame final geometry.

#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/minimizer/final_geometry.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/minimizer/opt_log.h"
#include "atomiverse/cpu/minimizer/trajectory_writer.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;

namespace {

cpu::Frame make_h2(double bond_length) {
    cpu::Frame f(2);
    f.atoms.x = {0.0, bond_length};
    f.atoms.y = {0.0, 0.0};
    f.atoms.z = {0.0, 0.0};
    f.atoms.mass = {1.008, 1.008};
    f.atoms.type = {1, 1};
    f.box = Mat3::cubic(100.0);
    return f;
}

int count_lines(const std::string& path) {
    std::ifstream is(path);
    if (!is)
        return -1;
    int n = 0;
    std::string line;
    while (std::getline(is, line))
        ++n;
    return n;
}

} // namespace

class MinimizerRecordersTest : public ::testing::Test {
  protected:
    std::filesystem::path tmpdir_;

    void SetUp() override {
        tmpdir_ = std::filesystem::temp_directory_path() / "atv_min_rec_test";
        std::filesystem::create_directories(tmpdir_);
    }
    void TearDown() override {
        std::error_code ec;
        std::filesystem::remove_all(tmpdir_, ec);
    }
};

TEST_F(MinimizerRecordersTest, TrajectoryFileIsAppendedEveryStep) {
    auto frame = make_h2(1.2);
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    const auto traj_path = tmpdir_ / "traj.xyz";
    cpu::minimizer::Minimize opt;
    opt.add(cpu::minimizer::LBFGS{5, 0.01, 10});
    opt.add_handler(cpu::ForceConverged{/*ftol=*/0.5});
    opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
    opt.on(minimizer::Phase::end_of_step)
        .add(
            cpu::minimizer::TrajectoryWriter{traj_path.string(), /*stride=*/1});
    opt.build();
    opt.run(frame, 4);

    // Each XYZ frame is `n_atoms + 2` lines (header + comment + per-atom).
    const int n_atoms = 2;
    int n_lines = count_lines(traj_path.string());
    ASSERT_GT(n_lines, 0) << "trajectory file not written";
    int n_frames = n_lines / (n_atoms + 2);
    EXPECT_EQ(n_frames, 4) << "expected 4 frames (1 per step) in " << traj_path;
}

TEST_F(MinimizerRecordersTest, TrajectoryRespectsStride) {
    auto frame = make_h2(1.2);
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    const auto traj_path = tmpdir_ / "traj_stride.xyz";
    cpu::minimizer::Minimize opt;
    opt.add(cpu::minimizer::LBFGS{5, 0.01, 10});
    opt.add_handler(cpu::ForceConverged{/*ftol=*/0.5});
    opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
    opt.on(minimizer::Phase::end_of_step)
        .add(
            cpu::minimizer::TrajectoryWriter{traj_path.string(), /*stride=*/2});
    opt.build();
    opt.run(frame, 4);

    // stride=2 over 4 steps writes frames at step 0,2 -> 2 frames.
    int n_atoms = 2;
    int n_lines = count_lines(traj_path.string());
    int n_frames = n_lines / (n_atoms + 2);
    EXPECT_EQ(n_frames, 2);
}

TEST_F(MinimizerRecordersTest, OptLogTSVHasOneRowPerStep) {
    auto frame = make_h2(1.2);
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    const auto log_path = tmpdir_ / "opt.tsv";
    cpu::minimizer::Minimize opt;
    opt.add(cpu::minimizer::LBFGS{5, 0.01, 10});
    opt.add_handler(cpu::ForceConverged{/*ftol=*/0.5});
    opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
    opt.on(minimizer::Phase::end_of_step)
        .add(cpu::minimizer::OptLog{log_path.string()});
    opt.build();
    opt.run(frame, 3);

    std::ifstream is(log_path);
    ASSERT_TRUE(is.good());
    std::string header;
    std::getline(is, header);
    EXPECT_NE(header.find("step"), std::string::npos);
    EXPECT_NE(header.find("energy_eV"), std::string::npos);
    int n_rows = 0;
    std::string line;
    while (std::getline(is, line))
        if (!line.empty())
            ++n_rows;
    EXPECT_EQ(n_rows, 3);
}

TEST_F(MinimizerRecordersTest, FinalGeometryContainsOneFrame) {
    auto frame = make_h2(1.2);
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    const auto final_path = tmpdir_ / "final.xyz";
    cpu::minimizer::Minimize opt;
    opt.add(cpu::minimizer::LBFGS{5, 0.01, 10});
    opt.add_handler(cpu::ForceConverged{/*ftol=*/0.5});
    opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
    opt.on(minimizer::Phase::end_of_step)
        .add(cpu::minimizer::FinalGeometry{final_path.string()});
    opt.build();
    opt.run(frame, 5);

    int n_atoms = 2;
    int n_lines = count_lines(final_path.string());
    // Exactly one frame = n_atoms + 2 lines.
    EXPECT_EQ(n_lines, n_atoms + 2);
}
