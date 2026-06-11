/// @file tests/cpu/recorder/test_density_cube.cpp
/// @brief Tests for cpu::recorder::DensityCubeWriter and OrbitalCubeWriter.
///
/// Core acceptance: the grid-integrated density ≈ n_electrons within a
/// grid-resolution-dependent tolerance.

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/recorder/density_cube_writer.h"
#include "atomiverse/cpu/recorder/orbital_cube_writer.h"
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

struct ParsedCube {
    int n_atoms;
    double origin[3];
    int dims[3];
    double spacing[3]; // diagonal
    std::vector<double> values;
    double voxel() const { return spacing[0] * spacing[1] * spacing[2]; }
};

// Parse a Gaussian cube file that our writer produced (diagonal axes).
ParsedCube parse_cube(const std::string& path) {
    std::ifstream is(path);
    if (!is)
        throw std::runtime_error("parse_cube: cannot open " + path);
    std::string line;
    std::getline(is, line); // title1
    std::getline(is, line); // title2
    ParsedCube c;
    std::getline(is, line);
    {
        std::istringstream ss(line);
        ss >> c.n_atoms >> c.origin[0] >> c.origin[1] >> c.origin[2];
    }
    for (int axis = 0; axis < 3; ++axis) {
        std::getline(is, line);
        std::istringstream ss(line);
        double ax, ay, az;
        ss >> c.dims[axis] >> ax >> ay >> az;
        c.spacing[axis] = (axis == 0) ? ax : (axis == 1) ? ay : az;
    }
    for (int a = 0; a < c.n_atoms; ++a)
        std::getline(is, line);
    long total = static_cast<long>(c.dims[0]) * c.dims[1] * c.dims[2];
    c.values.reserve(total);
    double v;
    while (is >> v)
        c.values.push_back(v);
    if (static_cast<long>(c.values.size()) != total)
        throw std::runtime_error("parse_cube: incomplete value grid");
    return c;
}

} // namespace

TEST(DensityCubeTest, H2DensityIntegratesToTwoElectrons) {
    auto frame = make_h2(0.74);
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::scf::SCFEvaluator scf{config};
    auto handle = scf.result_handle();

    auto tmp = std::filesystem::temp_directory_path() / "atv_rho_h2.cube";
    std::error_code ec;
    std::filesystem::remove(tmp, ec);

    cpu::minimizer::Minimize opt;
    opt.add(cpu::minimizer::LBFGS{5, 0.01, 0.5, 10});
    opt.on(minimizer::Phase::force).add(std::move(scf));
    opt.on(minimizer::Phase::end_of_step)
        .add(cpu::recorder::DensityCubeWriter{
            handle, tmp.string(), /*padding=*/6.0, /*spacing=*/0.2});
    opt.build();
    opt.run(frame, 1);

    ASSERT_TRUE(std::filesystem::exists(tmp));
    auto cube = parse_cube(tmp.string());
    EXPECT_EQ(cube.n_atoms, 2);
    double sum = 0;
    for (double v : cube.values)
        sum += v;
    double integrated = sum * cube.voxel();
    // 0.2-Bohr grid over a 6-Bohr-padded H2 box -- about 1% accurate.
    EXPECT_NEAR(integrated, 2.0, 0.1)
        << "grid-integrated density should approximate 2 electrons";

    std::filesystem::remove(tmp, ec);
}

TEST(OrbitalCubeTest, WritesOrbitalFileForHomoH2) {
    auto frame = make_h2(0.74);
    auto integrals =
        cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::scf::SCFEvaluator scf{config};
    auto handle = scf.result_handle();

    auto tmp = std::filesystem::temp_directory_path() / "atv_mo_h2.cube";
    std::error_code ec;
    std::filesystem::remove(tmp, ec);

    int homo_idx = config.n_alpha - 1;
    cpu::minimizer::Minimize opt;
    opt.add(cpu::minimizer::LBFGS{5, 0.01, 0.5, 10});
    opt.on(minimizer::Phase::force).add(std::move(scf));
    opt.on(minimizer::Phase::end_of_step)
        .add(cpu::recorder::OrbitalCubeWriter{
            handle, tmp.string(), homo_idx, cpu::recorder::OrbitalSpin::Alpha,
            /*padding=*/5.0, /*spacing=*/0.3});
    opt.build();
    opt.run(frame, 1);

    ASSERT_TRUE(std::filesystem::exists(tmp));
    auto cube = parse_cube(tmp.string());
    EXPECT_EQ(cube.n_atoms, 2);
    // Orbital should be normalized: integral of |psi|^2 ≈ 1.
    double norm = 0;
    for (double v : cube.values)
        norm += v * v;
    norm *= cube.voxel();
    EXPECT_NEAR(norm, 1.0, 0.1);
    std::filesystem::remove(tmp, ec);
}

TEST(DensityCubeTest, NullHandleThrows) {
    cpu::scf::SCFResultHandle empty{};
    EXPECT_THROW(cpu::recorder::DensityCubeWriter(empty, "/tmp/x.cube"),
                 std::invalid_argument);
    EXPECT_THROW(
        cpu::recorder::OrbitalCubeWriter(empty, "/tmp/x.cube", 0,
                                         cpu::recorder::OrbitalSpin::Alpha),
        std::invalid_argument);
}
