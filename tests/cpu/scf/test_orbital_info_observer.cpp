/// @file tests/cpu/scf/test_orbital_info_observer.cpp
/// @brief Tests for cpu::scf::OrbitalInfoObserver.

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/scf/orbital_info_observer.h"
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

} // namespace

TEST(OrbitalInfoObserverTest, WritesHeaderAndRowsForH2) {
    auto frame = make_h2(0.74);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    auto tmp = std::filesystem::temp_directory_path() / "atv_obs_h2.tsv";
    std::error_code ec; std::filesystem::remove(tmp, ec);

    cpu::scf::SCFEvaluator scf{config};
    auto handle = scf.result_handle();

    cpu::minimizer::LBFGS opt;
    opt.add(cpu::minimizer::LBFGSStep{5, 0.01, 0.5, 10});
    opt.on(minimizer::Phase::force).add(std::move(scf));
    opt.on(minimizer::Phase::end_of_step)
        .add(cpu::scf::OrbitalInfoObserver{handle, tmp.string()});
    opt.build();
    opt.run(frame, 2);

    std::ifstream is(tmp);
    ASSERT_TRUE(is.good());
    std::string header;
    std::getline(is, header);
    EXPECT_NE(header.find("homo_eV"), std::string::npos);
    EXPECT_NE(header.find("gap_eV"), std::string::npos);

    int rows = 0;
    double last_homo = 0, last_lumo = 0;
    std::string line;
    while (std::getline(is, line)) {
        if (line.empty()) continue;
        ++rows;
        int step;
        double homo, lumo, gap;
        std::sscanf(line.c_str(), "%d\t%lf\t%lf\t%lf", &step, &homo, &lumo, &gap);
        last_homo = homo;
        last_lumo = lumo;
    }
    EXPECT_EQ(rows, 2);
    // H2 / STO-3G HOMO sits near -11 eV and LUMO above 0 eV (gap ~ 20-25 eV).
    EXPECT_LT(last_homo, last_lumo);
    EXPECT_LT(last_homo, 0.0);
    std::filesystem::remove(tmp, ec);
}

TEST(OrbitalInfoObserverTest, NullHandleThrows) {
    cpu::scf::SCFResultHandle empty{};
    EXPECT_THROW(
        cpu::scf::OrbitalInfoObserver(empty, ""),
        std::invalid_argument);
}
