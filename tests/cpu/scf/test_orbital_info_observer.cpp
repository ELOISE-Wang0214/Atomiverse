#include "atomiverse/cpu/scf/methods.h"
/// @file tests/cpu/scf/test_orbital_info_observer.cpp
/// @brief Tests for cpu::scf::OrbitalInfoObserver.

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/scf/internal/scf_types.h"
#include "atomiverse/cpu/scf/orbital_info_observer.h"

using namespace ATV;

namespace {

cpu::Frame make_h2(double bond_length) {
    cpu::Frame f(2);
    f.span(ATV::kAtom_x)[0] = 0.0;
    f.span(ATV::kAtom_x)[1] = bond_length;
    f.span(ATV::kAtom_y)[0] = 0.0;
    f.span(ATV::kAtom_y)[1] = 0.0;
    f.span(ATV::kAtom_z)[0] = 0.0;
    f.span(ATV::kAtom_z)[1] = 0.0;
    f.span(ATV::kAtom_mass)[0] = 1.008;
    f.span(ATV::kAtom_mass)[1] = 1.008;
    f.span(ATV::kAtom_type)[0] = 1;
    f.span(ATV::kAtom_type)[1] = 1;
    f.simbox = ATV::cpu::Box::cubic(100.0);
    return f;
}

} // namespace

TEST(OrbitalInfoObserverTest, WritesHeaderAndRowsForH2) {
    auto frame = make_h2(0.74);
    auto tmp = std::filesystem::temp_directory_path() / "atv_obs_h2.tsv";
    std::error_code ec;
    std::filesystem::remove(tmp, ec);

    cpu::scf::RHF scf{};
    auto handle = scf.result_handle();

    cpu::minimizer::Minimize opt;
    opt.add(cpu::minimizer::LBFGS{5, 0.01, 10});
    opt.add_criterion(cpu::ForceConverged{/*ftol=*/0.5});
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
        if (line.empty())
            continue;
        ++rows;
        int step;
        double homo, lumo, gap;
        std::sscanf(line.c_str(), "%d\t%lf\t%lf\t%lf", &step, &homo, &lumo,
                    &gap);
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
    cpu::scf::MeanFieldOutputs empty{};
    EXPECT_THROW(cpu::scf::OrbitalInfoObserver(empty, ""),
                 std::invalid_argument);
}
