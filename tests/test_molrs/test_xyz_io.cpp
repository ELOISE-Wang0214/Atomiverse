/// @file tests/test_molrs/test_xyz_io.cpp
/// @brief Drift fence for molrs XYZ / single-frame IO doors ATV uses.

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/io/extxyz.h"
#include "atomiverse/cpu/io/frame_writer.h"
#include "atomiverse/molrs.hpp"

namespace {

std::string tmp_path(const char* name) {
    return (std::filesystem::temp_directory_path() / name).string();
}

void remove_path(const std::string& path) {
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(MolrsXyzIo, WriteFrameXyzThenXyzReadFirstFrame) {
    const std::string path = tmp_path("atv_molrs_xyz_io.xyz");
    remove_path(path);
    const std::vector<std::int32_t> type = {1, 1};
    const std::vector<double> x = {0.0, 0.74};
    const std::vector<double> y = {0.0, 0.0};
    const std::vector<double> z = {0.0, 0.0};
    const std::vector<double> box;
    molrs::write_frame_xyz(
        {path.data(), path.size()},
        ::rust::Slice<const std::int32_t>(type.data(), type.size()),
        ::rust::Slice<const double>(x.data(), x.size()),
        ::rust::Slice<const double>(y.data(), y.size()),
        ::rust::Slice<const double>(z.data(), z.size()),
        ::rust::Slice<const double>(box.data(), box.size()),
        /*append=*/false);

    auto fr = molrs::xyz_read_first_frame({path.data(), path.size()});
    auto xr = molrs::frame_column_f64(*fr, "atoms", "x");
    ASSERT_EQ(xr.size(), 2u);
    EXPECT_NEAR(xr[1], 0.74, 1e-6);
    remove_path(path);
}

TEST(MolrsXyzIo, WriteFrameXyzTypedAcceptsEmptyMeta) {
    const std::string path = tmp_path("atv_molrs_xyz_typed.xyz");
    remove_path(path);
    const std::vector<std::int32_t> type = {8};
    const std::vector<double> x = {0.0};
    const std::vector<double> y = {0.0};
    const std::vector<double> z = {0.0};
    const std::vector<double> box = {5, 0, 0, 0, 5, 0, 0, 0, 5};
    rust::Vec<molrs::MetaEntry> meta;
    molrs::write_frame_xyz_typed(
        {path.data(), path.size()},
        ::rust::Slice<const std::int32_t>(type.data(), type.size()),
        ::rust::Slice<const double>(x.data(), x.size()),
        ::rust::Slice<const double>(y.data(), y.size()),
        ::rust::Slice<const double>(z.data(), z.size()),
        ::rust::Slice<const double>(box.data(), box.size()), std::move(meta),
        /*append=*/false);
    EXPECT_TRUE(std::filesystem::exists(path));
    remove_path(path);
}

TEST(MolrsXyzIo, CpuWriteFrameXyzWrapper) {
    const std::string path = tmp_path("atv_molrs_cpu_write_xyz.xyz");
    remove_path(path);
    const std::array<double, 9> box = {10, 0, 0, 0, 10, 0, 0, 0, 10};
    ATV::cpu::write_frame_xyz(path, {1}, {0.0}, {0.0}, {0.0}, box,
                              /*append=*/false);
    EXPECT_TRUE(std::filesystem::exists(path));
    remove_path(path);
}

TEST(MolrsXyzIo, ExtxyzMissingFileThrows) {
    EXPECT_THROW(ATV::cpu::load_extxyz("/no/such/atv_molrs_missing.xyz"),
                 std::runtime_error);
}

} // namespace
