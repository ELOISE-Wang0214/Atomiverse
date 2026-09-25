/// @file tests/test_molrs/test_mrec.cpp
/// @brief Drift fence for molrs TrajectoryWriter + write_frame / read_first_frame.

#include <cstdint>
#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/molrs.hpp"

#if !defined(ATV_HAS_MOLRS_TRAJECTORY_WRITER)
#error "test_mrec.cpp requires ATV_HAS_MOLRS_TRAJECTORY_WRITER"
#endif

namespace {

std::string tmp_mrec(const char* name) {
    auto p = std::filesystem::temp_directory_path() / name;
    std::error_code ec;
    std::filesystem::remove_all(p, ec);
    return p.string();
}

::rust::Box<molrs::FrameRef> one_atom() {
    auto fr = molrs::frame_new();
    const std::vector<double> x = {0.0};
    const std::vector<double> y = {0.0};
    const std::vector<double> z = {0.0};
    const std::vector<std::uint64_t> znum = {18};
    const std::vector<double> mass = {39.948};
    molrs::frame_set_column_f64(
        *fr, "atoms", "x", ::rust::Slice<const double>(x.data(), x.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "y", ::rust::Slice<const double>(y.data(), y.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "z", ::rust::Slice<const double>(z.data(), z.size()));
    molrs::frame_set_column_u32(
        *fr, "atoms", "atomic_number",
        ::rust::Slice<const std::uint64_t>(znum.data(), znum.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "mass",
        ::rust::Slice<const double>(mass.data(), mass.size()));
    const std::vector<double> h = {20, 0, 0, 0, 20, 0, 0, 0, 20};
    molrs::frame_set_box(*fr, ::rust::Slice<const double>(h.data(), h.size()));
    return fr;
}

TEST(MolrsMrec, TrajectoryWriterCreateAppendFlushClose) {
    const std::string path = tmp_mrec("atv_molrs_traj.mrec");
    auto schema = one_atom();
    auto writer = molrs::trajectory_writer_create({path.data(), path.size()},
                                                  *schema, /*flush_every=*/1,
                                                  /*durable=*/false);
    auto frame = one_atom();
    molrs::trajectory_writer_append(*writer, *frame, /*step=*/0, /*time=*/0.0,
                                    /*has_time=*/true);
    molrs::trajectory_writer_flush(*writer);
    EXPECT_GE(molrs::trajectory_writer_committed(*writer), 1u);
    molrs::trajectory_writer_close(std::move(writer));

    auto back = molrs::read_first_frame({path.data(), path.size()});
    auto z = molrs::frame_column_u32(*back, "atoms", "atomic_number");
    ASSERT_EQ(z.size(), 1u);
    EXPECT_EQ(z[0], 18u);
    std::error_code ec;
    std::filesystem::remove_all(path, ec);
}

TEST(MolrsMrec, WriteFrameSingleRecordThenRead) {
    const std::string path = tmp_mrec("atv_molrs_single.mrec");
    const std::vector<std::int32_t> type = {1};
    const std::vector<double> x = {0.1};
    const std::vector<double> y = {0.2};
    const std::vector<double> z = {0.3};
    const std::vector<double> box = {5, 0, 0, 0, 5, 0, 0, 0, 5};
    rust::Vec<rust::String> names;
    molrs::write_frame(
        {path.data(), path.size()},
        ::rust::Slice<const std::int32_t>(type.data(), type.size()),
        ::rust::Slice<const double>(x.data(), x.size()),
        ::rust::Slice<const double>(y.data(), y.size()),
        ::rust::Slice<const double>(z.data(), z.size()),
        ::rust::Slice<const double>(box.data(), box.size()), std::move(names),
        ::rust::Slice<const double>(nullptr, 0));
    auto fr = molrs::read_first_frame({path.data(), path.size()});
    auto xr = molrs::frame_column_f64(*fr, "atoms", "x");
    ASSERT_EQ(xr.size(), 1u);
    EXPECT_NEAR(xr[0], 0.1, 1e-9);
    std::error_code ec;
    std::filesystem::remove_all(path, ec);
}

} // namespace
