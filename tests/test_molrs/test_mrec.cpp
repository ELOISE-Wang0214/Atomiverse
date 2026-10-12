/// @file tests/test_molrs/test_mrec.cpp
/// @brief Drift fence for molrs MrecWriterRef / read_mrec_trajectory_frame and
///        write_mrec_frame / read_mrec_frame.

#include <cstdint>
#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/molrs.hpp"
#include "rust/cxx.h"

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
    molrs::frame_set_column_u64(
        *fr, "atoms", "atomic_number",
        ::rust::Slice<const std::uint64_t>(znum.data(), znum.size()));
    molrs::frame_set_column_f64(
        *fr, "atoms", "mass",
        ::rust::Slice<const double>(mass.data(), mass.size()));
    const std::vector<double> h = {20, 0, 0, 0, 20, 0, 0, 0, 20};
    molrs::frame_set_box_h(*fr,
                           ::rust::Slice<const double>(h.data(), h.size()));
    return fr;
}

TEST(MolrsMrec, MrecWriterCreateAppendFlushClose) {
    const std::string path = tmp_mrec("atv_molrs_traj.mrec");
    auto schema = one_atom();
    auto writer = molrs::mrec_writer_create({path.data(), path.size()}, *schema,
                                            /*flush_every=*/1,
                                            /*durable=*/false);
    auto frame = one_atom();
    molrs::mrec_writer_append(*writer, *frame, /*step=*/0, /*time=*/0.0,
                              /*has_time=*/true);
    molrs::mrec_writer_flush(*writer);
    EXPECT_GE(molrs::mrec_writer_committed(*writer), 1u);
    molrs::mrec_writer_close(std::move(writer));

    auto back = molrs::read_mrec_trajectory_frame({path.data(), path.size()},
                                                  /*index=*/0);
    auto z = molrs::frame_column_u64(*back, "atoms", "atomic_number");
    ASSERT_EQ(z.size(), 1u);
    EXPECT_EQ(z[0], 18u);
    EXPECT_THROW(molrs::read_mrec_trajectory_frame({path.data(), path.size()},
                                                   /*index=*/1),
                 rust::Error);
    std::error_code ec;
    std::filesystem::remove_all(path, ec);
}

TEST(MolrsMrec, WriteMrecFrameThenReadMrecFrame) {
    const std::string path = tmp_mrec("atv_molrs_single.mrec");
    const std::vector<std::int32_t> atomic_number = {1};
    const std::vector<double> x = {0.1};
    const std::vector<double> y = {0.2};
    const std::vector<double> z = {0.3};
    const std::vector<double> box = {5, 0, 0, 0, 5, 0, 0, 0, 5};
    rust::Vec<rust::String> names;
    molrs::write_mrec_frame({path.data(), path.size()},
                            ::rust::Slice<const std::int32_t>(
                                atomic_number.data(), atomic_number.size()),
                            ::rust::Slice<const double>(x.data(), x.size()),
                            ::rust::Slice<const double>(y.data(), y.size()),
                            ::rust::Slice<const double>(z.data(), z.size()),
                            ::rust::Slice<const double>(box.data(), box.size()),
                            std::move(names),
                            ::rust::Slice<const double>(nullptr, 0));
    auto fr = molrs::read_mrec_frame({path.data(), path.size()});
    auto xr = molrs::frame_column_f64(*fr, "atoms", "x");
    ASSERT_EQ(xr.size(), 1u);
    EXPECT_NEAR(xr[0], 0.1, 1e-9);
    std::error_code ec;
    std::filesystem::remove_all(path, ec);
}

} // namespace
