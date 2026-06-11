/// @file tests/mpi/test_gather.cpp
/// @brief Phase 4 -- `gather_to_root` correctness across rank counts.
///
/// Verifies that:
///   - Single-rank gather copies the local owner range verbatim.
///   - Multi-rank gather collects every owner exactly once and sorts the
///     result by global tag (so I/O is decomposition-independent).

#include <gtest/gtest.h>

#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/driver/backend.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/comm/domain.h"
#include "atomiverse/cpu/comm/gather.h"

using namespace ATV;

namespace {

// Lay 27 atoms on a 3x3x3 grid in a cubic box. Each atom gets a deterministic
// position so we can assert the gathered field values match exactly.
cpu::Frame make_grid_frame(double box_len, int per_axis) {
    const int n = per_axis * per_axis * per_axis;
    cpu::Frame frame(n);
    const double dx = box_len / per_axis;
    int idx = 0;
    for (int ix = 0; ix < per_axis; ++ix) {
        for (int iy = 0; iy < per_axis; ++iy) {
            for (int iz = 0; iz < per_axis; ++iz) {
                frame.atoms.x[idx] = (ix + 0.5) * dx;
                frame.atoms.y[idx] = (iy + 0.5) * dx;
                frame.atoms.z[idx] = (iz + 0.5) * dx;
                frame.atoms.vx[idx] = idx * 0.01;
                frame.atoms.vy[idx] = idx * 0.02;
                frame.atoms.vz[idx] = idx * 0.03;
                frame.atoms.mass[idx] = 1.0 + 0.001 * idx;
                frame.atoms.type[idx] = idx % 3;
                ++idx;
            }
        }
    }
    frame.box = Mat3::cubic(box_len);
    return frame;
}

} // namespace

TEST(Gather, SingleRankCopiesOwners) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size != 1)
        GTEST_SKIP() << "single-rank guard";

    auto frame = make_grid_frame(30.0, 3); // 27 atoms
    cpu::Backend backend;
    backend.init();
    backend.set_frame(frame);

    auto* dd = backend.context().dd;
    ASSERT_NE(dd, nullptr);

    auto g = cpu::md::mpi::gather_to_root(frame, *dd);
    EXPECT_EQ(g.n_total, 27);
    ASSERT_EQ(static_cast<int>(g.x.size()), 27);
    for (int i = 0; i < 27; ++i) {
        EXPECT_DOUBLE_EQ(g.x[i], frame.atoms.x[i]);
        EXPECT_DOUBLE_EQ(g.vx[i], frame.atoms.vx[i]);
        EXPECT_EQ(g.tag[i], frame.atoms.tag[i]);
        EXPECT_EQ(g.type[i], frame.atoms.type[i]);
    }

    backend.shutdown();
}

TEST(Gather, MultiRankSortsByTag) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size <= 1)
        GTEST_SKIP() << "multi-rank test";

    auto frame = make_grid_frame(30.0, 3); // 27 atoms
    cpu::Backend backend;
    backend.init();
    backend.set_frame(frame);

    auto* dd = backend.context().dd;
    ASSERT_NE(dd, nullptr);

    int my_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    auto g = cpu::md::mpi::gather_to_root(frame, *dd);
    EXPECT_EQ(g.n_total, 27);

    if (my_rank == 0) {
        ASSERT_EQ(static_cast<int>(g.tag.size()), 27);
        // Tags must cover [0, 27) once each, in ascending order.
        for (int i = 0; i < 27; ++i) {
            EXPECT_EQ(g.tag[i], static_cast<std::int64_t>(i))
                << "tags not sorted ascending at index " << i;
        }
        // Type field index aligns with tag (since type[i] = i % 3 was the
        // original assignment, and tags are global atom IDs == initial idx).
        for (int i = 0; i < 27; ++i) {
            EXPECT_EQ(g.type[i], i % 3);
        }
    } else {
        // Non-root receives empty buffers but a valid n_total.
        EXPECT_TRUE(g.x.empty());
    }

    backend.shutdown();
}

int main(int argc, char** argv) {
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    ::testing::InitGoogleTest(&argc, argv);
    const int rc = RUN_ALL_TESTS();
    MPI_Finalize();
    return rc;
}
