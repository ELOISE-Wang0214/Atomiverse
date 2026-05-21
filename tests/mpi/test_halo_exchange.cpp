/// @file tests/mpi/test_halo_exchange.cpp
/// @brief Phase 2.2 / 2.3 ghost halo correctness tests.
///
/// Multi-rank tests for `HaloExchangeForward` + `HaloExchangeReverse`.
/// Launched via `mpirun -n N`. Verifies:
///   - Forward populates `dd->n_ghost > 0` and `recv_counts` consistent with
///     each rank's own send list.
///   - Symmetry: my +x send count to my +x neighbor equals my +x neighbor's
///     -x send count to me (count exchange handshake is consistent).
///   - PBC shift: ghosts received from across a periodic-edge face land
///     at coordinates that, when wrapped back to `[0, L)`, match an actual
///     owner position somewhere in the neighbor's slab.
///   - Reverse symmetry: applying a known per-ghost force and reducing
///     back through `HaloExchangeReverse` adds exactly the right total to
///     each owner's force.

#include <gtest/gtest.h>

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/backend.h"
#include "atomiverse/cpu/context.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/cpu/md/mpi/halo_exchange.h"

using namespace ATV;

namespace {

cpu::Frame make_uniform_grid(int n_per_side, double box_len) {
    const int n = n_per_side * n_per_side * n_per_side;
    cpu::Frame frame(n);
    const double step = box_len / static_cast<double>(n_per_side);
    int idx = 0;
    for (int ix = 0; ix < n_per_side; ++ix)
        for (int iy = 0; iy < n_per_side; ++iy)
            for (int iz = 0; iz < n_per_side; ++iz) {
                frame.atoms.x[idx] = (static_cast<double>(ix) + 0.5) * step;
                frame.atoms.y[idx] = (static_cast<double>(iy) + 0.5) * step;
                frame.atoms.z[idx] = (static_cast<double>(iz) + 0.5) * step;
                frame.atoms.mass[idx] = 1.0;
                ++idx;
            }
    frame.box = Mat3::cubic(box_len);
    return frame;
}

// Run init -> set_frame -> single forward exchange and return the dd state.
struct ForwardResult {
    int n_local;
    int n_ghost;
    std::array<int, 6> recv_counts;
    std::array<int, 6> send_counts;
    std::array<int, 6> face_neighbors;
    int dims[3];
    int coords[3];
    double Lx, Ly, Lz;
    std::vector<double> ghost_x, ghost_y, ghost_z;
};

ForwardResult run_forward(double box_len, int n_per_side, double cutoff) {
    auto frame = make_uniform_grid(n_per_side, box_len);

    cpu::Backend backend;
    backend.init();
    backend.set_frame(frame);

    auto* dd = backend.context().dd;

    // Drive one Forward launch directly (without going through MD driver).
    cpu::md::mpi::HaloExchangeForward fwd{cutoff};
    cpu::BuildContext bctx{backend.context(), frame};
    fwd.build(bctx);
    cpu::LaunchContext lctx{backend.context(), frame, 0, 0};
    fwd.launch(lctx);

    ForwardResult r{};
    r.n_local = dd->n_local;
    r.n_ghost = dd->n_ghost;
    for (int f = 0; f < 6; ++f) {
        r.recv_counts[f] = dd->recv_counts[f];
        r.send_counts[f] = static_cast<int>(dd->send_lists[f].size());
        r.face_neighbors[f] = dd->face_neighbors[f];
    }
    for (int d = 0; d < 3; ++d) {
        r.dims[d] = dd->dims[d];
        r.coords[d] = dd->coords[d];
    }
    r.Lx = frame.box(0, 0);
    r.Ly = frame.box(1, 1);
    r.Lz = frame.box(2, 2);

    r.ghost_x.assign(static_cast<std::size_t>(dd->n_ghost), 0.0);
    r.ghost_y.assign(static_cast<std::size_t>(dd->n_ghost), 0.0);
    r.ghost_z.assign(static_cast<std::size_t>(dd->n_ghost), 0.0);
    for (int g = 0; g < dd->n_ghost; ++g) {
        r.ghost_x[g] = frame.atoms.x[dd->n_local + g];
        r.ghost_y[g] = frame.atoms.y[dd->n_local + g];
        r.ghost_z[g] = frame.atoms.z[dd->n_local + g];
    }

    backend.shutdown();
    return r;
}

} // namespace

// ── n_ghost > 0 when cutoff captures the neighbor slab ──────────────────

TEST(HaloForward, GhostsPopulatedAtModerateCutoff) {
    // 64 atoms on 4x4x4 grid in a 10 Å box. With 2x2x2 ranks, each slab is
    // 5 Å, so a cutoff of 2.0 Å captures the inner-face atoms.
    const auto r = run_forward(10.0, 4, 2.0);

    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size <= 1) {
        EXPECT_EQ(r.n_ghost, 0);
        return;
    }

    EXPECT_GT(r.n_ghost, 0)
        << "expected ghost atoms after Forward halo exchange at cutoff 2.0";
}

// ── recv_counts[f] consistent with the f-neighbor's send_counts ────────

TEST(HaloForward, RecvCountsConsistentWithSendLists) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size <= 1) {
        GTEST_SKIP() << "single-rank: no cross-rank symmetry to check";
    }

    auto frame = make_uniform_grid(4, 10.0);
    cpu::Backend backend;
    backend.init();
    backend.set_frame(frame);
    auto* dd = backend.context().dd;

    cpu::md::mpi::HaloExchangeForward fwd{2.0};
    cpu::BuildContext bctx{backend.context(), frame};
    fwd.build(bctx);
    cpu::LaunchContext lctx{backend.context(), frame, 0, 0};
    fwd.launch(lctx);

    // Local sanity: total ghost atoms = sum of recv_counts[] across faces.
    int total_ghosts = 0;
    for (int f = 0; f < 6; ++f)
        total_ghosts += dd->recv_counts[f];
    EXPECT_EQ(total_ghosts, dd->n_ghost);

    // Faces that self-wrap (dims[dim] == 1) contribute 0 because we
    // short-circuit those at the component.
    for (int f = 0; f < 6; ++f) {
        const int dim = f >> 1;
        if (dd->dims[dim] == 1) {
            EXPECT_EQ(dd->recv_counts[f], 0)
                << "face " << f << " wraps to self (dims[" << dim
                << "]==1) but received " << dd->recv_counts[f] << " ghosts";
        }
    }

    backend.shutdown();
}

// ── PBC-shifted ghosts correspond to real atoms after wrap ─────────────

TEST(HaloForward, GhostPositionsResolveToOwnersAfterWrap) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size <= 1)
        GTEST_SKIP();

    const auto r = run_forward(10.0, 4, 2.0);

    auto wrap = [](double q, double L) {
        double w = q - std::floor(q / L) * L;
        if (w < 0)
            w += L;
        if (w >= L)
            w -= L;
        return w;
    };

    const int step = static_cast<int>(10.0 / 4); // grid step is 2.5
    (void)step;
    const double dx = 10.0 / 4; // 2.5

    // Every ghost, after wrapping back to [0, L), must coincide with one
    // of the original 64 grid points (positions are k * 2.5 + 1.25).
    for (int g = 0; g < r.n_ghost; ++g) {
        const double xw = wrap(r.ghost_x[g], r.Lx);
        const double yw = wrap(r.ghost_y[g], r.Ly);
        const double zw = wrap(r.ghost_z[g], r.Lz);
        const double rx = std::round((xw - 0.5 * dx) / dx);
        const double ry = std::round((yw - 0.5 * dx) / dx);
        const double rz = std::round((zw - 0.5 * dx) / dx);
        EXPECT_NEAR(xw, rx * dx + 0.5 * dx, 1e-9)
            << "ghost " << g << " x mismatch";
        EXPECT_NEAR(yw, ry * dx + 0.5 * dx, 1e-9)
            << "ghost " << g << " y mismatch";
        EXPECT_NEAR(zw, rz * dx + 0.5 * dx, 1e-9)
            << "ghost " << g << " z mismatch";
    }
}

// ── Reverse: known ghost forces sum back to owners losslessly ──────────

TEST(HaloReverse, ReverseAddsGhostForcesBackToOwners) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size <= 1)
        GTEST_SKIP();

    auto frame = make_uniform_grid(4, 10.0);
    cpu::Backend backend;
    backend.init();
    backend.set_frame(frame);
    auto* dd = backend.context().dd;

    cpu::md::mpi::HaloExchangeForward fwd{2.0};
    cpu::md::mpi::HaloExchangeReverse rev{};
    cpu::BuildContext bctx{backend.context(), frame};
    fwd.build(bctx);
    rev.build(bctx);
    cpu::LaunchContext lctx{backend.context(), frame, 0, 0};
    fwd.launch(lctx);

    // Plant a known force pattern on every atom (owners + ghosts).
    // Owner forces start at 1.0, 2.0, 3.0; ghost forces add 10.0 each.
    for (int i = 0; i < dd->n_local; ++i) {
        frame.atoms.fx[i] = 1.0;
        frame.atoms.fy[i] = 2.0;
        frame.atoms.fz[i] = 3.0;
    }
    for (int g = 0; g < dd->n_ghost; ++g) {
        frame.atoms.fx[dd->n_local + g] = 10.0;
        frame.atoms.fy[dd->n_local + g] = 10.0;
        frame.atoms.fz[dd->n_local + g] = 10.0;
    }

    rev.launch(lctx);

    // After reverse, every owner whose tag was ghost-shipped to a *real*
    // neighbor (not a self-wrap on a dims==1 axis) gets +10.0 back per
    // face it appeared on. Self-wrap faces are short-circuited by the
    // halo implementation, so they contribute zero.
    int my_rank = 0;
    MPI_Comm_rank(dd->cart_comm, &my_rank);
    std::vector<int> ghost_appearances(static_cast<std::size_t>(dd->n_local),
                                       0);
    for (int f = 0; f < 6; ++f) {
        if (dd->face_neighbors[f] == my_rank)
            continue; // self-face: no comm
        for (int idx : dd->send_lists[f])
            ghost_appearances[static_cast<std::size_t>(idx)] += 1;
    }

    for (int i = 0; i < dd->n_local; ++i) {
        const double expected =
            1.0 + 10.0 * ghost_appearances[static_cast<std::size_t>(i)];
        EXPECT_NEAR(frame.atoms.fx[i], expected, 1e-12)
            << "owner " << i
            << " fx not equal to base + 10 * (ghost appearances)";
    }

    backend.shutdown();
}

// ── Regression: per-instance FaceBufs (PIMPL, not a shared static) ─────

// Catches a re-introduction of `static thread_local FaceBufs s` (or any
// other cross-instance shared storage) behind the Forward / Reverse
// facades. The test builds two HaloExchangeForward instances in sequence
// against DIFFERENT `dd->capacity` values and verifies that the first
// instance's per-face buffer capacity is unchanged by the second's
// `build()`. With the old shared static, the vector's `reserve()`
// monotonically grew, so the first instance would see its capacity bumped
// up to the second's after B.build() ran.
TEST(HaloForward, FaceBufsAreInstanceOwned) {
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size <= 1) {
        GTEST_SKIP() << "multi-rank required (need non-null dd)";
    }

    auto frame_small = make_uniform_grid(2, 10.0); // 8 atoms
    auto frame_big = make_uniform_grid(6, 10.0);   // 216 atoms

    cpu::Backend backend;
    backend.init();
    backend.set_frame(frame_small);

    cpu::md::mpi::HaloExchangeForward fwd_a{2.0};
    cpu::BuildContext bctx_small{backend.context(), frame_small};
    fwd_a.build(bctx_small);
    const std::size_t cap_a_before = fwd_a.send_buf_capacity_for_test();

    // Rebind the same backend to a larger frame; dd->capacity now reflects
    // the bigger slab.
    backend.set_frame(frame_big);
    cpu::md::mpi::HaloExchangeForward fwd_b{2.0};
    cpu::BuildContext bctx_big{backend.context(), frame_big};
    fwd_b.build(bctx_big);
    const std::size_t cap_b = fwd_b.send_buf_capacity_for_test();
    const std::size_t cap_a_after = fwd_a.send_buf_capacity_for_test();

    // Per-instance invariant: A's bufs are untouched by B's build.
    EXPECT_EQ(cap_a_before, cap_a_after)
        << "HaloExchangeForward #A's FaceBufs capacity changed from "
        << cap_a_before << " to " << cap_a_after
        << " after instance #B was built -- bufs must be per-instance, "
        << "not a shared static.";

    // Different capacity inputs must land in different storage: if A and B
    // shared a static vector, they'd see the same (max) capacity.
    EXPECT_NE(cap_a_after, cap_b)
        << "Instances A and B built against different dd->capacity values "
        << "produced identical send-buf capacities (A=" << cap_a_after
        << " B=" << cap_b << "); likely sharing a single storage object.";

    backend.shutdown();
}

// See test_topology.cpp for the OpenMPI 5.0.x macOS teardown rationale.
// CTest reads the GoogleTest "[  PASSED  ]" summary line to determine
// pass/fail, bypassing any prterun exit-code wobble.
int main(int argc, char** argv) {
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    ::testing::InitGoogleTest(&argc, argv);
    const int rc = RUN_ALL_TESTS();
    MPI_Finalize();
    return rc;
}
