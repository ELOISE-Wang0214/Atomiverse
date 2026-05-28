/// @file tests/cpu/md/test_topology_resolved_locals.cpp
/// @brief Unit tests for the per-step local-index resolver shared by
///        Bond / Angle / Dihedral bases (CPU). Verifies the
///        owner-filter + per-step compaction contract on (a) a
///        single-rank pass-through and (b) a multi-rank-simulated
///        context with controlled `dd->tag_to_local`.

#include <gtest/gtest.h>

#include <cstdint>
#include <memory>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/context.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/cpu/potentials/angle.h"
#include "atomiverse/cpu/potentials/bond.h"
#include "atomiverse/cpu/potentials/dihedral.h"

using namespace ATV;

namespace {

// Minimal concrete subclasses that expose the protected resolver methods
// for testing. compute() is a no-op — the test never calls launch().

class TestBond : public cpu::Bond {
  public:
    OwnedBondsView resolve(const cpu::LaunchContext& ctx, int n_tuples,
                           const int* ai, const int* aj) {
        return resolve_owned_bonds(ctx, n_tuples, ai, aj);
    }
    void compute(const cpu::LaunchContext&) override {}
    std::string name() const override { return "test_bond"; }
};

class TestAngle : public cpu::Angle {
  public:
    OwnedAnglesView resolve(const cpu::LaunchContext& ctx, int n_tuples,
                            const int* ai, const int* aj, const int* ak) {
        return resolve_owned_angles(ctx, n_tuples, ai, aj, ak);
    }
    void compute(const cpu::LaunchContext&) override {}
    std::string name() const override { return "test_angle"; }
};

class TestDihedral : public cpu::Dihedral {
  public:
    OwnedDihedralsView resolve(const cpu::LaunchContext& ctx, int n_tuples,
                               const int* ai, const int* aj, const int* ak,
                               const int* al) {
        return resolve_owned_dihedrals(ctx, n_tuples, ai, aj, ak, al);
    }
    void compute(const cpu::LaunchContext&) override {}
    std::string name() const override { return "test_dihedral"; }
};

cpu::Frame make_chain(int n_atoms) {
    cpu::Frame frame(n_atoms);
    for (int i = 0; i < n_atoms; ++i) {
        frame.atoms.x[i] = static_cast<double>(i);
        frame.atoms.mass[i] = 14.0;
        frame.atoms.tag[i] = static_cast<std::int64_t>(i);
    }
    frame.box = Mat3::cubic(100.0);
    return frame;
}

} // namespace

// ── single-rank: identity mapping ────────────────────────────────────

TEST(BondResolvedLocals, SingleRankIsIdentity) {
    auto frame = make_chain(6);
    // 5 bonds along the chain: (0-1), (1-2), (2-3), (3-4), (4-5).
    frame.bonds.atomi = xt::xtensor<int, 1>{0, 1, 2, 3, 4};
    frame.bonds.atomj = xt::xtensor<int, 1>{1, 2, 3, 4, 5};
    frame.bonds.type = xt::zeros<int>({5});

    cpu::Context ctx;
    cpu::LaunchContext lctx{ctx, frame, 0, 0};
    TestBond b;
    auto v = b.resolve(lctx, 5, frame.bonds.atomi.data(),
                       frame.bonds.atomj.data());

    EXPECT_EQ(v.n_work, 5);
    for (int k = 0; k < 5; ++k) {
        EXPECT_EQ(v.tuple_idx[k], k);
        EXPECT_EQ(v.local_i[k], k);
        EXPECT_EQ(v.local_j[k], k + 1);
    }
}

TEST(AngleResolvedLocals, SingleRankIsIdentity) {
    auto frame = make_chain(6);
    // 4 angles along the chain: (0-1-2), (1-2-3), (2-3-4), (3-4-5).
    frame.angles.atomi = xt::xtensor<int, 1>{0, 1, 2, 3};
    frame.angles.atomj = xt::xtensor<int, 1>{1, 2, 3, 4};
    frame.angles.atomk = xt::xtensor<int, 1>{2, 3, 4, 5};
    frame.angles.type = xt::zeros<int>({4});

    cpu::Context ctx;
    cpu::LaunchContext lctx{ctx, frame, 0, 0};
    TestAngle a;
    auto v = a.resolve(lctx, 4, frame.angles.atomi.data(),
                       frame.angles.atomj.data(), frame.angles.atomk.data());

    EXPECT_EQ(v.n_work, 4);
    for (int k = 0; k < 4; ++k) {
        EXPECT_EQ(v.local_i[k], k);
        EXPECT_EQ(v.local_j[k], k + 1);
        EXPECT_EQ(v.local_k[k], k + 2);
    }
}

TEST(DihedralResolvedLocals, SingleRankIsIdentity) {
    auto frame = make_chain(6);
    // 3 dihedrals along the chain: (0-1-2-3), (1-2-3-4), (2-3-4-5).
    frame.dihedrals.atomi = xt::xtensor<int, 1>{0, 1, 2};
    frame.dihedrals.atomj = xt::xtensor<int, 1>{1, 2, 3};
    frame.dihedrals.atomk = xt::xtensor<int, 1>{2, 3, 4};
    frame.dihedrals.atoml = xt::xtensor<int, 1>{3, 4, 5};
    frame.dihedrals.type = xt::zeros<int>({3});

    cpu::Context ctx;
    cpu::LaunchContext lctx{ctx, frame, 0, 0};
    TestDihedral d;
    auto v = d.resolve(lctx, 3, frame.dihedrals.atomi.data(),
                       frame.dihedrals.atomj.data(),
                       frame.dihedrals.atomk.data(),
                       frame.dihedrals.atoml.data());

    EXPECT_EQ(v.n_work, 3);
    for (int k = 0; k < 3; ++k) {
        EXPECT_EQ(v.local_i[k], k);
        EXPECT_EQ(v.local_j[k], k + 1);
        EXPECT_EQ(v.local_k[k], k + 2);
        EXPECT_EQ(v.local_l[k], k + 3);
    }
}

// ── multi-rank simulated: partial tag_to_local ────────────────────────
//
// Configure a DD with mpi_size = 2, n_local = 4, and a tag_to_local that
// covers tags {0,1,2,3,4} (atoms 0..3 owned locally; atom 4 is a halo
// ghost slot at local row 4). Tags 5..9 are off-rank or out-of-halo.

TEST(BondResolvedLocals, MultiRankFilterAndCompact) {
    // Topology spans tags 0..9; this rank owns tags 0..3 (n_local=4)
    // and additionally has ghost slot for tag 4.
    auto frame = make_chain(10);
    // 5 bonds spanning various ownership/halo cases:
    //   0-1  owner min=0 (local) | partner 1 in halo  → kept
    //   2-3  owner min=2 (local) | partner 3 in halo  → kept
    //   3-4  owner min=3 (local) | partner 4 in halo  → kept (ghost)
    //   3-5  owner min=3 (local) | partner 5 NOT in halo → dropped (partner miss)
    //   4-5  owner min=4 NOT local (4 >= n_local)      → dropped (owner)
    //   6-7  owner min=6 NOT mapped                    → dropped (owner)
    frame.bonds.atomi = xt::xtensor<int, 1>{0, 2, 3, 3, 4, 6};
    frame.bonds.atomj = xt::xtensor<int, 1>{1, 3, 4, 5, 5, 7};
    frame.bonds.type = xt::zeros<int>({6});

    cpu::md::mpi::DomainDecomposition dd;
    dd.n_local = 4;
    dd.n_ghost = 1; // one halo ghost at row 4
    dd.capacity = 10;
    dd.neighbor_generation = 7;
    dd.tag_to_local = {{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}};

    cpu::Context ctx;
    ctx.mpi_size = 2;
    ctx.dd = &dd;
    cpu::LaunchContext lctx{ctx, frame, 0, 0};
    TestBond b;
    auto v = b.resolve(lctx, 6, frame.bonds.atomi.data(),
                       frame.bonds.atomj.data());

    ASSERT_EQ(v.n_work, 3) << "expected 0-1, 2-3, 3-4 to survive";
    // 0-1
    EXPECT_EQ(v.tuple_idx[0], 0);
    EXPECT_EQ(v.local_i[0], 0);
    EXPECT_EQ(v.local_j[0], 1);
    // 2-3
    EXPECT_EQ(v.tuple_idx[1], 1);
    EXPECT_EQ(v.local_i[1], 2);
    EXPECT_EQ(v.local_j[1], 3);
    // 3-4 (4 is a halo ghost — local row 4)
    EXPECT_EQ(v.tuple_idx[2], 2);
    EXPECT_EQ(v.local_i[2], 3);
    EXPECT_EQ(v.local_j[2], 4);
}

TEST(BondResolvedLocals, GenerationKeyCachesOwnerFilter) {
    auto frame = make_chain(8);
    frame.bonds.atomi = xt::xtensor<int, 1>{0, 2};
    frame.bonds.atomj = xt::xtensor<int, 1>{1, 3};
    frame.bonds.type = xt::zeros<int>({2});

    cpu::md::mpi::DomainDecomposition dd;
    dd.n_local = 4;
    dd.n_ghost = 0;
    dd.capacity = 8;
    dd.neighbor_generation = 1;
    dd.tag_to_local = {{0, 0}, {1, 1}, {2, 2}, {3, 3}};

    cpu::Context ctx;
    ctx.mpi_size = 2;
    ctx.dd = &dd;
    cpu::LaunchContext lctx{ctx, frame, 0, 0};
    TestBond b;
    auto v1 = b.resolve(lctx, 2, frame.bonds.atomi.data(),
                        frame.bonds.atomj.data());
    EXPECT_EQ(v1.n_work, 2);

    // Now mutate tag_to_local to drop tag 3. Since the owner filter is
    // cached by generation, the second call (same generation) keeps the
    // (2-3) bond in the owner list, but the per-step resolution will
    // drop it because partner-lookup misses.
    dd.tag_to_local.erase(3);
    auto v2 = b.resolve(lctx, 2, frame.bonds.atomi.data(),
                        frame.bonds.atomj.data());
    EXPECT_EQ(v2.n_work, 1) << "(2-3) drops at per-step resolve, not cache";
    EXPECT_EQ(v2.tuple_idx[0], 0);

    // Bumping the generation forces the owner filter to rebuild. Owner
    // for bond (2-3) is tag 2 which is still local, so it stays in
    // owned_; but per-step resolve still drops it.
    dd.neighbor_generation = 2;
    auto v3 = b.resolve(lctx, 2, frame.bonds.atomi.data(),
                        frame.bonds.atomj.data());
    EXPECT_EQ(v3.n_work, 1);
}
