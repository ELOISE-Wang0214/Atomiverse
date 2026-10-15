/// @file tests/cpu/semiempirical/test_am1_bcc_charges.cpp
/// @brief Unit gates for AM1BCCCharges: Frame bonds → local graph, the
///        missing-topology contract, the parameter-set names. No SCF runs.

#include <array>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/am1_bcc_charges.h"
#include "atomiverse/cpu/semiempirical/internal/bcc_graph.h"

using namespace ATV;
using namespace ATV::cpu::semiempirical;

namespace {

/// Four atoms with non-contiguous global tags 10..13 (row i has tag 10 + i).
cpu::Frame tagged_chain() {
    cpu::Frame f(4);
    for (int i = 0; i < 4; ++i) {
        f.span(kAtom_type)[static_cast<std::size_t>(i)] = 6;
        f.span(kAtom_tag)[static_cast<std::size_t>(i)] = 10 + i;
    }
    f.clear_tag_index();
    f.rebuild_tag_index(0, f.nlocal());
    return f;
}

} // namespace

TEST(AM1BCCChargesLocalBonds, MapsTagsToRowsOncePerBond) {
    auto f = tagged_chain();
    f.add_bond(12, 13, 0);
    f.add_bond(10, 11, 0);
    f.add_bond(11, 12, 0);
    f.add_bond(11, 10, 0); // the same bond again: must come back once

    const auto bonds = internal::local_bonds(f);
    const std::vector<std::array<int, 2>> want = {{0, 1}, {1, 2}, {2, 3}};
    EXPECT_EQ(bonds, want);
}

TEST(AM1BCCChargesLocalBonds, StoredOnBothEndsComesBackOnce) {
    auto f = tagged_chain();
    f.add_bond(10, 11, 0); // owner row 0 → partner tag 11
    // Newton-off storage: the same bond also listed under row 1.
    f.grow_bond_per_atom(2);
    f.span(kAtom_bond_atom)(1, 0) = 10;
    f.span(kAtom_num_bond)[1] = 1;

    const auto bonds = internal::local_bonds(f);
    ASSERT_EQ(bonds.size(), 1u);
    EXPECT_EQ(bonds[0], (std::array<int, 2>{0, 1}));
}

TEST(AM1BCCChargesLocalBonds, NoTopologyIsEmpty) {
    const auto f = tagged_chain();
    EXPECT_TRUE(internal::local_bonds(f).empty());
}

TEST(AM1BCCChargesLocalBonds, PartnerOffFrameThrows) {
    auto f = tagged_chain();
    f.add_bond(10, 99, 0);
    EXPECT_THROW(internal::local_bonds(f), std::invalid_argument);
}

TEST(AM1BCCCharges, MultiAtomFrameWithoutBondsThrows) {
    auto f = tagged_chain();
    AM1BCCCharges model;
    EXPECT_THROW(model.assign(f), std::invalid_argument);
}

TEST(AM1BCCCharges, ParameterSetNamesAreTheMolrsNames) {
    EXPECT_EQ(internal::molrs_name(BccParameterSet::Bcc), "bcc");
    EXPECT_EQ(internal::molrs_name(BccParameterSet::Abcg2), "abcg2");
    EXPECT_EQ(AM1BCCCharges{}.parameter_set(), BccParameterSet::Bcc);
    EXPECT_EQ((AM1BCCCharges{AM1{}, BccParameterSet::Abcg2}.parameter_set()),
              BccParameterSet::Abcg2);
}
