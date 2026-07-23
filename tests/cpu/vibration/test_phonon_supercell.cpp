/// @file tests/cpu/vibration/test_phonon_supercell.cpp
/// @brief Unit tests for expand_supercell and cpu::vib::Phonon wiring.

#include <cmath>
#include <stdexcept>

#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/sw.h"
#include "atomiverse/cpu/vibration/phonon.h"
#include "atomiverse/cpu/vibration/supercell.h"
#include "atomiverse/molrs.hpp"

using namespace ATV;
using namespace ATV::cpu::vib;

namespace {

cpu::Frame diamond_conventional(double a, int n_cell) {
    const double frac[8][3] = {
        {0, 0, 0},     {0.25, 0.25, 0.25}, {0.5, 0.5, 0}, {0.75, 0.75, 0.25},
        {0.5, 0, 0.5}, {0.75, 0.25, 0.75}, {0, 0.5, 0.5}, {0.25, 0.75, 0.75}};
    const int n = 8 * n_cell * n_cell * n_cell;
    cpu::Frame f(n);
    int idx = 0;
    for (int ix = 0; ix < n_cell; ++ix)
        for (int iy = 0; iy < n_cell; ++iy)
            for (int iz = 0; iz < n_cell; ++iz)
                for (int k = 0; k < 8; ++k) {
                    f.span(kAtom_x)[idx] = static_cast<StorageReal>(
                        (double(ix) + frac[k][0]) * a);
                    f.span(kAtom_y)[idx] = static_cast<StorageReal>(
                        (double(iy) + frac[k][1]) * a);
                    f.span(kAtom_z)[idx] = static_cast<StorageReal>(
                        (double(iz) + frac[k][2]) * a);
                    f.span(kAtom_type)[idx] = 14;
                    f.span(kAtom_mass)[idx] = 28.0855;
                    f.span(kAtom_tag)[idx] = static_cast<TagInt>(idx + 1);
                    ++idx;
                }
    f.simbox = cpu::Box::cubic(a * n_cell);
    return f;
}

cpu::Frame diamond_primitive(double a) {
    cpu::Frame f(2);
    f.span(kAtom_type)[0] = 14;
    f.span(kAtom_type)[1] = 14;
    f.span(kAtom_mass)[0] = 28.0855;
    f.span(kAtom_mass)[1] = 28.0855;
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_x)[1] = a / 4.0;
    f.span(kAtom_y)[1] = a / 4.0;
    f.span(kAtom_z)[1] = a / 4.0;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    const double h = a / 2.0;
    Mat3 cell;
    cell(0, 0) = 0;
    cell(0, 1) = h;
    cell(0, 2) = h;
    cell(1, 0) = h;
    cell(1, 1) = 0;
    cell(1, 2) = h;
    cell(2, 0) = h;
    cell(2, 1) = h;
    cell(2, 2) = 0;
    f.simbox.set_cell(cell);
    return f;
}

} // namespace

TEST(PhononSupercell, ExpandCountsAndCell) {
    cpu::Frame unit = diamond_primitive(5.431);
    cpu::Frame sc = expand_supercell(unit, 2, 2, 2);
    EXPECT_EQ(sc.nlocal(), 16);
    EXPECT_NEAR(sc.simbox.volume(), unit.simbox.volume() * 8.0, 1e-9);
    const Mat3& c = unit.simbox.cell();
    const double ax = c(0, 0), ay = c(0, 1), az = c(0, 2);
    bool found = false;
    for (int i = 0; i < sc.nlocal(); ++i) {
        const double dx = static_cast<double>(sc.span(kAtom_x)[i]) - ax;
        const double dy = static_cast<double>(sc.span(kAtom_y)[i]) - ay;
        const double dz = static_cast<double>(sc.span(kAtom_z)[i]) - az;
        if (std::abs(dx) < 1e-12 && std::abs(dy) < 1e-12 && std::abs(dz) < 1e-12)
            found = true;
    }
    EXPECT_TRUE(found);
}

TEST(PhononSupercell, ExpandRejectsBadFactors) {
    cpu::Frame unit = diamond_primitive(5.431);
    EXPECT_THROW(expand_supercell(unit, 0, 1, 1), std::invalid_argument);
}

TEST(PhononSupercell, BuildRequiresForceKind) {
    Phonon phonon(1, 1, 1);
    EXPECT_THROW(phonon.build(), std::runtime_error);
}

TEST(PhononSupercell, SwAcousticNearZero) {
    // 2×2×2 conventional diamond (64 atoms): box L=10.86 Å > 2*r_cut for SW.
    // Unit-scale: only check expand + acoustic sum rule (not full e2e pin).
    const double a0 = 5.431;
    cpu::Frame unit = diamond_conventional(a0, /*n_cell=*/2);
    ASSERT_EQ(unit.nlocal(), 64);
    cpu::SWPotential pot(molrs::Element::Si, "stillinger-weber");
    const double cutoff = pot.cutoff();
    ASSERT_GT(unit.simbox.length(0), 2.0 * cutoff);

    Phonon phonon(1, 1, 1);
    // always rebuild (check=false): displacements must re-enumerate pairs.
    phonon.set_neighbor(cpu::md::NeighborList<>{cutoff, 0.3, 1, 0, false});
    phonon.on(ATV::vib::Phase::force).add(cpu::StillingerWeber{pot});
    phonon.build();
    const PhononResult r = phonon.run(unit);

    ASSERT_EQ(r.n_atoms, 64);
    ASSERT_EQ(static_cast<int>(r.frequencies.size()), 192);

    double max_ac = 0.0;
    for (int k = 0; k < 3; ++k)
        max_ac = std::max(max_ac, std::abs(r.frequencies(k)));
    EXPECT_LT(max_ac, 50.0) << "acoustic |ω| max=" << max_ac;
    EXPECT_GT(r.frequencies(191), 100.0);
}
