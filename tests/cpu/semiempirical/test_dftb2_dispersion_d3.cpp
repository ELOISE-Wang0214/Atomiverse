/// @file tests/cpu/semiempirical/test_dftb2_dispersion_d3.cpp
/// @brief Unit: DFTB2::d3 folds D3(BJ) into total energy (water dimer).
#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"

using namespace ATV;

namespace {

cpu::Frame water_dimer() {
    // Spec geometry (Å); cubic box 20 Å (isolated).
    cpu::Frame f(6);
    const double xyz[6][3] = {
        {0.000, 0.000, 0.000},  {0.757, 0.000, 0.586},  {-0.757, 0.000, 0.586},
        {0.000, 0.000, 2.900},  {0.757, 0.000, 3.486},  {-0.757, 0.000, 3.486},
    };
    const int Z[6] = {8, 1, 1, 8, 1, 1};
    for (int i = 0; i < 6; ++i) {
        f.span(kAtom_type)[i] = Z[i];
        f.span(kAtom_x)[i] = static_cast<StorageReal>(xyz[i][0]);
        f.span(kAtom_y)[i] = static_cast<StorageReal>(xyz[i][1]);
        f.span(kAtom_z)[i] = static_cast<StorageReal>(xyz[i][2]);
        f.span(kAtom_mass)[i] = (Z[i] == 8) ? 15.999 : 1.008;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    }
    f.box = cpu::Box::cubic(20.0);
    f.meta.get(kEnergy) = 0;
    return f;
}

} // namespace

TEST(DFTB2DispersionD3, OffMeansZeroDisp) {
    auto f = water_dimer();
    cpu::Context ctx;
    cpu::semiempirical::DFTB2 dftb;
    dftb.d3 = false;
    dftb.with_forces = false;
    dftb.build({ctx, f});
    dftb.launch({ctx, f, 0, 0});
    EXPECT_NEAR(dftb.last_e_disp(), 0.0, 1e-15);
    EXPECT_TRUE(std::isfinite(dftb.last_energy()));
}

TEST(DFTB2DispersionD3, OnAddsAttractiveDisp) {
    auto f = water_dimer();
    cpu::Context ctx;
    cpu::semiempirical::DFTB2 dftb;
    dftb.d3 = true;
    dftb.with_forces = false;
    dftb.build({ctx, f});
    dftb.launch({ctx, f, 0, 0});
    EXPECT_LT(dftb.last_e_disp(), 0.0);
    EXPECT_TRUE(std::isfinite(dftb.last_e_disp()));
    EXPECT_NEAR(f.meta.get(kEnergy), dftb.last_energy(), 1e-12);
}

TEST(DFTB2DispersionD3, TotDiffEqualsEDisp) {
    cpu::Context ctx;
    auto f_off = water_dimer();
    cpu::semiempirical::DFTB2 off;
    off.d3 = false;
    off.with_forces = false;
    off.build({ctx, f_off});
    off.launch({ctx, f_off, 0, 0});
    const double e_off = off.last_energy();

    auto f_on = water_dimer();
    cpu::semiempirical::DFTB2 on;
    on.d3 = true;
    on.with_forces = false;
    on.build({ctx, f_on});
    on.launch({ctx, f_on, 0, 0});
    const double e_on = on.last_energy();
    const double ed = on.last_e_disp();

    EXPECT_NEAR(e_on - e_off, ed, 1e-9);
}
