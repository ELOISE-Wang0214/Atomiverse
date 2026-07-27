/// @file tests/cpu/communication/test_halo_refresh_charges.cpp
/// @brief Unit tests for `comm::HaloRefresh{Field::charges}` (serial paths).
///
/// Multi-rank transport correctness lives in
/// `tests/mpi/test_halo_refresh_charges.cpp` (launched under mpirun). Here:
/// naming, and the single-rank launch being a strict no-op on the charge
/// column.

#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/communication/halo_exchange.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"

using namespace ATV;

namespace {

cpu::Frame two_charged_atoms() {
    cpu::Frame f(2);
    for (int i = 0; i < 2; ++i) {
        f.span(kAtom_x)[i] = 1.0 + i;
        f.span(kAtom_y)[i] = 2.0;
        f.span(kAtom_z)[i] = 3.0;
        f.span(kAtom_charge)[i] = (i == 0) ? 0.25 : -0.25;
        f.span(kAtom_type)[i] = 1;
        f.span(kAtom_mass)[i] = 1.0;
        f.span(kAtom_tag)[i] = i + 1;
    }
    f.box = cpu::Box::cubic(10.0);
    return f;
}

} // namespace

TEST(HaloRefreshCharges, Name) {
    using F = cpu::comm::HaloRefresh::Field;
    EXPECT_EQ(cpu::comm::HaloRefresh{F::charges}.name(), "halo_refresh_q");
    EXPECT_EQ(cpu::comm::HaloRefresh{F::positions}.name(), "halo_refresh_x");
    EXPECT_EQ(cpu::comm::HaloRefresh{F::velocities}.name(), "halo_refresh_v");
}

TEST(HaloRefreshCharges, SingleRankLaunchIsNoOp) {
    auto f = two_charged_atoms();
    cpu::Context ctx; // no DD bound; mpi_size = 1 default
    cpu::comm::HaloRefresh refresh{cpu::comm::HaloRefresh::Field::charges};
    refresh.build({ctx, f});
    EXPECT_NO_THROW(refresh.launch({ctx, f, /*stage=*/0, /*step=*/0}));
    EXPECT_DOUBLE_EQ(static_cast<double>(f.span(kAtom_charge)[0]), 0.25);
    EXPECT_DOUBLE_EQ(static_cast<double>(f.span(kAtom_charge)[1]), -0.25);
}
