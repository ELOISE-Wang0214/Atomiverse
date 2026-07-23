#include <gtest/gtest.h>
#include "atomiverse/cpu/minimize/anisotropic_cell_relax.h"
using namespace ATV;
TEST(AnisotropicCellRelax, Name) {
    cpu::minimizer::AnisotropicCellRelax a(0.0, 0.05);
    EXPECT_EQ(a.name(), "anisotropic_cell_relax");
}
