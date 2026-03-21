/// @file cpu/recorder/trajectory.cpp
/// @brief CPU plain-XYZ trajectory frame writer.

#include "atomiverse/cpu/recorder/trajectory.h"

#include <iomanip>

namespace ATV {
namespace cpu {

void trajectory_append(std::ofstream& os,
                       const Int* type_id,
                       const Real* x, const Real* y, const Real* z,
                       Int n_atoms, int step)
{
    os << n_atoms << "\n";
    os << "Step " << step << "\n";
    os << std::fixed << std::setprecision(6);
    for (Int i = 0; i < n_atoms; ++i) {
        os << type_id[i] << " "
           << x[i] << " " << y[i] << " " << z[i] << "\n";
    }
}

} // namespace cpu
} // namespace ATV
