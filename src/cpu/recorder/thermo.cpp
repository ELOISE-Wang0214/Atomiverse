/// @file cpu/recorder/thermo.cpp
/// @brief CPU thermodynamic output: temperature calculation + formatted printing.

#include "atomiverse/cpu/recorder/thermo.h"
#include "atomiverse/core/units.h"

#include <cstdio>

namespace ATV {
namespace cpu {

void thermo_print(int step, Real pe, Real ke, Int n_atoms,
                  bool& header_printed, std::ofstream& file)
{
    // Temperature: T = 2*KE / (ndof * kB),  ndof = 3N - 3
    Int ndof = 3 * n_atoms - 3;
    Real temp = (ndof > 0) ? Real(2) * ke / (Real(ndof) * kBoltzmann) : Real(0);
    Real total = pe + ke;

    if (!header_printed) {
        std::printf("%8s %12s %14s %14s %14s\n",
                    "Step", "Temp", "PotEng", "KinEng", "TotEng");
        header_printed = true;
    }
    std::printf("%8d %12.4f %14.6f %14.6f %14.6f\n",
                step,
                static_cast<double>(temp),
                static_cast<double>(pe),
                static_cast<double>(ke),
                static_cast<double>(total));

    if (file.is_open()) {
        char buf[128];
        std::snprintf(buf, sizeof(buf), "%d %.6f %.6f %.6f %.6f\n",
                      step,
                      static_cast<double>(temp),
                      static_cast<double>(pe),
                      static_cast<double>(ke),
                      static_cast<double>(total));
        file << buf;
    }
}

} // namespace cpu
} // namespace ATV
