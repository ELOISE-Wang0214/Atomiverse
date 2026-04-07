/// @file cpu/recorder/thermo.cpp
/// @brief CPU thermodynamic output: temperature calculation + formatted printing.

#include "atomiverse/cpu/recorder/thermo.h"
#include "atomiverse/core/log.h"
#include "atomiverse/core/units.h"

#include <spdlog/fmt/fmt.h>

namespace ATV {
namespace cpu {

void thermo_print(int step, double pe, double ke, int n_atoms,
                  bool& header_printed, std::ofstream& file)
{
    // Temperature: T = 2*KE / (ndof * kB),  ndof = 3N - 3
    const int ndof = 3 * n_atoms - 3;
    const double temp = (ndof > 0)
        ? 2.0 * ke / (static_cast<double>(ndof) * kBoltzmann)
        : 0.0;
    const double total = pe + ke;

    if (!header_printed) {
        ATV_INFO("{:>8s} {:>12s} {:>14s} {:>14s} {:>14s}",
                 "Step", "Temp", "PotEng", "KinEng", "TotEng");
        header_printed = true;
    }
    ATV_INFO("{:8d} {:12.4f} {:14.6f} {:14.6f} {:14.6f}",
             step,
             static_cast<double>(temp),
             static_cast<double>(pe),
             static_cast<double>(ke),
             static_cast<double>(total));

    if (file.is_open()) {
        file << fmt::format("{} {:.6f} {:.6f} {:.6f} {:.6f}\n",
                            step,
                            static_cast<double>(temp),
                            static_cast<double>(pe),
                            static_cast<double>(ke),
                            static_cast<double>(total));
    }
}

} // namespace cpu
} // namespace ATV
