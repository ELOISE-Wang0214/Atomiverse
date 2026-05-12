/// @file cpu/pw/internal/upf_reader.h
/// @brief UPF v1 pseudopotential file parser.
#pragma once

#include <string>
#include <vector>
#include "atomiverse/core/real.h"

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

struct UPFData {
    int z_val = 0;
    int l_max = 0;
    std::vector<Real> r;          // radial mesh [Bohr]
    std::vector<Real> v_loc;      // local potential [Hartree]
    std::vector<Real> rho_atom;   // atomic pseudocharge density
    struct KBProjector {
        int l;
        std::vector<Real> beta;   // radial beta function
        Real d;                    // D coefficient for this l
    };
    std::vector<KBProjector> kb;  // Kleinman-Bylander projectors
};

UPFData read_upf_v1(const std::string& path);

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
