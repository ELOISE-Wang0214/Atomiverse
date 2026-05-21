/// @file cpu/pw/internal/upf_reader.h
/// @brief UPF v1 pseudopotential file parser.
#pragma once

#include "atomiverse/core/real.h"
#include <string>
#include <vector>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

struct UPFData {
    int z_val = 0;
    int l_max = 0;
    std::vector<Real> r;        // radial mesh [Bohr]
    std::vector<Real> v_loc;    // local potential [Hartree]
    std::vector<Real> rho_atom; // atomic pseudocharge density
    struct KBProjector {
        int l;
        std::vector<Real> beta; // radial beta function
        Real d;                 // D coefficient for this l
    };
    std::vector<KBProjector> kb; // Kleinman-Bylander projectors
    std::vector<Real> d_ij; // Full D_ij matrix (n_proj × n_proj, row-major)
    std::vector<Real>
        q_ij; // Augmentation charge matrix (n_proj × n_proj, row-major)
    int n_proj() const { return static_cast<int>(kb.size()); }

    // Atomic wavefunctions for initial guess (atomic_wfc_mod.f90)
    std::vector<std::vector<Real>> chi; // chi[nb][ir]: atomic wfc radial part
    std::vector<int> lchi;              // angular momentum for each chi
    std::vector<Real> oc;               // occupation numbers
    int nwfc = 0;                       // number of atomic wavefunctions
};

UPFData read_upf_v1(const std::string& path);

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
