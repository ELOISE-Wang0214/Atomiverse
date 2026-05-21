/// @file cpu/scf/basis_data/pople_631gs.h
/// @brief 6-31g* basis-set table (auto-generated; do not edit).
///
/// Emitted by `tools/gen_basis_data.py` from PySCF's `mol._basis`
/// dictionary, so numerical values match the PySCF ground-truth path
/// bit-for-bit. Contraction coefficients use the unnormalized
/// convention (physical GTO normalization is applied by
/// `build_shell_list`).
///
/// Supported elements: H, C, N, O, F.
///
/// Regenerate:
///   python tools/gen_basis_data.py 6-31g* pople_631gs H C N O F
#pragma once

#include <vector>

#include "atomiverse/cpu/scf/basis_data/basis_table.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace basis_data {

/// @brief Get the 6-31g* table, indexed by atomic number Z.
/// @return Reference to a static table with AtomBasis entries for
///         Z=1..9; position 0 and any unsupported Z are empty.
inline const std::vector<AtomBasis>& pople_631gs() {
    static const std::vector<AtomBasis> data = []() {
        return std::vector<AtomBasis>{
            AtomBasis{}, // Z=0 unused
            // Z=1 (H): 2 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {18.731137, 2.8253937, 0.6401217},
                                {0.0334946, 0.23472695, 0.81375733}},
                ContractedShell{/*l=*/0 /* s */, {0.1612778}, {1.0}},
            }},
            AtomBasis{}, // Z=2 unsupported
            AtomBasis{}, // Z=3 unsupported
            AtomBasis{}, // Z=4 unsupported
            AtomBasis{}, // Z=5 unsupported
            // Z=6 (C): 6 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {3047.5249, 457.36951, 103.94869, 29.210155,
                                 9.286663, 3.163927},
                                {0.0018347, 0.0140373, 0.0688426, 0.2321844,
                                 0.4679413, 0.362312}},
                ContractedShell{/*l=*/0 /* s */,
                                {7.8682724, 1.8812885, 0.5442493},
                                {-0.1193324, -0.1608542, 1.1434564}},
                ContractedShell{/*l=*/0 /* s */, {0.1687144}, {1.0}},
                ContractedShell{/*l=*/1 /* p */,
                                {7.8682724, 1.8812885, 0.5442493},
                                {0.0689991, 0.316424, 0.7443083}},
                ContractedShell{/*l=*/1 /* p */, {0.1687144}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {0.8}, {1.0}},
            }},
            // Z=7 (N): 6 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {4173.511, 627.4579, 142.9021, 40.23433,
                                 12.82021, 4.390437},
                                {0.0018348, 0.013995, 0.068587, 0.232241,
                                 0.46907, 0.360455}},
                ContractedShell{/*l=*/0 /* s */,
                                {11.626358, 2.71628, 0.772218},
                                {-0.114961, -0.169118, 1.145852}},
                ContractedShell{/*l=*/0 /* s */, {0.2120313}, {1.0}},
                ContractedShell{/*l=*/1 /* p */,
                                {11.626358, 2.71628, 0.772218},
                                {0.06758, 0.323907, 0.740895}},
                ContractedShell{/*l=*/1 /* p */, {0.2120313}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {0.8}, {1.0}},
            }},
            // Z=8 (O): 6 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {5484.6717, 825.23495, 188.04696, 52.9645,
                                 16.89757, 5.7996353},
                                {0.0018311, 0.0139501, 0.0684451, 0.2327143,
                                 0.470193, 0.3585209}},
                ContractedShell{/*l=*/0 /* s */,
                                {15.539616, 3.5999336, 1.0137618},
                                {-0.1107775, -0.1480263, 1.130767}},
                ContractedShell{/*l=*/0 /* s */, {0.2700058}, {1.0}},
                ContractedShell{/*l=*/1 /* p */,
                                {15.539616, 3.5999336, 1.0137618},
                                {0.0708743, 0.3397528, 0.7271586}},
                ContractedShell{/*l=*/1 /* p */, {0.2700058}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {0.8}, {1.0}},
            }},
            // Z=9 (F): 6 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {7001.71309, 1051.36609, 239.28569, 67.3974453,
                                 21.5199573, 7.4031013},
                                {0.0018196169, 0.0139160796, 0.0684053245,
                                 0.23318576, 0.471267439, 0.356618546}},
                ContractedShell{/*l=*/0 /* s */,
                                {20.8479528, 4.80830834, 1.34406986},
                                {-0.108506975, -0.146451658, 1.12868858}},
                ContractedShell{/*l=*/0 /* s */, {0.358151393}, {1.0}},
                ContractedShell{/*l=*/1 /* p */,
                                {20.8479528, 4.80830834, 1.34406986},
                                {0.0716287243, 0.345912103, 0.722469957}},
                ContractedShell{/*l=*/1 /* p */, {0.358151393}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {0.8}, {1.0}},
            }},
        };
    }();
    return data;
}

} // namespace basis_data
} // namespace scf
} // namespace cpu
} // namespace ATV
