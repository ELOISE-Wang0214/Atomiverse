/// @file cpu/scf/basis_data/def2_svp.h
/// @brief def2-svp basis-set table (auto-generated; do not edit).
///
/// Emitted by `tools/gen_basis_data.py` from PySCF's `mol._basis`
/// dictionary, so numerical values match the PySCF ground-truth path
/// bit-for-bit. Contraction coefficients use the unnormalized
/// convention (physical GTO normalization is applied by
/// `ShellList` constructor).
///
/// Supported elements: H, C, N, O, F.
///
/// Regenerate:
///   python tools/gen_basis_data.py def2-svp def2_svp H C N O F
#pragma once

#include <vector>

#include "atomiverse/cpu/scf/basis_data/basis_table.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace basis_data {

/// @brief Get the def2-svp table, indexed by atomic number Z.
/// @return Reference to a static table with AtomBasis entries for
///         Z=1..9; position 0 and any unsupported Z are empty.
inline const std::vector<AtomBasis>& def2_svp() {
    static const std::vector<AtomBasis> data = []() {
        return std::vector<AtomBasis>{
            AtomBasis{}, // Z=0 unused
            // Z=1 (H): 3 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {13.010701, 1.9622572, 0.44453796},
                                {0.019682158, 0.13796524, 0.47831935}},
                ContractedShell{/*l=*/0 /* s */, {0.12194962}, {1.0}},
                ContractedShell{/*l=*/1 /* p */, {0.8}, {1.0}},
            }},
            AtomBasis{}, // Z=2 unsupported
            AtomBasis{}, // Z=3 unsupported
            AtomBasis{}, // Z=4 unsupported
            AtomBasis{}, // Z=5 unsupported
            // Z=6 (C): 6 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {1238.4016938, 186.29004992, 42.251176346,
                                 11.676557932, 3.5930506482},
                                {0.0054568832082, 0.040638409211, 0.18025593888,
                                 0.46315121755, 0.44087173314}},
                ContractedShell{/*l=*/0 /* s */, {0.40245147363}, {1.0}},
                ContractedShell{/*l=*/0 /* s */, {0.13090182668}, {1.0}},
                ContractedShell{/*l=*/1 /* p */,
                                {9.4680970621, 2.0103545142, 0.54771004707},
                                {0.038387871728, 0.21117025112, 0.51328172114}},
                ContractedShell{/*l=*/1 /* p */, {0.15268613795}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {0.8}, {1.0}},
            }},
            // Z=7 (N): 6 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {1712.8415853, 257.64812677, 58.458245853,
                                 16.198367905, 5.0052600809},
                                {-0.0053934125305, -0.040221581118,
                                 -0.1793114499, -0.46376317823,
                                 -0.44171422662}},
                ContractedShell{/*l=*/0 /* s */, {0.58731856571}, {1.0}},
                ContractedShell{/*l=*/0 /* s */, {0.18764592253}, {1.0}},
                ContractedShell{
                    /*l=*/1 /* p */,
                    {13.571470233, 2.9257372874, 0.79927750754},
                    {-0.040072398852, -0.21807045028, -0.51294466049}},
                ContractedShell{/*l=*/1 /* p */, {0.21954348034}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {1.0}, {1.0}},
            }},
            // Z=8 (O): 6 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {2266.1767785, 340.87010191, 77.363135167,
                                 21.47964494, 6.6589433124},
                                {-0.0053431809926, -0.03989003923,
                                 -0.17853911985, -0.46427684959,
                                 -0.44309745172}},
                ContractedShell{/*l=*/0 /* s */, {0.80975975668}, {1.0}},
                ContractedShell{/*l=*/0 /* s */, {0.25530772234}, {1.0}},
                ContractedShell{/*l=*/1 /* p */,
                                {17.721504317, 3.863550544, 1.0480920883},
                                {0.043394573193, 0.23094120765, 0.51375311064}},
                ContractedShell{/*l=*/1 /* p */, {0.27641544411}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {1.2}, {1.0}},
            }},
            // Z=9 (F): 6 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {2894.832599, 435.4193912, 98.843328866,
                                 27.485198001, 8.5405498171},
                                {-0.0053408255515, -0.039904258866,
                                 -0.17912768038, -0.46758090825,
                                 -0.4465313102}},
                ContractedShell{/*l=*/0 /* s */, {1.0654578038}, {1.0}},
                ContractedShell{/*l=*/0 /* s */, {0.33247346748}, {1.0}},
                ContractedShell{
                    /*l=*/1 /* p */,
                    {22.696633924, 4.9872339257, 1.3491613954},
                    {-0.045212874436, -0.23754317067, -0.51287353587}},
                ContractedShell{/*l=*/1 /* p */, {0.34829881977}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {1.4}, {1.0}},
            }},
        };
    }();
    return data;
}

} // namespace basis_data
} // namespace scf
} // namespace cpu
} // namespace ATV
