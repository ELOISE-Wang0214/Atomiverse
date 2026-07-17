/// @file cpu/scf/basis_data/def2_tzvp.h
/// @brief def2-tzvp basis-set table (auto-generated; do not edit).
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
///   python tools/gen_basis_data.py def2-tzvp def2_tzvp H C N O F
#pragma once

#include <vector>

#include "atomiverse/cpu/scf/basis_data/basis_table.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace basis_data {

/// @brief Get the def2-tzvp table, indexed by atomic number Z.
/// @return Reference to a static table with AtomBasis entries for
///         Z=1..9; position 0 and any unsupported Z are empty.
inline const std::vector<AtomBasis>& def2_tzvp() {
    static const std::vector<AtomBasis> data = []() {
        return std::vector<AtomBasis>{
            AtomBasis{}, // Z=0 unused
            // Z=1 (H): 4 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {34.061341, 5.1235746, 1.1646626},
                                {0.0060251978, 0.045021094, 0.20189726}},
                ContractedShell{/*l=*/0 /* s */, {0.32723041}, {1.0}},
                ContractedShell{/*l=*/0 /* s */, {0.10307241}, {1.0}},
                ContractedShell{/*l=*/1 /* p */, {0.8}, {1.0}},
            }},
            AtomBasis{}, // Z=2 unsupported
            AtomBasis{}, // Z=3 unsupported
            AtomBasis{}, // Z=4 unsupported
            AtomBasis{}, // Z=5 unsupported
            // Z=6 (C): 11 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {13575.349682, 2035.233368, 463.22562359,
                                 131.20019598, 42.853015891, 15.584185766},
                                {0.00022245814352, 0.0017232738252,
                                 0.0089255715314, 0.035727984502, 0.11076259931,
                                 0.24295627626}},
                ContractedShell{/*l=*/0 /* s */,
                                {6.2067138508, 2.5764896527},
                                {0.41440263448, 0.23744968655}},
                ContractedShell{/*l=*/0 /* s */, {0.57696339419}, {1.0}},
                ContractedShell{/*l=*/0 /* s */, {0.22972831358}, {1.0}},
                ContractedShell{/*l=*/0 /* s */, {0.095164440028}, {1.0}},
                ContractedShell{
                    /*l=*/1 /* p */,
                    {34.697232244, 7.9582622826, 2.3780826883, 0.81433208183},
                    {0.0053333657805, 0.035864109092, 0.14215873329,
                     0.34270471845}},
                ContractedShell{/*l=*/1 /* p */, {0.28887547253}, {1.0}},
                ContractedShell{/*l=*/1 /* p */, {0.10056823671}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {1.097}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {0.318}, {1.0}},
                ContractedShell{/*l=*/3 /* f */, {0.761}, {1.0}},
            }},
            // Z=7 (N): 11 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {19730.800647, 2957.8958745, 673.22133595,
                                 190.68249494, 62.295441898, 22.654161182},
                                {0.00021887984991, 0.0016960708803,
                                 0.0087954603538, 0.035359382605, 0.11095789217,
                                 0.24982972552}},
                ContractedShell{/*l=*/0 /* s */,
                                {8.9791477428, 3.686300237},
                                {0.40623896148, 0.24338217176}},
                ContractedShell{/*l=*/0 /* s */, {0.84660076805}, {1.0}},
                ContractedShell{/*l=*/0 /* s */, {0.33647133771}, {1.0}},
                ContractedShell{/*l=*/0 /* s */, {0.13647653675}, {1.0}},
                ContractedShell{
                    /*l=*/1 /* p */,
                    {49.20038051, 11.346790537, 3.4273972411, 1.1785525134},
                    {0.0055552416751, 0.038052379723, 0.14953671029,
                     0.3494930523}},
                ContractedShell{/*l=*/1 /* p */, {0.41642204972}, {1.0}},
                ContractedShell{/*l=*/1 /* p */, {0.14260826011}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {1.654}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {0.469}, {1.0}},
                ContractedShell{/*l=*/3 /* f */, {1.093}, {1.0}},
            }},
            // Z=8 (O): 11 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {27032.382631, 4052.3871392, 922.3272271,
                                 261.24070989, 85.354641351, 31.035035245},
                                {0.00021726302465, 0.0016838662199,
                                 0.0087395616265, 0.035239968808, 0.11153519115,
                                 0.25588953961}},
                ContractedShell{/*l=*/0 /* s */,
                                {12.260860728, 4.9987076005},
                                {0.39768730901, 0.2462784943}},
                ContractedShell{/*l=*/0 /* s */, {1.1703108158}, {1.0}},
                ContractedShell{/*l=*/0 /* s */, {0.46474740994}, {1.0}},
                ContractedShell{/*l=*/0 /* s */, {0.18504536357}, {1.0}},
                ContractedShell{
                    /*l=*/1 /* p */,
                    {63.274954801, 14.627049379, 4.4501223456, 1.5275799647},
                    {0.0060685103418, 0.041912575824, 0.16153841088,
                     0.35706951311}},
                ContractedShell{/*l=*/1 /* p */, {0.52935117943}, {1.0}},
                ContractedShell{/*l=*/1 /* p */, {0.1747842127}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {2.314}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {0.645}, {1.0}},
                ContractedShell{/*l=*/3 /* f */, {1.428}, {1.0}},
            }},
            // Z=9 (F): 11 shells
            AtomBasis{{
                ContractedShell{/*l=*/0 /* s */,
                                {35479.100441, 5318.4728983, 1210.4810975,
                                 342.8551814, 112.01943181, 40.714740248},
                                {0.00021545014888, 0.0016700686527,
                                 0.0086733211476, 0.035049933175, 0.11165320133,
                                 0.25988506647}},
                ContractedShell{/*l=*/0 /* s */,
                                {16.039678111, 6.503818674},
                                {0.3942296688, 0.24998238551}},
                ContractedShell{/*l=*/0 /* s */, {1.5440477509}, {1.0}},
                ContractedShell{/*l=*/0 /* s */, {0.61223452862}, {1.0}},
                ContractedShell{/*l=*/0 /* s */, {0.24027979698}, {1.0}},
                ContractedShell{
                    /*l=*/1 /* p */,
                    {80.233900483, 18.594010743, 5.6867902653, 1.9511006294},
                    {0.0063685999134, 0.04430314353, 0.16867248708,
                     0.36166346255}},
                ContractedShell{/*l=*/1 /* p */, {0.66970211298}, {1.0}},
                ContractedShell{/*l=*/1 /* p */, {0.2165130041}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {3.107}, {1.0}},
                ContractedShell{/*l=*/2 /* d */, {0.855}, {1.0}},
                ContractedShell{/*l=*/3 /* f */, {1.917}, {1.0}},
            }},
        };
    }();
    return data;
}

} // namespace basis_data
} // namespace scf
} // namespace cpu
} // namespace ATV
