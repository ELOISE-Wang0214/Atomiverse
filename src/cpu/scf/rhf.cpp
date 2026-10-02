/// @file cpu/scf/rhf.cpp
/// @brief RHF constructor defaults and occupation resolution. The shared
///        driver loop lives in mean_field_method.cpp. Seeger–Pople stability
///        is `Stability<RHF>` (`response.cpp`), not a method here.

#include "atomiverse/cpu/scf/rhf.h"

namespace ATV {
namespace cpu {
namespace scf {

RHF::RHF() {
    spin_ = SpinType::Restricted;
    method_ = MeanFieldType::HF;
    init_guess = InitGuess::MINAO;
    damping_factor = 0.0;
    tol = 1e-9;
}

RHF::RHF(BasisSet basis_in) : RHF() { basis = basis_in; }

RHF::RHF(BasisSet basis_in, int charge_in) : RHF(basis_in) {
    charge = charge_in;
}

RHF::RHF(BasisSet basis_in, int charge_in, double tol_in, int max_iter_in)
    : RHF(basis_in, charge_in) {
    tol = tol_in;
    max_iter = max_iter_in;
}
void RHF::resolve_occupation(int n_electrons) {
    n_alpha_ = n_electrons / 2;
    n_beta_ = n_electrons / 2;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
