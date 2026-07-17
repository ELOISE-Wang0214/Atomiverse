/// @file examples/11_koopmans.cpp
/// @brief Koopmans IP from RHF HOMO (concrete method API).

#include <iomanip>
#include <iostream>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;

static cpu::Frame make_water() {
    cpu::Frame f(3);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_x)[1] = 0.7572;
    f.span(kAtom_x)[2] = -0.7572;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_y)[2] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_z)[1] = -0.5860;
    f.span(kAtom_z)[2] = -0.5860;
    f.span(kAtom_type)[0] = 8;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_mass)[0] = 15.999;
    f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_mass)[2] = 1.008;
    f.simbox = cpu::Box::cubic(20.0);
    return f;
}

int main() {
    auto frame = make_water();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-10;
    const auto& r = rhf.run(frame);
    // HOMO = last occupied (n_alpha - 1)
    const int n_occ = r.n_alpha;
    const double homo_ha = r.eigenvalues_alpha(n_occ - 1);
    const double ip_ev = -homo_ha * units::kHartreeInEv;
    std::cout << std::setprecision(6) << "H2O RHF HOMO = " << homo_ha << " Ha\n"
              << "Koopmans IP  = " << ip_ev << " eV\n";
    return r.converged ? 0 : 1;
}
