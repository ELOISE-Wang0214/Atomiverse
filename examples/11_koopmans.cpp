/// @file examples/11_koopmans.cpp
/// @brief Koopmans' theorem: vertical ionization potential from the RHF HOMO.
///
/// Koopmans' theorem (T. Koopmans, Physica 1, 104 (1934)) states that, in
/// the frozen-orbital approximation, the vertical ionization potential of a
/// closed-shell molecule is the negative of its highest-occupied canonical
/// RHF orbital energy: IP ≈ −ε_HOMO.
///
/// This example sweeps five molecules × two basis sets, runs a single-point
/// RHF via `cpu::Evaluate + cpu::scf::SCF`, reads the HOMO through
/// `cpu::scf::OrbitalInfoObserver` (bound at `evaluate::Phase::finalize`,
/// sharing state via `scf.result_handle()`), and tabulates:
///
///   system | basis | −ε_HOMO [Ha] | IP_Koopmans [eV] | experimental IP [eV]
///
/// Units: orbital energies in Hartree; IP in eV via 1 Ha = 27.211386245988 eV
/// (`ATV::units::kHartreeInEv`). The experimental column is a *qualitative* reference
/// only — Koopmans neglects orbital relaxation and correlation, so no numeric
/// tolerance is implied. The reference −ε_HOMO values here are regression-
/// pinned in benchmarks/qc/rhf/validation/test_koopmans_parity.cpp.
///
/// Build:
///   cmake -B build -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
///   cmake --build build --target 11_koopmans
///
/// Run:
///   ./build/examples/11_koopmans

#ifdef ATV_HAS_RHF

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/driver/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/orbital_info_observer.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_types.h"

using namespace ATV;
using ATV::cpu::scf::BasisSet;

namespace {

struct Atom {
    int z;           ///< Atomic number.
    double x, y, zc; ///< Position [Angstrom].
};

struct Molecule {
    std::string name;        ///< Label.
    std::vector<Atom> atoms; ///< Geometry [Å], type = Z.
    double exp_ip;           ///< Experimental IP [eV] (qualitative).
    std::string exp_label;   ///< "adiabatic" / "vertical" tag.
};

// Reference geometries [Å]; experimental IPs from NIST (qualitative only).
const std::vector<Molecule> kMolecules = {
    {"He", {{2, 0, 0, 0}}, 24.587, "vertical"},
    {"H2", {{1, 0, 0, 0}, {1, 0, 0, 0.7414}}, 15.426, "adiabatic"},
    {"H2O",
     {{8, 0, 0, 0.1173}, {1, 0, 0.7572, -0.4692}, {1, 0, -0.7572, -0.4692}},
     12.62,
     "vertical"},
    {"NH3",
     {{7, 0, 0, 0},
      {1, 0, -0.9377, -0.3816},
      {1, 0.8121, 0.4689, -0.3816},
      {1, -0.8121, 0.4689, -0.3816}},
     10.8,
     "vertical"},
    {"CH4",
     {{6, 0, 0, 0},
      {1, 0.627, 0.627, 0.627},
      {1, 0.627, -0.627, -0.627},
      {1, -0.627, 0.627, -0.627},
      {1, -0.627, -0.627, 0.627}},
     13.6,
     "vertical (Jahn-Teller)"},
};

cpu::Frame make_frame(const Molecule& mol) {
    cpu::Frame f(static_cast<int>(mol.atoms.size()));
    for (std::size_t i = 0; i < mol.atoms.size(); ++i) {
        f.span(ATV::kAtom_type)[i] = mol.atoms[i].z;
        f.span(ATV::kAtom_x)[i] = mol.atoms[i].x;
        f.span(ATV::kAtom_y)[i] = mol.atoms[i].y;
        f.span(ATV::kAtom_z)[i] = mol.atoms[i].zc;
    }
    f.simbox = ATV::cpu::Box::cubic(100.0);
    return f;
}

} // namespace

int main() {
    std::cout << "=== Koopmans' theorem: IP ≈ −ε_HOMO (RHF, CPU) ===\n\n"
              << std::left << std::setw(6) << "system" << std::setw(9)
              << "basis" << std::right << std::setw(15) << "−ε_HOMO [Ha]"
              << std::setw(18) << "IP_Koopmans [eV]" << "   " << std::left
              << "experimental IP [eV]\n"
              << std::string(78, '-') << "\n";

    const BasisSet bases[] = {BasisSet::STO3G, BasisSet::Pople631G};

    for (const auto& mol : kMolecules) {
        for (BasisSet basis : bases) {
            cpu::Frame frame = make_frame(mol);

            // Size an SCF config from a representative integral pack, then let
            // the SCF rebuild integrals from the frame each launch.
            auto integrals = cpu::scf::build_integrals(frame, basis);
            auto config = cpu::scf::rhf_config(integrals, /*tol=*/1e-10);

            cpu::scf::SCF scf{config, basis};
            scf.with_gradient(false); // energy + orbitals only (single point)
            auto handle = scf.result_handle();

            // Evaluate runs the SCF once; OrbitalInfoObserver reads the HOMO
            // from the shared handle at finalize (and can mirror to stdout /
            // a TSV — here we keep it quiet and read the handle ourselves).
            cpu::Evaluate eval;
            eval.on(evaluate::Phase::evaluate).add(std::move(scf));
            eval.on(evaluate::Phase::finalize)
                .add(cpu::scf::OrbitalInfoObserver{handle});
            eval.build();
            eval.run(frame);

            // The handle carries exactly what OrbitalInfoObserver::homo_eV()
            // reads: ε_HOMO = eigenvalues_alpha[n_alpha − 1] [Hartree].
            const auto& res = *handle.result;
            const double eps_homo_ha = res.eigenvalues_alpha(res.n_alpha - 1);
            const double neg_eps_homo_ha = -eps_homo_ha;
            const double ip_koopmans_eV = neg_eps_homo_ha * units::kHartreeInEv;

            std::cout << std::left << std::setw(6) << mol.name << std::setw(9)
                      << cpu::scf::basis_set_name(basis) << std::right
                      << std::fixed << std::setprecision(8) << std::setw(15)
                      << neg_eps_homo_ha << std::setprecision(3)
                      << std::setw(18) << ip_koopmans_eV << "   " << std::left
                      << std::setprecision(3) << mol.exp_ip << " ("
                      << mol.exp_label << ")\n";
        }
    }

    std::cout << "\nNote: He/STO-3G has a single basis function and no virtual "
                 "orbital,\n"
                 "so only the HOMO / IP are meaningful (no HOMO-LUMO gap).\n"
                 "Koopmans neglects orbital relaxation and correlation — the "
                 "experimental\ncolumn is a qualitative guide, not a "
                 "numerical target.\n";
    return 0;
}

#else // ATV_HAS_RHF

#include <iostream>

int main() {
    std::cout << "11_koopmans requires -DATV_BUILD_RHF=ON\n";
    return 0;
}

#endif // ATV_HAS_RHF
