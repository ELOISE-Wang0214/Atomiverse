/// @file src/cpu/semiempirical/nddo_system.cpp
/// @brief NDDO reference-shaped NDDO system indexing and population helpers.

#include "atomiverse/cpu/semiempirical/nddo_system.h"

#include <cmath>
#include <stdexcept>

#include "atomiverse/cpu/semiempirical/am1_parameters.h"

namespace ATV {
namespace cpu {
namespace semiempirical {

int am1_orbital_count(int z) noexcept {
    const auto* p = am1_parameter(z);
    return p == nullptr ? 0 : p->n_orbitals;
}

double am1_core_charge(int z) noexcept {
    const auto* p = am1_parameter(z);
    return p == nullptr ? 0.0 : p->core_charge;
}

NDDOSystem build_am1_system(const cpu::Frame& frame, int charge,
                            int multiplicity) {
    if (multiplicity < 1) {
        throw std::runtime_error("AM1 multiplicity must be positive");
    }

    NDDOSystem system;
    system.charge = charge;
    system.atoms.reserve(static_cast<std::size_t>(frame.atoms.size()));

    int next_ao = 0;
    double electron_count = -static_cast<double>(charge);
    for (int i = 0; i < frame.atoms.size(); ++i) {
        const int z = frame.atoms.type(i);
        if (!am1_has_parameter(z)) {
            throw std::runtime_error(
                "AM1 parameters unavailable for atomic number "
                + std::to_string(z));
        }
        const int n_orb = am1_orbital_count(z);
        const double core = am1_core_charge(z);
        if (n_orb <= 0 || core == 0.0) {
            throw std::runtime_error(
                "AM1 orbital/core data unavailable for atomic number "
                + std::to_string(z));
        }
        NDDOAtom atom;
        atom.atomic_number = z;
        atom.first_ao = next_ao;
        atom.n_orbitals = n_orb;
        atom.last_ao = next_ao + n_orb - 1;
        atom.core_charge = core;
        system.atoms.push_back(atom);
        next_ao += n_orb;
        electron_count += core;
    }

    system.n_orbitals = next_ao;
    system.n_electrons = static_cast<int>(std::lround(electron_count));
    if (system.n_electrons < 0 || system.n_electrons > 2 * system.n_orbitals) {
        throw std::runtime_error(
            "AM1 electron count is outside the AO capacity");
    }

    const int unpaired = multiplicity - 1;
    if ((system.n_electrons + unpaired) % 2 != 0) {
        throw std::runtime_error(
            "AM1 multiplicity is inconsistent with the electron count");
    }
    system.n_alpha = (system.n_electrons + unpaired) / 2;
    system.n_beta = system.n_electrons - system.n_alpha;
    if (system.n_beta < 0 || system.n_alpha > system.n_orbitals) {
        throw std::runtime_error(
            "AM1 spin electron count is outside AO capacity");
    }
    return system;
}

std::size_t packed_index(int i, int j) noexcept {
    if (i < j)
        return packed_index(j, i);
    return static_cast<std::size_t>(i * (i + 1) / 2 + j);
}

std::vector<double> atomic_populations_from_packed_density(
    const NDDOSystem& system, const std::vector<double>& packed_density) {
    const std::size_t expected =
        static_cast<std::size_t>(system.n_orbitals)
        * static_cast<std::size_t>(system.n_orbitals + 1) / 2u;
    if (packed_density.size() != expected) {
        throw std::runtime_error(
            "packed density size does not match NDDOSystem");
    }

    std::vector<double> populations(system.atoms.size(), 0.0);
    for (std::size_t a = 0; a < system.atoms.size(); ++a) {
        const auto& atom = system.atoms[a];
        double q = 0.0;
        for (int ao = atom.first_ao; ao <= atom.last_ao; ++ao) {
            q += packed_density[packed_index(ao, ao)];
        }
        populations[a] = q;
    }
    return populations;
}

std::vector<double>
atomic_charges_from_packed_density(const NDDOSystem& system,
                                   const std::vector<double>& packed_density) {
    auto populations =
        atomic_populations_from_packed_density(system, packed_density);
    for (std::size_t i = 0; i < populations.size(); ++i) {
        populations[i] = system.atoms[i].core_charge - populations[i];
    }
    return populations;
}

std::vector<double> initial_diagonal_density_guess(const NDDOSystem& system) {
    std::vector<double> diag(static_cast<std::size_t>(system.n_orbitals), 0.0);
    const double yy = static_cast<double>(system.charge)
                      / (static_cast<double>(system.n_orbitals) + 1.0e-10);
    for (const auto& atom : system.atoms) {
        if (atom.n_orbitals == 1) {
            diag[static_cast<std::size_t>(atom.first_ao)] =
                atom.core_charge - yy;
        } else if (atom.n_orbitals == 4) {
            const double w = atom.core_charge * 0.25 - yy;
            for (int ao = atom.first_ao; ao <= atom.last_ao; ++ao) {
                diag[static_cast<std::size_t>(ao)] = w;
            }
        } else if (atom.n_orbitals == 9) {
            const double w = atom.core_charge * 0.25 - yy;
            for (int ao = atom.first_ao; ao < atom.first_ao + 4; ++ao) {
                diag[static_cast<std::size_t>(ao)] = w;
            }
            for (int ao = atom.first_ao + 4; ao <= atom.last_ao; ++ao) {
                diag[static_cast<std::size_t>(ao)] = -yy;
            }
        }
    }
    return diag;
}

} // namespace semiempirical
} // namespace cpu
} // namespace ATV
