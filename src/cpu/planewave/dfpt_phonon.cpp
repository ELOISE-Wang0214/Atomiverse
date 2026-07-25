/// @file cpu/planewave/dfpt_phonon.cpp
/// @brief DfptPhonon — mass weighting and diagonalization over the internal
///        Γ-point DFPT force-constant assembly.

#include "atomiverse/cpu/planewave/dfpt_phonon.h"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/planewave/internal/dfpt.h"
#include "atomiverse/cpu/planewave/internal/pw_subspace.h"
#include "atomiverse/cpu/vibration/atomic_masses.h"
#include "atomiverse/cpu/vibration/nist_constants.h"

namespace ATV {
namespace cpu {
namespace pw {

DfptPhononResult DfptPhonon::compute(const PW& method,
                                     const cpu::Frame& frame) const {
    const auto out = internal::dfpt_gamma_force_constants(
        frame, method, AccReal(response_tol), response_max_iter,
        AccReal(cg_tol), cg_max_iter, AccReal(mix_beta));

    const int n = frame.nlocal();
    const int dim = 3 * n;

    DfptPhononResult result;
    result.n_atoms = n;
    result.scf_converged = out.scf_converged;
    result.response_converged = out.response_converged;

    // Force constants: Ha/Bohr² → eV/Å² for the public field (vib units).
    const double to_ev_ang2 =
        1.0 / vib::nist::EV_PER_ANG2_TO_HARTREE_PER_BOHR2;
    result.max_asymmetry =
        static_cast<double>(out.max_asymmetry) * to_ev_ang2;
    result.force_constants = xt::zeros<double>(
        {static_cast<std::size_t>(dim), static_cast<std::size_t>(dim)});
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            result.force_constants(i, j) =
                static_cast<double>(
                    out.force_constants[static_cast<std::size_t>(i) * dim
                                        + j])
                * to_ev_ang2;

    // Mass-weighted dynamical matrix in atomic units [Ha/(Bohr² amu)],
    // the same convention as cpu::vib::MassWeighting.
    const int* types = frame.span(kAtom_type);
    std::vector<double> inv_sqrt_m(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        inv_sqrt_m[static_cast<std::size_t>(i)] =
            1.0 / std::sqrt(vib::mass_amu(types[i]));

    std::vector<AccReal> dmat(static_cast<std::size_t>(dim) * dim);
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            dmat[static_cast<std::size_t>(i) * dim + j] =
                out.force_constants[static_cast<std::size_t>(i) * dim + j]
                * static_cast<AccReal>(inv_sqrt_m[static_cast<std::size_t>(
                      i / 3)])
                * static_cast<AccReal>(inv_sqrt_m[static_cast<std::size_t>(
                      j / 3)]);

    std::vector<AccReal> eig(static_cast<std::size_t>(dim));
    std::vector<AccReal> vecs(static_cast<std::size_t>(dim) * dim);
    internal::jacobi_eigh(dim, dmat.data(), eig.data(), vecs.data(),
                          /*max_sweep=*/50);

    const double au2cm = vib::nist::au2wavenumber();
    result.frequencies = xt::zeros<double>({static_cast<std::size_t>(dim)});
    for (int k = 0; k < dim; ++k) {
        const double lambda = static_cast<double>(eig[static_cast<
            std::size_t>(k)]);
        result.frequencies(k) = (lambda < 0.0)
                                    ? -std::sqrt(-lambda) * au2cm
                                    : std::sqrt(lambda) * au2cm;
    }
    return result;
}

} // namespace pw
} // namespace cpu
} // namespace ATV
