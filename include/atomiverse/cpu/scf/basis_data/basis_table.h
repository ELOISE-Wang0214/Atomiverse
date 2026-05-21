/// @file cpu/scf/basis_data/basis_table.h
/// @brief Generic container types for checked-in basis-set data.
///
/// Each basis-set header under `basis_data/` declares one inline
/// function returning a `const std::vector<AtomBasis>&` indexed by
/// atomic number Z (position 0 is unused; empty AtomBasis marks an
/// unsupported element).
///
/// The types model the same "raw" / "unnormalized" contraction
/// convention used by PySCF's `mol._basis`: `coeff[p]` is the Pople
/// contraction coefficient for the primitive with exponent `alpha[p]`,
/// and the physical GTO normalization (`phys_norm_s(alpha)` etc.) is
/// applied by the shell-list builder.
///
/// @see cpu::scf::build_shell_list for the builder that consumes these
///      tables.
/// @see tools/gen_basis_data.py for the generator that emits them.
#pragma once

#include <vector>

namespace ATV {
namespace cpu {
namespace scf {
namespace basis_data {

/// @brief One contracted Gaussian shell.
///
/// For `l = 1` (p), the same `alpha/coeff` pair is reused across the
/// three Cartesian components px/py/pz. For `l = 2` (d), likewise
/// across the six Cartesian components d_xx..d_yz (Stage B2b).
struct ContractedShell {
    int l;                     ///< Angular momentum (0 = s, 1 = p, 2 = d).
    std::vector<double> alpha; ///< Gaussian exponents, length `n_prim`.
    std::vector<double>
        coeff; ///< Contraction coefficients (unnormalized), length `n_prim`.
};

/// @brief All shells belonging to one element.
///
/// `shells` is empty for atomic numbers that the basis set does not
/// cover. The shell-list builder throws `std::runtime_error` when it
/// encounters such an atom.
struct AtomBasis {
    std::vector<ContractedShell> shells;
};

} // namespace basis_data
} // namespace scf
} // namespace cpu
} // namespace ATV
