/// @file cpu/gpw/internal/multigrid.h
/// @brief CP2K-style CUTOFF/REL_CUTOFF multigrid for GPW collocation.
///
/// A single uniform grid forces every Gaussian onto the finest spacing: diffuse
/// Gaussians then span huge bounding boxes (many points) — slow — and a single
/// grid leaves an O(mHa) discretization residual vs CP2K's multigrid. CP2K maps
/// each Gaussian (by its exponent) to the COARSEST grid level that still
/// resolves it (REL_CUTOFF), collocates there (compact support, far fewer points
/// for diffuse functions), interpolates the per-level densities to a common fine
/// grid for the (nonlinear) Hartree + XC, then restricts the fine potential back
/// to each level to integrate the KS matrix.
///
/// Transfer between levels is real-space trilinear interpolation (prolong) and
/// its exact transpose (restrict). The two being adjoints (R = c·Pᵀ) keeps the
/// integrated KS matrix the exact gradient of the collocated energy, so the SCF
/// stays variational/convergent; the discretization is faithful (energy → the
/// fine-grid limit as the grids refine), not an uncontrolled approximation.
///
/// Gated: `ngrids == 1` bypasses all of this and uses the single-grid path
/// unchanged (bit-identical).
#pragma once

#include <vector>

#include "atomiverse/cpu/gpw/internal/collocate.h"

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief A multigrid: `levels[0]` is the finest (at CUTOFF); coarser levels
/// follow at geometrically smaller cutoffs. The finest grid is the reference on
/// which Poisson + XC run.
struct MultiGrid {
    std::vector<UniformGrid> levels;             ///< finest first.
    std::vector<double> level_cutoff_ry;         ///< per-level density cutoff [Ry].
};

/// @brief Build `ngrids` geometrically-spaced grid levels for an orthorhombic
///        cell. Finest cutoff = `cutoff_ry`; each coarser level divides the
///        cutoff by a fixed ratio so the whole REL_CUTOFF..CUTOFF range is
///        covered. `ngrids == 1` returns the single CUTOFF grid.
/// @param cell_abc Box lengths a,b,c [Bohr].
/// @param cutoff_ry Finest (plane-wave) density cutoff [Ry].
/// @param ngrids Number of grid levels (>= 1).
/// @param rel_cutoff Reference cutoff [Ry] for the exponent->level mapping; sets
///        the coarsest level's cutoff ~ rel_cutoff and the geometric ratio.
MultiGrid make_multigrid(const double cell_abc[3], double cutoff_ry, int ngrids,
                         double rel_cutoff);

/// @brief Assign each basis function to a grid level by its most-diffuse
///        exponent. A Gaussian of exponent α needs an effective density cutoff
///        ~ rel_cutoff·(α/α_ref); it goes on the coarsest level whose cutoff
///        still meets that need (tight α -> fine level 0; diffuse α -> coarse).
/// @return level index per basis function, length basis.size().
std::vector<int> assign_levels(const MultiGrid& mg,
                               const std::vector<GridBasisFunction>& basis,
                               double rel_cutoff);

/// @brief Prolong (interpolate) a coarse-grid field onto the fine grid by
///        periodic trilinear interpolation. @p periodic wraps indices.
std::vector<double> prolong_trilinear(const UniformGrid& coarse,
                                      const UniformGrid& fine,
                                      const std::vector<double>& f_coarse,
                                      bool periodic);

/// @brief Restrict a fine-grid field onto the coarse grid as the EXACT transpose
///        of `prolong_trilinear` (the same trilinear weights, scattered) times
///        the cell-volume-per-point ratio so that ⟨V_fine, P f_coarse⟩_fine =
///        ⟨R V_fine, f_coarse⟩_coarse — i.e. integration is the adjoint of
///        collocation across the level transfer.
std::vector<double> restrict_transpose(const UniformGrid& fine,
                                       const UniformGrid& coarse,
                                       const std::vector<double>& f_fine,
                                       bool periodic);

} // namespace gpw
} // namespace cpu
} // namespace ATV
