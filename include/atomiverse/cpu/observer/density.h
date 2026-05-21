/// @file cpu/observer/density.h
/// @brief Electron density grid evaluation on the CPU (STO-3G basis).
///
/// Called by cuda::DensityRecorder after async D2H of the density matrix P.
/// The CUDA component owns no computation -- it lives here.
///
/// **Phase 1 note:** unlike `cpu::observer::RDF`, this is **not** wrapped as
/// a `cpu::AsyncComponentBase` subclass. The evaluation requires a
/// density matrix P from SCF, which is not a `cpu::Frame` field and is not
/// carried in the default `FrameSnapshot`. The natural async-observer
/// integration point is an SCF-driver observer wrapper that uses
/// `FrameSnapshot::custom_payload` to shuttle P from the SCF solver to the
/// analysis thread. That hook depends on `cpu::scf::SCFEvaluator` being
/// extended to publish P into custom_payload, plus a concrete SCF driver
/// class to auto-add the matching `SnapshotPublisher`. Deferred to Phase 1b
/// or to the AIMD density-observation workstream.
///
/// @see cuda::DensityRecorder (cuda/recorder/density_recorder.cuh)
/// @see .claude/specs/parallel-md-engine.md  § B.1 custom_payload
#pragma once

#include <vector>

namespace ATV {
namespace cpu {

// -- Result type ---------------------------------------------------------------

/// @brief Uniform real-space grid of electron density rho(r).
///
/// The grid is axis-aligned and non-periodic.
/// Indexing: rho[ix * ny * nz + iy * nz + iz]  (C order).
struct DensityGrid {
    int nx = 0, ny = 0, nz = 0;    ///< Grid dimensions.
    double ox = 0, oy = 0, oz = 0; ///< Grid origin [Angstrom].
    double lx = 0, ly = 0, lz = 0; ///< Box extent [Angstrom].
    std::vector<double> rho;       ///< Density values [nx*ny*nz, e/Angstrom^3].

    bool empty() const { return rho.empty(); }
};

// -- Free function -------------------------------------------------------------

/// @brief Evaluate RHF electron density on a uniform grid (STO-3G, H/He/Li only).
///
/// rho(r) = sum_{mu,nu} P_{mu,nu} phi_mu(r) phi_nu(r)
///
/// The grid bounds are determined from atom positions +/- padding.
/// Returns an empty DensityGrid (rho.empty() == true) if any atom's Z
/// is outside the supported range [1, 3].
///
/// @param type_id   Atomic numbers [n_atoms].
/// @param x         Atom x positions [n_atoms, Angstrom].
/// @param y         Atom y positions [n_atoms, Angstrom].
/// @param z         Atom z positions [n_atoms, Angstrom].
/// @param P         Density matrix [n_basis * n_basis], row-major (P is symmetric).
/// @param n_basis   Number of STO-3G basis functions.
/// @param nx        Grid points along x.
/// @param ny        Grid points along y.
/// @param nz        Grid points along z.
/// @param padding   Extra padding around atom bounds [Angstrom].
/// @return Populated DensityGrid, or empty on failure.
DensityGrid eval_density_grid(const std::vector<int>& type_id,
                              const std::vector<double>& x,
                              const std::vector<double>& y,
                              const std::vector<double>& z, const double* P,
                              int n_basis, int nx, int ny, int nz,
                              double padding);

} // namespace cpu
} // namespace ATV
