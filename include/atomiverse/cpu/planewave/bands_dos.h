/// @file cpu/planewave/bands_dos.h
/// @brief Band structure and density-of-states post-processing (CPU PW).
///
/// After a self-consistent density is obtained from `cpu::pw::PW` members,
/// non-self-consistent diagonalization on a frozen \(V_\mathrm{eff}\) yields
/// \(\varepsilon_n(\mathbf{k})\) (bands) and a Gaussian-broadened DOS.
/// Parameters are members of `Bands` / `Dos` — no public `*Config` bag.
///
/// @code
///   cpu::pw::PW pw;
///   pw.upf_path = "Si.upf";
///   pw.ecut_ry = 15.0;
///   pw.k_mesh = 2;
///   pw.xc = cpu::scf::XCType::LSDA;
///
///   cpu::pw::Bands bands;
///   bands.n_bands = 8;
///   bands.set_diamond_path(4);
///   auto br = bands.compute(pw, frame);   // gap_ry, eigenvalues_ry
///
///   cpu::pw::Dos dos;
///   dos.k_mesh = 4;
///   dos.n_bands = 8;
///   dos.sigma_ry = 0.05;
///   auto dr = dos.compute(pw, frame);     // dos, n_elec_to_vbm
/// @endcode
#pragma once

#include <array>
#include <vector>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/planewave/pw_solver.h"

namespace ATV {
namespace cpu {
namespace pw {

/// @brief Result of a band-structure (NSCF path) calculation.
struct BandsResult {
    /// Fractional k-points used for the path (same order as eigenvalues).
    std::vector<std::array<double, 3>> kfrac;
    /// Eigenvalues \(\varepsilon_n(\mathbf{k})\) [Ry], layout `[ik][ib]`.
    std::vector<std::vector<double>> eigenvalues_ry;
    double e_vbm_ry = 0; ///< Valence-band maximum [Ry].
    double e_cbm_ry = 0; ///< Conduction-band minimum [Ry].
    double gap_ry = 0;   ///< \(E_\mathrm{cbm}-E_\mathrm{vbm}\) [Ry].
    int vbm_k = 0;       ///< k-index of VBM.
    int cbm_k = 0;       ///< k-index of CBM.
    int vbm_band = 0;    ///< Band index of VBM (0-based).
    int cbm_band = 0;    ///< Band index of CBM (0-based).
    bool scf_converged = false; ///< Whether the preceding SCF converged.
    double energy_ry = 0;       ///< SCF total energy [Ry].
};

/// @brief Result of a Gaussian-broadened density of states.
struct DosResult {
    std::vector<double> energy_ry; ///< Bin centres [Ry].
    /// DOS \(D(E)\) [states / Ry / cell]. Closed-shell: each band carries
    /// weight 2 when integrating the occupied spectrum.
    std::vector<double> dos;
    /// \(\int_{-\infty}^{E_\mathrm{vbm}} D(E)\,dE\) (target \(\approx n_e\)).
    double n_elec_to_vbm = 0;
    double e_vbm_ry = 0;
    double e_cbm_ry = 0;
    double gap_ry = 0;
    bool scf_converged = false;
};

/// @brief High-symmetry-path band structure (NSCF on frozen \(V_\mathrm{eff}\)).
///
/// SCF uses `method` members (`ecut_ry`, `k_mesh`, `xc`, …). Path k-points
/// are `kfrac` (filled by `set_diamond_path` when empty at `compute` time).
class Bands {
  public:
    /// @brief Total bands to diagonalize (valence + conduction).
    ///        Silicon has 4 occupied valence bands (\(n_e=8\)).
    int n_bands = 8;
    /// @brief Fractional (crystal) k-list. Empty → `set_diamond_path(4)`.
    std::vector<std::array<double, 3>> kfrac;

    /// @brief FCC/diamond high-symmetry path Γ–X–W–L–Γ–K.
    /// @param n_seg Number of samples per segment (including endpoints once).
    void set_diamond_path(int n_seg = 4);

    /// @brief SCF + NSCF bands for the geometry in `frame`.
    /// @param method Plane-wave method (members only; not mutated).
    /// @param frame  Atom positions [Å], types (Z), simulation cell.
    /// @return Owned `BandsResult` with eigenvalues [Ry] and gap.
    BandsResult compute(const PW& method, const cpu::Frame& frame) const;
};

/// @brief Gaussian-broadened DOS from denser Monkhorst–Pack eigenvalues.
class Dos {
  public:
    /// @brief DOS k-mesh per dimension (MP \(k_\mathrm{mesh}^3\)).
    int k_mesh = 4;
    /// @brief Bands per k-point (must cover valence + some conduction).
    int n_bands = 8;
    /// @brief Gaussian broadening width [Ry].
    double sigma_ry = 0.05;
    /// @brief Number of energy bins for \(D(E)\).
    int n_bins = 200;
    /// @brief Energy window [Ry]; if `e_max_ry <= e_min_ry`, auto from evals.
    double e_min_ry = 0;
    double e_max_ry = 0;

    /// @brief SCF (method.k_mesh) + denser-mesh NSCF + Gaussian DOS.
    /// @param method Plane-wave method (members; SCF uses method.k_mesh).
    /// @param frame  Geometry [Å].
    /// @return Owned `DosResult`.
    DosResult compute(const PW& method, const cpu::Frame& frame) const;
};

/// @brief Locate VBM/CBM from a spectrum with `n_occ` occupied bands.
///
/// Pure helper for unit tests and both post-processors. Eigenvalues must be
/// sorted ascending per k-point. Units are those of the input (Ry or Ha).
///
/// @param evals  `[ik][ib]` eigenvalues (sorted ascending in `ib`).
/// @param n_occ  Number of occupied bands (\(n_e/2\) for closed shell).
/// @param e_vbm  Output valence-band maximum.
/// @param e_cbm  Output conduction-band minimum.
/// @param vbm_k  Output k-index of VBM.
/// @param cbm_k  Output k-index of CBM.
/// @param vbm_b  Output band index of VBM.
/// @param cbm_b  Output band index of CBM.
/// @return gap = e_cbm − e_vbm.
double gap_from_spectrum(const std::vector<std::vector<double>>& evals,
                         int n_occ, double& e_vbm, double& e_cbm, int& vbm_k,
                         int& cbm_k, int& vbm_b, int& cbm_b);

/// @brief Gaussian DOS on a uniform energy grid (unit-testable free math).
///
/// \(D(E_i)=\sum_{n k} w_{nk}\, (2/\sqrt{2\pi}\sigma)\,
/// \exp(-(E_i-\varepsilon)^2/(2\sigma^2))\). The factor 2 is spinless
/// closed-shell band degeneracy.
///
/// @param evals     Flat list of eigenvalues [energy units].
/// @param weights   Per-eigenvalue k-weight (same length as evals).
/// @param e_min     Grid start.
/// @param e_max     Grid end.
/// @param n_bins    Number of bins.
/// @param sigma     Gaussian width (same energy units).
/// @param energy    Output bin centres.
/// @param dos       Output \(D(E)\).
void gaussian_dos(const std::vector<double>& evals,
                  const std::vector<double>& weights, double e_min,
                  double e_max, int n_bins, double sigma,
                  std::vector<double>& energy, std::vector<double>& dos);

} // namespace pw
} // namespace cpu
} // namespace ATV
