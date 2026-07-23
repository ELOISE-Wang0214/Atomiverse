/// @file cpu/multiscale/qmmm_polarizable.h
/// @brief Polarizable QM/MM embedding: RHF + induced MM dipoles (CPU).
///
/// MM sites carry isotropic polarizability \(\alpha\) [Å³]. Induced dipoles
/// obey the Applequist/Thole self-consistent field
/// \(\boldsymbol{\mu}_i=\alpha_i\mathbf{E}_i^{\mathrm{tot}}\) with permanent
/// field from QM nuclei (Z as point charges) and permanent MM charges, and
/// mutual induction via a Thole-damped dipole tensor. Polarization energy
/// \(E_{\mathrm{pol}}=-\frac12\sum_i\boldsymbol{\mu}_i\cdot\mathbf{E}_i^{\mathrm{perm}}\)
/// [eV] is added to the RHF electrostatic embedding energy + LJ.
///
/// Thompson & Schenter, *J. Phys. Chem.* **99**, 6374 (1995);
/// Thole, *Chem. Phys.* **59**, 341 (1981). Units: e, Å, eV.
#pragma once

#include <string>
#include <unordered_set>
#include <vector>

#include "atomiverse/core/precision.h"
#include "atomiverse/cpu/component/potential.h"

namespace ATV {
namespace cpu {

/// @brief One iteration of induced-dipole SCC (Gaussian electrostatics).
///
/// Permanent field \(\mathbf{E}^{\mathrm{perm}}\) [e/Å²] and mutual \(\mathbf{T}\)
/// use Gaussian units; \(\boldsymbol{\mu}=\alpha\mathbf{E}\) with \(\alpha\) in
/// Å³ yields \(\boldsymbol{\mu}\) in e·Å. Energy conversion uses
/// `units::coulomb` so \(E_{\mathrm{pol}}\) is in eV.
///
/// @param n Number of polarizable sites.
/// @param pos Site positions, length 3n, packed xyz [Å].
/// @param alpha Isotropic polarizability per site [Å³].
/// @param e_perm Permanent field per site, length 3n [e/Å²].
/// @param thole_a Thole damping factor (dimensionless; 0 disables damping).
/// @param tol Max |\Deltaμ| convergence threshold [e·Å].
/// @param max_iter Maximum SCC iterations.
/// @param[out] mu Induced dipoles, length 3n [e·Å] (resized if needed).
/// @param[out] iters_used Optional; written with SCC iterations used.
/// @return Polarization energy [eV].
/// @throws std::runtime_error if SCC does not converge within max_iter.
double solve_induced_dipoles(int n, const double* pos, const double* alpha,
                             const double* e_perm, double thole_a, double tol,
                             int max_iter, std::vector<double>& mu,
                             int* iters_used = nullptr);

/// @brief Electrostatic QM/MM with isotropic induced MM dipoles.
///
/// QM region: RHF/STO-3G with permanent MM point charges as embedding.
/// MM region: uniform permanent charge `mm_charge` [e], isotropic α [Å³],
/// Thole-damped mutual induction, LJ mechanical on non-QM–QM pairs.
class QmmmPolarizable : public Potential {
  public:
    /// @param qm_tags Tags treated as QM.
    /// @param mm_charge Permanent charge [e] on every MM atom (toy uniform).
    /// @param alpha Isotropic polarizability [Å³] on every MM atom.
    /// @param thole_a Thole damping factor (dimensionless; default 0.39).
    /// @param epsilon LJ ε [eV].
    /// @param sigma LJ σ [Å].
    /// @param cutoff LJ cutoff [Å].
    /// @param pol_tol Induced-dipole SCC threshold [e·Å] (default 1e-8).
    /// @param pol_max_iter Max SCC iterations (default 50).
    QmmmPolarizable(std::vector<TagInt> qm_tags, double mm_charge, double alpha,
                    double thole_a, double epsilon, double sigma, double cutoff,
                    double pol_tol = 1e-8, int pol_max_iter = 50);

    std::string name() const override { return "qmmm_polarizable"; }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;

    double last_qm_energy() const { return e_qm_; }
    double last_mm_energy() const { return e_mm_; }
    /// @brief Polarization energy of the last launch [eV] (≤ 0 for α ≥ 0).
    double last_pol_energy() const { return e_pol_; }
    /// @brief SCC iterations used on the last launch.
    int last_pol_iters() const { return last_pol_iters_; }
    /// @brief Whether induced-dipole SCC converged on the last launch.
    bool last_pol_converged() const { return last_pol_converged_; }
    /// @brief Induced dipoles of the last launch (packed xyz, e·Å).
    const std::vector<double>& last_mu() const { return last_mu_; }

  private:
    std::unordered_set<TagInt> qm_;
    double mm_q_ = 0;
    double alpha_ = 0;
    double thole_a_ = 0.39;
    double eps_, sig_, rcut_;
    double pol_tol_ = 1e-8;
    int pol_max_iter_ = 50;
    double e_qm_ = 0, e_mm_ = 0, e_pol_ = 0;
    int last_pol_iters_ = 0;
    bool last_pol_converged_ = false;
    std::vector<double> last_mu_;
};

} // namespace cpu
} // namespace ATV
