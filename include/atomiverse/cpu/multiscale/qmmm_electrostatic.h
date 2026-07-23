/// @file cpu/multiscale/qmmm_electrostatic.h
/// @brief Electrostatic QM/MM: RHF + MM point charges + LJ (CPU).
#pragma once

#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "atomiverse/core/precision.h"
#include "atomiverse/cpu/component/potential.h"
#include "atomiverse/cpu/multiscale/link_atom.h"
#include "atomiverse/cpu/scf/methods.h"

namespace ATV {
namespace cpu {

/// @brief Electrostatic embedding: QM = RHF/STO-3G on QM tags (+ optional H
///        link atoms); MM atoms as `point_charges` on the QM Hamiltonian;
///        LJ mechanical on MM/cross pairs.
///
/// Link bonds (optional): each `LinkBond{qm_tag, mm_tag}` places a private H
/// cap on the QM geometry with scale `link_scale` [dimensionless]. The MM
/// partner is excluded from electrostatic embedding (charge zeroed at the cut).
/// Link H is never written to the user Frame. Forces on the H (when requested)
/// are projected back onto the real QM host and MM partner via the chain rule.
class QmmmElectrostatic : public Potential {
  public:
    /// @param qm_tags Tags treated as QM.
    /// @param mm_charge Charge [e] assigned to every MM atom (uniform toy).
    /// @param epsilon LJ ε [eV].
    /// @param sigma LJ σ [Å].
    /// @param cutoff LJ cutoff [Å].
    /// @param links Covalent QM–MM cuts (empty = no link atoms).
    /// @param link_scale Dimensionless H placement scale g (default 0.709).
    QmmmElectrostatic(std::vector<TagInt> qm_tags, double mm_charge,
                      double epsilon, double sigma, double cutoff,
                      std::vector<LinkBond> links = {},
                      double link_scale = 0.709);

    std::string name() const override { return "qmmm_electrostatic"; }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;

    double last_qm_energy() const { return e_qm_; }
    double last_mm_energy() const { return e_mm_; }
    /// @brief Number of link hydrogens placed on the last launch.
    int last_n_links() const { return last_n_links_; }

  private:
    std::unordered_set<TagInt> qm_;
    std::vector<LinkBond> links_;
    double mm_q_ = 0;
    double eps_, sig_, rcut_;
    double link_scale_ = 0.709;
    double e_qm_ = 0, e_mm_ = 0;
    int last_n_links_ = 0;
};

} // namespace cpu
} // namespace ATV
