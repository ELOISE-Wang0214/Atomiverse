/// @file cpu/multiscale/qmmm_electrostatic.h
/// @brief Electrostatic QM/MM: RHF + MM point charges + LJ (CPU).
#pragma once

#include <string>
#include <unordered_set>
#include <vector>

#include "atomiverse/core/precision.h"
#include "atomiverse/cpu/component/potential.h"
#include "atomiverse/cpu/scf/methods.h"

namespace ATV {
namespace cpu {

/// @brief Electrostatic embedding: QM = RHF/STO-3G on QM tags; MM atoms as
///        `point_charges` on the QM Hamiltonian; LJ mechanical on MM/cross.
class QmmmElectrostatic : public Potential {
  public:
    /// @param qm_tags Tags treated as QM.
    /// @param mm_charge Charge [e] assigned to every MM atom (uniform toy).
    /// @param epsilon LJ ε [eV].
    /// @param sigma LJ σ [Å].
    /// @param cutoff LJ cutoff [Å].
    QmmmElectrostatic(std::vector<TagInt> qm_tags, double mm_charge,
                      double epsilon, double sigma, double cutoff);

    std::string name() const override { return "qmmm_electrostatic"; }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;

    double last_qm_energy() const { return e_qm_; }
    double last_mm_energy() const { return e_mm_; }

  private:
    std::unordered_set<TagInt> qm_;
    double mm_q_ = 0;
    double eps_, sig_, rcut_;
    double e_qm_ = 0, e_mm_ = 0;
};

} // namespace cpu
} // namespace ATV
