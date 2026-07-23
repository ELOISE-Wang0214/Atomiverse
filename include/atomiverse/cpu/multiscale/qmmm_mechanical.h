/// @file cpu/multiscale/qmmm_mechanical.h
/// @brief Mechanical QM/MM: QM energy on tag set + LJ on environment (CPU).
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "atomiverse/core/precision.h"
#include "atomiverse/cpu/component/potential.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"

namespace ATV {
namespace cpu {

/// @brief Mechanical embedding Potential: DFTB2 on H QM tags + LJ elsewhere.
///
/// v1 textbook: QM region is H-only (DFTB2 skeleton); MM is LJ between
/// non-QM pairs and QM–MM pairs (mechanical — no electrostatic embedding).
class QmmmMechanical : public Potential {
  public:
    /// @param qm_tags Global atom tags in the QM region (must be H for DFTB2).
    /// @param epsilon LJ ε [eV] for MM / cross pairs.
    /// @param sigma LJ σ [Å].
    /// @param cutoff LJ cutoff [Å].
    QmmmMechanical(std::vector<TagInt> qm_tags, double epsilon, double sigma,
                   double cutoff);

    std::string name() const override { return "qmmm_mechanical"; }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;

    double last_qm_energy() const { return e_qm_; }
    double last_mm_energy() const { return e_mm_; }

  private:
    std::unordered_set<TagInt> qm_;
    double eps_, sig_, rcut_;
    double e_qm_ = 0, e_mm_ = 0;
    semiempirical::DFTB2 dftb_;
};

} // namespace cpu
} // namespace ATV
