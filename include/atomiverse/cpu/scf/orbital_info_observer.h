/// @file cpu/scf/orbital_info_observer.h
/// @brief OrbitalInfoObserver -- per-step HOMO/LUMO/gap summary from SCFResult.
///
/// Reads eigenvalues from the latest SCFResult published by a sibling
/// `SCFEvaluator` via `SCFResultHandle`. Writes a TSV log and optionally
/// mirrors to stdout.
///
/// Usage:
/// @code
///   cpu::scf::SCFEvaluator scf{config};
///   auto handle = scf.result_handle();                       // share state
///   opt.on(minimizer::Phase::force).add(std::move(scf));
///   opt.on(minimizer::Phase::end_of_step)
///      .add(cpu::scf::OrbitalInfoObserver{handle, "orb.tsv"});
/// @endcode
///
/// Quantities logged (always written in eV):
///   HOMO = highest-occupied orbital energy (alpha channel)
///   LUMO = lowest-unoccupied orbital energy (alpha channel)
///   gap  = LUMO - HOMO
#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Orbital-energy observer component.
///
/// No default bindings: user chooses `end_of_step` of whichever driver hosts
/// the SCFEvaluator (MD, Minimize, or Evaluate).
class OrbitalInfoObserver : public cpu::ComponentBase {
  public:
    /// @brief Construct.
    /// @param handle        Shared view of the SCFEvaluator's latest state.
    /// @param filename      Output TSV path. Empty string disables file output.
    /// @param also_stdout   If true, mirror each row to stdout.
    OrbitalInfoObserver(SCFResultHandle handle, std::string filename = "",
                        bool also_stdout = false);

    std::string name() const override { return "orbital_info_observer"; }

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;
    void destroy(const cpu::DestroyContext& ctx) override;

    /// @brief Number of rows written (excluding header).
    int n_rows() const { return n_rows_; }

    /// @brief HOMO energy [eV] from the last launch (undefined before first launch).
    double homo_eV() const { return homo_eV_; }

    /// @brief LUMO energy [eV] from the last launch.
    double lumo_eV() const { return lumo_eV_; }

    /// @brief HOMO-LUMO gap [eV] from the last launch.
    double gap_eV() const { return lumo_eV_ - homo_eV_; }

  private:
    SCFResultHandle handle_;
    std::string filename_;
    bool also_stdout_;
    std::ofstream os_;
    int n_rows_ = 0;
    double homo_eV_ = 0;
    double lumo_eV_ = 0;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
