/// @file cpu/vib/hessian.h
/// @brief cpu::vib::Hessian — concrete numerical-Hessian driver.
///
/// `cpu::vib::Hessian` is the user-facing CPU driver of the Response
/// category. It auto-registers the four structural components that
/// assemble a finite-difference (central-difference) molecular Hessian
/// and exposes the result via `run()` / `last_result()`.
///
/// Usage — register any `Capability::Force` evaluator on
/// `vib::Phase::force` (an FF evaluator or `cpu::scf::SCFEvaluator`),
/// exactly as AIMD registers a force evaluator into `cpu::md::MD`:
///
/// @code{.cpp}
///   auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
///   auto cfg  = cpu::scf::rhf_config(ints);
///
///   cpu::vib::Hessian driver;
///   driver.on(ATV::vib::Phase::force)
///         .add(cpu::scf::SCFEvaluator{cfg});
///   driver.build();
///   cpu::vib::HessianResult H = driver.run(frame);   // [3N,3N] eV/Angstrom^2
/// @endcode
///
/// Layering: this header includes no `cpu/scf/` header. The force
/// evaluator is supplied entirely by user code through
/// `on(vib::Phase::force).add(...)`; no domain-to-domain header edge is
/// created.
///
/// Inheritance: Hessian → ResponseDriver → CPUDriver<vib::Phase> →
/// Driver<vib::Phase, cpu::device_tag>.
#pragma once

#include <memory>
#include <string>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/driver/response_driver.h"
#include "atomiverse/cpu/vib/displace_step.h"
#include "atomiverse/cpu/vib/hessian_accumulator.h"
#include "atomiverse/cpu/vib/hessian_result.h"
#include "atomiverse/cpu/vib/hessian_workspace.h"
#include "atomiverse/cpu/vib/reference_geometry.h"
#include "atomiverse/cpu/vib/symmetrize.h"
#include "atomiverse/cpu/vib/thermochem.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace vib {

/// @brief Finite-difference numerical-Hessian driver.
class Hessian : public cpu::ResponseDriver {
  public:
    /// @brief Tag type to construct a driver WITHOUT the HessianAccumulator.
    ///
    /// Exists only so tests can verify `ResponseDriver::validate_components()`
    /// throws when no `Capability::Hessian` component is present. Production
    /// code uses the default constructor.
    struct WithoutAccumulator {};

    /// @brief Output paths for the two auto-injected Recorder components.
    ///
    /// `frequency_tsv` is the `FrequencyReport` TSV path;
    /// `mode_xyz_prefix` is the `NormalModeWriter` file-name prefix
    /// (mode `k` is written to `"<prefix>_mode<k>.xyz"`).
    struct OutputConfig {
        std::string frequency_tsv = "frequencies.tsv"; ///< TSV path.
        std::string mode_xyz_prefix = "mode";          ///< XYZ prefix.
    };

    /// @brief Construct a fully-wired numerical-Hessian driver.
    ///
    /// Auto-registers `ReferenceGeometry`, `DisplaceStep`,
    /// `HessianAccumulator`, `Symmetrize`, `MassWeighting`,
    /// `EigenSolver`, plus the thermochemistry chain (`Thermochem`,
    /// `FrequencyReport`, `NormalModeWriter`), all bound to a
    /// driver-owned shared `HessianWorkspace`. The caller must still
    /// register a `Capability::Force` evaluator on `vib::Phase::force`.
    ///
    /// Thermochemistry uses standard conditions (298.15 K, 1 atm,
    /// closed-shell singlet, symmetry number 1). Use the
    /// `Thermochem::Config` overload to override.
    Hessian();

    /// @brief Construct a numerical-Hessian driver with custom
    ///        thermochemistry inputs.
    ///
    /// Identical wiring to the default constructor, but the auto-injected
    /// `Thermochem` component uses `tc` (temperature / pressure /
    /// multiplicity / symmetry number) instead of standard conditions.
    /// @param tc Thermochemistry inputs forwarded to `cpu::vib::Thermochem`.
    explicit Hessian(const Thermochem::Config& tc);

    /// @brief Construct with custom thermochemistry inputs and custom
    ///        Recorder output paths.
    ///
    /// Identical wiring to the default constructor; the auto-injected
    /// `Thermochem` uses `tc` and the `FrequencyReport` /
    /// `NormalModeWriter` Recorders use the paths in `out`.
    /// @param tc  Thermochemistry inputs forwarded to `cpu::vib::Thermochem`.
    /// @param out Output paths for the FrequencyReport / NormalModeWriter.
    Hessian(const Thermochem::Config& tc, const OutputConfig& out);

    /// @brief Construct a driver WITHOUT the HessianAccumulator (test-only).
    ///
    /// Registers only `ReferenceGeometry`, `DisplaceStep`, and
    /// `Symmetrize`. `build()` then throws from
    /// `ResponseDriver::validate_components()` for the missing
    /// `Capability::Hessian` component.
    /// @param tag Disambiguating tag (`WithoutAccumulator{}`).
    explicit Hessian(WithoutAccumulator tag);

    /// @brief Run the `6N`-displacement finite-difference Hessian assembly.
    ///
    /// Sizes the shared workspace from `frame.atoms.size()`, runs the
    /// Response-category loop (`6N` displacement passes + one finalize
    /// pass), and returns the symmetrized Hessian. The result is also
    /// cached for `last_result()`. On return the frame carries the
    /// equilibrium geometry again.
    ///
    /// @param frame Host-resident system state (positions in Angstrom).
    /// @return `HessianResult` with the `[3N,3N]` Hessian in eV/Angstrom^2.
    HessianResult run(cpu::Frame& frame);

    /// @brief Most recent `HessianResult` produced by `run()`.
    /// @pre `run()` has completed at least once.
    /// @return The cached Hessian result.
    const HessianResult& last_result() const { return last_result_; }

  private:
    /// @brief Shared scratch for the four structural components.
    ///
    /// Heap-allocated so its address is stable across any driver move —
    /// the components hold raw pointers into it.
    std::unique_ptr<HessianWorkspace> workspace_{
        std::make_unique<HessianWorkspace>()};

    /// @brief Cached result of the most recent `run()`.
    HessianResult last_result_;
};

} // namespace vib
} // namespace cpu
} // namespace ATV
