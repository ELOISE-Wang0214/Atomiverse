/// @file cpu/kspace/pme.h
/// @brief Smooth Particle-Mesh Ewald reciprocal sum (LAMMPS `kspace_style
/// pppm`/`pme`) — CPU. Requires `ATV_BUILD_KSPACE` (FFTW3).
///
/// A drop-in replacement for `KSpaceEwald`: same real-space `CoulLong`,
/// same self / net-charge terms, but the O(N²) reciprocal sum is replaced
/// by B-spline charge spreading → FFT → influence-function multiply →
/// inverse FFT → force interpolation, giving O(N log N) scaling.
///
/// Smooth PME of Essmann et al., J. Chem. Phys. 103, 8577 (1995):
///   E_recip = kCoulomb/(2πV) · Σ_{m≠0} exp(−π²m²/α²)/m² · B(m) · |F(Q)(m)|²
/// with Q the spline-interpolated charge grid and B(m) the Euler spline
/// modulus.
///
/// **Distributed (MPI) operation.** Each rank spreads its owned charges onto
/// a replicated grid; the grids are `MPI_Allreduce`-summed, every rank runs
/// the (redundant) FFT and interpolates forces on its own atoms. Reciprocal
/// energy / net-charge are added on rank 0; self-energy partitions per rank.
/// (Replicated-grid PME — correct under domain decomposition; a distributed
/// parallel FFT is a future scalability refinement.)
#pragma once

#include <array>
#include <complex>
#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/core/real.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {

/// @brief Smooth Particle-Mesh Ewald reciprocal-space electrostatics.
class KSpacePME : public cpu::ComponentBase {
  public:
    /// @param alpha        Ewald splitting parameter α [1/Å]. Must match the
    ///                     real-space `CoulLong` α.
    /// @param grid_spacing Target real-space grid spacing [Å] (grid sizes are
    ///                     derived from the box and rounded to FFT-friendly
    ///                     values). Default 1.0 Å.
    /// @param order        B-spline interpolation order (≥ 3). Default 6.
    explicit KSpacePME(double alpha, double grid_spacing = 1.0, int order = 6);
    ~KSpacePME() override;

    std::string name() const override { return "pme"; }
    ComponentRole role() const override { return ComponentRole::Evaluator; }
    Capability capabilities() const override {
        return Capability::Energy | Capability::Force | Capability::Stress;
    }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::force), 100}};
    }

    void build(const BuildContext& ctx) override;
    void launch(const LaunchContext& ctx) override;
    void destroy(const DestroyContext& ctx) override;

  private:
    void ensure_grid(const cpu::Frame& frame);

    double alpha_;
    double spacing_;
    int order_;

    int K_[3] = {0, 0, 0}; ///< Grid dimensions.
    std::array<double, 9> cell_cache_{};
    bool have_grid_ = false;

    // Euler spline moduli per grid line, |b_d(m)|².
    std::vector<double> bsp_mod_[3];

    // Charge / convolution grid (complex, size K0·K1·K2) and FFTW plans
    // (opaque to keep fftw3.h out of this header).
    std::vector<std::complex<double>> grid_;
    void* plan_fwd_ = nullptr;
    void* plan_bwd_ = nullptr;
    int planned_n_ = 0;
};

} // namespace cpu
} // namespace ATV
