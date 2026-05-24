/// @file cuda/observer/mulliken.cuh
/// @brief MullikenObs -- header-only Mulliken population observer (async D2H bridge).
///
/// Pure observer: computes and stores per-atom charges.  Does NOT write to
/// MolRec or any file format — that is the recorder's job.
///
/// Registers on a `cuda::Evaluate` driver alongside a
/// `cuda::scf::SCFEvaluator`; it runs in the `evaluate::Phase::finalize`
/// pass, after the inner `cuda::scf::SCFSolver` has converged the density.
///
/// The converged density matrix is no longer carried on `cuda::Frame`; the
/// observer holds a reference to the `cuda::scf::SCFEvaluator` and reads the
/// density from its SCF workspace. Construct the evaluator as a stable
/// lvalue and register it borrowed (lvalue `add`), so the reference the
/// observer holds stays valid for the whole run.
///
/// @code{.cpp}
/// cuda::scf::SCFEvaluator scf{integrals, tol, max_iter};
/// MullikenObs mulliken(integrals, scf);
/// cuda::Evaluate eval;
/// eval.add(scf);                          // lvalue -> borrowed
/// eval.on(evaluate::Phase::finalize).add(mulliken);
/// eval.build();
/// eval.run(frame);
/// auto q = mulliken.charges();            // observer: pure data
/// rec.add_field("charge", q);             // recorder: I/O
/// @endcode
#pragma once

#include <memory>
#include <string>
#include <vector>

#include <cuda_runtime.h>

#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/observer/mulliken.h"
#include "atomiverse/cuda/async_bridge.cuh"

#ifdef ATV_HAS_RHF
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cuda/scf/scf_evaluator.cuh"
#endif

namespace ATV {
namespace cuda {

#ifdef ATV_HAS_RHF

/// @brief Mulliken population observer for the CUDA RHF SCF pipeline.
///
/// Header-only — no .cu needed.  Pure observer (stores charges, no I/O).
/// Shares the pure C++ Mulliken implementation with the CPU path; the
/// only CUDA-specific responsibility is bringing the device density
/// matrix back to host via the async D2H bridge.
class MullikenObs : public AsyncBridge {
  public:
    /// @brief Construct a Mulliken observer.
    /// @param ints      Integral pack carrying atoms, shells, and overlap S.
    /// @param evaluator SCF evaluator whose converged workspace supplies the
    ///                  density matrix (borrowed -- must outlive this observer).
    /// @param every     Sampling stride in steps (1 = every step).
    MullikenObs(const cpu::scf::IntegralPack& ints,
                const scf::SCFEvaluator& evaluator, int every = 1)
        : integrals_(ints), evaluator_(&evaluator), every_(every) {}

    ~MullikenObs() override = default;

    std::string name() const override { return "mulliken"; }
    /// @brief Observer — derives Mulliken charges from SCF density (no I/O).
    ComponentRole role() const override { return ComponentRole::Observer; }
    int frequency() const override { return every_; }

    const std::vector<double>& charges() const { return charges_; }
    double energy() const { return energy_val_; }
    int iterations() const { return iterations_; }
    int n_occ() const { return n_occ_; }

  protected:
    void on_build(const BuildContext&) override {
        n_basis_ = evaluator_->workspace().n;
        cudaMallocHost(&h_P_, n_basis_ * n_basis_ * sizeof(double));
        cudaMallocHost(&h_energy_, sizeof(double));
        *h_energy_ = 0.0;
    }

    bool can_run(const LaunchContext&) const override {
        return evaluator_->workspace().P != nullptr;
    }

    void issue_d2h(const LaunchContext& ctx) override {
        auto ds = static_cast<cudaStream_t>(ctx.runtime.defer_stream);
        auto& f = ctx.frame;
        const int n2 = n_basis_ * n_basis_;
        cudaMemcpyAsync(h_P_, evaluator_->workspace().P, n2 * sizeof(double),
                        cudaMemcpyDeviceToHost, ds);
        if (f.scalars)
            cudaMemcpyAsync(h_energy_, &f.scalars->energy, sizeof(double),
                            cudaMemcpyDeviceToHost, ds);
    }

    void on_data_ready(const LaunchContext&) override {
        // Adapt the pinned host buffer as a non-owning xtensor view.
        const std::size_t nb = static_cast<std::size_t>(n_basis_);
        auto P_view = xt::adapt(h_P_, nb * nb, xt::no_ownership(),
                                std::vector<std::size_t>{nb, nb});
        const xt::xtensor<double, 2> P = P_view;

        charges_ = cpu::mulliken_charge(integrals_, P);
        energy_val_ = *h_energy_;
        iterations_ = evaluator_->workspace().iteration;
        n_occ_ = static_cast<int>(evaluator_->workspace().n_occ);
    }

    void on_destroy(const DestroyContext&) override {
        if (h_P_) {
            cudaFreeHost(h_P_);
            h_P_ = nullptr;
        }
        if (h_energy_) {
            cudaFreeHost(h_energy_);
            h_energy_ = nullptr;
        }
    }

  private:
    cpu::scf::IntegralPack integrals_;
    const scf::SCFEvaluator* evaluator_ =
        nullptr; ///< Borrowed SCF evaluator (density source).
    int every_;

    std::vector<double> charges_;
    double energy_val_ = 0;
    int iterations_ = 0;
    int n_occ_ = 0;

    int n_basis_ = 0;
    double* h_P_ = nullptr;
    double* h_energy_ = nullptr;
};

#endif // ATV_HAS_RHF

} // namespace cuda
} // namespace ATV
