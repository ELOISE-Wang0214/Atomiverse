/// @file cuda/observer/mulliken.cuh
/// @brief MullikenObs -- header-only Mulliken population observer (async D2H bridge).
///
/// Pure observer: computes and stores per-atom charges.  Does NOT write to
/// MolRec or any file format — that is the recorder's job.
///
/// @code{.cpp}
/// MullikenObs mulliken(integrals);
/// solver.on(scf::Phase::end_of_step).add(mulliken);
/// solver.run(frame, 100);
/// auto q = mulliken.charges();            // observer: pure data
/// rec.add_field("charge", q);             // recorder: I/O
/// @endcode
#pragma once

#include <memory>
#include <string>
#include <vector>

#include <cuda_runtime.h>

#include <xtensor/containers/xtensor.hpp>
#include <xtensor/containers/xadapt.hpp>

#include "atomiverse/cpu/observer/mulliken.h"
#include "atomiverse/cuda/async_bridge.cuh"

#ifdef ATV_HAS_RHF
#include "atomiverse/cpu/scf/integrals.h"
#endif

namespace ATV {

#ifdef ATV_HAS_RHF

/// @brief Mulliken population observer for the CUDA RHF SCF pipeline.
///
/// Header-only — no .cu needed.  Pure observer (stores charges, no I/O).
/// Shares the pure C++ Mulliken implementation with the CPU path; the
/// only CUDA-specific responsibility is bringing the device density
/// matrix back to host via the async D2H bridge.
class MullikenObs : public AsyncBridge {
public:
    explicit MullikenObs(const cpu::scf::IntegralPack& ints, int every = 1)
        : integrals_(ints), every_(every) {}

    ~MullikenObs() override = default;

    std::string name() const override { return "mulliken"; }
    int frequency() const override { return every_; }

    const std::vector<double>& charges() const { return charges_; }
    double energy() const { return energy_val_; }
    int iterations() const { return iterations_; }
    int n_occ() const { return n_occ_; }

protected:
    void on_build(const BuildContext& ctx) override {
        n_basis_ = ctx.frame.scf.n;
        cudaMallocHost(&h_P_,      n_basis_ * n_basis_ * sizeof(double));
        cudaMallocHost(&h_energy_, sizeof(double));
        *h_energy_ = 0.0;
    }

    bool can_run(const LaunchContext& ctx) const override {
        return ctx.frame.scf.P != nullptr;
    }

    void issue_d2h(const LaunchContext& ctx) override {
        auto ds    = static_cast<cudaStream_t>(ctx.runtime.defer_stream);
        auto& f    = ctx.frame;
        const int n2 = n_basis_ * n_basis_;
        cudaMemcpyAsync(h_P_, f.scf.P, n2 * sizeof(double), cudaMemcpyDeviceToHost, ds);
        if (f.energy)
            cudaMemcpyAsync(h_energy_, f.energy, sizeof(double), cudaMemcpyDeviceToHost, ds);
    }

    void on_data_ready(const LaunchContext& ctx) override {
        // Adapt the pinned host buffer as a non-owning xtensor view.
        const std::size_t nb = static_cast<std::size_t>(n_basis_);
        auto P_view = xt::adapt(h_P_, nb * nb, xt::no_ownership(),
                                std::vector<std::size_t>{nb, nb});
        const xt::xtensor<double, 2> P = P_view;

        charges_    = cpu::mulliken_charge(integrals_, P);
        energy_val_ = *h_energy_;
        iterations_ = ctx.frame.scf.iteration;
        n_occ_      = static_cast<int>(ctx.frame.scf.n_occ);
    }

    void on_destroy(const DestroyContext&) override {
        if (h_P_)      { cudaFreeHost(h_P_);      h_P_      = nullptr; }
        if (h_energy_) { cudaFreeHost(h_energy_); h_energy_ = nullptr; }
    }

private:
    cpu::scf::IntegralPack integrals_;
    int every_;

    std::vector<double> charges_;
    double energy_val_ = 0;
    int iterations_  = 0;
    int n_occ_       = 0;

    int     n_basis_  = 0;
    double* h_P_      = nullptr;
    double* h_energy_ = nullptr;
};

#endif // ATV_HAS_RHF

} // namespace ATV
