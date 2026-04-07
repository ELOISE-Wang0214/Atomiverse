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

#include "atomiverse/core/mat3.h"

#include "atomiverse/cpu/frame.h"
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
        auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
        auto& frame = ctx.frame;
        const int n = frame.n_atoms();
        n_basis_    = frame.scf.n;

        // Synchronous D2H of geometry at build time.
        std::vector<int>  tid(n);
        std::vector<double> xh(n), yh(n), zh(n);
        cudaMemcpyAsync(tid.data(), frame.atoms.type, n * sizeof(int),  cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(xh.data(),  frame.atoms.x,       n * sizeof(double), cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(yh.data(),  frame.atoms.y,       n * sizeof(double), cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(zh.data(),  frame.atoms.z,       n * sizeof(double), cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);

        cpu_frame_ = std::make_unique<cpu::Frame>(n);
        for (int i = 0; i < n; ++i) {
            cpu_frame_->type_id[i] = tid[i];
            cpu_frame_->x[i] = xh[i];
            cpu_frame_->y[i] = yh[i];
            cpu_frame_->z[i] = zh[i];
        }
        cpu_frame_->box = Mat3::cubic(double(100));

        cudaMallocHost(&h_P_,      n_basis_ * n_basis_ * sizeof(double));
        cudaMallocHost(&h_energy_, sizeof(double));
        *h_energy_ = double(0);
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
        charges_    = cpu::mulliken_charge(*cpu_frame_, h_P_, integrals_.S.data(), n_basis_);
        energy_val_ = *h_energy_;
        iterations_ = ctx.frame.scf.iteration;
        n_occ_      = static_cast<int>(ctx.frame.scf.n_occ);
    }

    void on_destroy(const DestroyContext&) override {
        if (h_P_)      { cudaFreeHost(h_P_);      h_P_      = nullptr; }
        if (h_energy_) { cudaFreeHost(h_energy_); h_energy_ = nullptr; }
        cpu_frame_.reset();
    }

private:
    cpu::scf::IntegralPack integrals_;
    int every_;

    std::unique_ptr<cpu::Frame> cpu_frame_;
    std::vector<double> charges_;
    double energy_val_ = 0;
    int iterations_  = 0;
    int n_occ_       = 0;

    int   n_basis_  = 0;
    double* h_P_      = nullptr;
    double* h_energy_ = nullptr;
};

#endif // ATV_HAS_RHF

} // namespace ATV
