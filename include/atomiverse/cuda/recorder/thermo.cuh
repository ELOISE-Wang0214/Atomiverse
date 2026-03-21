/// @file cuda/recorder/thermo.cuh
/// @brief Thermo -- header-only thermodynamic output recorder (async D2H bridge).
///
/// @see cpu::thermo_print (cpu/recorder/thermo.h)
#pragma once

#include <fstream>
#include <string>
#include <utility>

#include <cuda_runtime.h>

#include "atomiverse/cpu/recorder/thermo.h"
#include "atomiverse/cuda/async_bridge.cuh"

namespace ATV {

/// @brief Thermodynamic output (LAMMPS thermo_style one).
///
/// @code{.cpp}
/// sim.on(md::Phase::end_of_step).add(Thermo{100});               // stdout only
/// sim.on(md::Phase::end_of_step).add(Thermo{100, "thermo.dat"}); // stdout + file
/// @endcode
class Thermo : public AsyncBridge {
public:
    explicit Thermo(int every = 1, std::string filename = "")
        : every_(every), filename_(std::move(filename)) {}

    std::string name() const override { return "thermo"; }
    int frequency() const override { return every_; }

protected:
    void on_build(const BuildContext& ctx) override {
        n_atoms_        = ctx.frame.n_atoms();
        header_printed_ = false;
        cudaMallocHost(&h_pe_, sizeof(Real));
        cudaMallocHost(&h_ke_, sizeof(Real));
        *h_pe_ = Real(0);
        *h_ke_ = Real(0);
        if (!filename_.empty()) {
            file_.open(filename_);
            if (file_.is_open())
                file_ << "# Step Temp PotEng KinEng TotEng\n";
        }
    }

    void issue_d2h(const LaunchContext& ctx) override {
        auto ds = static_cast<cudaStream_t>(ctx.runtime.defer_stream);
        cudaMemcpyAsync(h_pe_, ctx.frame.energy,  sizeof(Real), cudaMemcpyDeviceToHost, ds);
        cudaMemcpyAsync(h_ke_, ctx.frame.kinetic, sizeof(Real), cudaMemcpyDeviceToHost, ds);
    }

    void on_data_ready(const LaunchContext& ctx) override {
        cpu::thermo_print(ctx.step, *h_pe_, *h_ke_, n_atoms_, header_printed_, file_);
    }

    void on_destroy(const DestroyContext&) override {
        if (file_.is_open()) file_.close();
        if (h_pe_) { cudaFreeHost(h_pe_); h_pe_ = nullptr; }
        if (h_ke_) { cudaFreeHost(h_ke_); h_ke_ = nullptr; }
    }

private:
    int every_;
    std::string filename_;
    int  n_atoms_         = 0;
    bool header_printed_  = false;
    std::ofstream file_;
    Real* h_pe_ = nullptr;
    Real* h_ke_ = nullptr;
};

} // namespace ATV
