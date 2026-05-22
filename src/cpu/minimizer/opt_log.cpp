/// @file cpu/minimizer/opt_log.cpp
/// @brief OptLog -- TSV writer for geometry-optimization progress.

#include "atomiverse/cpu/minimizer/opt_log.h"

#include <cmath>
#include <cstdio>

namespace ATV {
namespace cpu {
namespace minimizer {

OptLog::OptLog(std::string filename, bool also_stdout)
    : filename_(std::move(filename)), also_stdout_(also_stdout) {}

void OptLog::build(const cpu::BuildContext& ctx) {
    if (ctx.runtime.mpi_rank != 0)
        return;
    os_.open(filename_, std::ios::out | std::ios::trunc);
    os_ << "step\tenergy_eV\tmax_force\trms_force\tstep_norm\n";
    if (also_stdout_)
        std::printf("%6s  %15s  %12s  %12s  %12s\n", "step", "E[eV]", "max_F",
                    "rms_F", "|dx|");
}

void OptLog::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_rank != 0)
        return;

    const auto& f = ctx.frame;
    const int n = f.atoms.size();
    double sumsq = 0, maxsq = 0;
    for (int i = 0; i < n; ++i) {
        double fx = f.atoms.fx(i);
        double fy = f.atoms.fy(i);
        double fz = f.atoms.fz(i);
        double m2 = fx * fx + fy * fy + fz * fz;
        sumsq += m2;
        if (m2 > maxsq)
            maxsq = m2;
    }
    double max_force = std::sqrt(maxsq);
    double rms_force = (n > 0) ? std::sqrt(sumsq / n) : 0.0;

    double step_norm = 0;
    if (has_prev_) {
        double acc = 0;
        for (int i = 0; i < n; ++i) {
            double dx = f.atoms.x(i) - x_prev_[i];
            double dy = f.atoms.y(i) - y_prev_[i];
            double dz = f.atoms.z(i) - z_prev_[i];
            acc += dx * dx + dy * dy + dz * dz;
        }
        step_norm = std::sqrt(acc);
    }

    x_prev_.resize(n);
    y_prev_.resize(n);
    z_prev_.resize(n);
    for (int i = 0; i < n; ++i) {
        x_prev_[i] = f.atoms.x(i);
        y_prev_[i] = f.atoms.y(i);
        z_prev_[i] = f.atoms.z(i);
    }
    has_prev_ = true;

    if (os_.is_open()) {
        os_ << ctx.step << '\t' << f.scalars.energy << '\t' << max_force << '\t'
            << rms_force << '\t' << step_norm << '\n';
        os_.flush();
    }
    if (also_stdout_) {
        std::printf("%6d  %15.8e  %12.6e  %12.6e  %12.6e\n", ctx.step,
                    f.scalars.energy, max_force, rms_force, step_norm);
    }
    ++n_rows_;
}

void OptLog::destroy(const cpu::DestroyContext& /*ctx*/) {
    if (os_.is_open())
        os_.close();
}

} // namespace minimizer
} // namespace cpu
} // namespace ATV
