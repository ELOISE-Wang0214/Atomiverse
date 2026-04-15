/// @file cpu/scf/orbital_info_observer.cpp
/// @brief OrbitalInfoObserver -- HOMO/LUMO/gap TSV writer.

#include "atomiverse/cpu/scf/orbital_info_observer.h"

#include <cstdio>
#include <stdexcept>

#include "atomiverse/core/units.h"

namespace ATV {
namespace cpu {
namespace scf {

OrbitalInfoObserver::OrbitalInfoObserver(SCFResultHandle handle,
                                         std::string filename,
                                         bool also_stdout)
    : handle_(std::move(handle)),
      filename_(std::move(filename)),
      also_stdout_(also_stdout) {
    if (!handle_.result || !handle_.ready)
        throw std::invalid_argument(
            "OrbitalInfoObserver: null SCFResultHandle (did you forget "
            "SCFEvaluator::result_handle()?)");
}

void OrbitalInfoObserver::build(const cpu::BuildContext& ctx) {
    if (ctx.runtime.mpi_rank != 0) return;
    if (!filename_.empty()) {
        os_.open(filename_, std::ios::out | std::ios::trunc);
        os_ << "step\thomo_eV\tlumo_eV\tgap_eV\n";
    }
    if (also_stdout_)
        std::printf("%6s  %12s  %12s  %12s\n", "step", "HOMO[eV]", "LUMO[eV]", "gap[eV]");
}

void OrbitalInfoObserver::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_rank != 0) return;
    if (!*handle_.ready) return;  // no SCF completed yet this step

    const SCFResult& r = *handle_.result;
    const int n_occ_a = r.n_alpha;
    const auto& eps = r.eigenvalues_alpha;
    const int n_mo = static_cast<int>(eps.size());

    if (n_occ_a <= 0 || n_occ_a > n_mo)
        throw std::runtime_error(
            "OrbitalInfoObserver: invalid n_alpha=" + std::to_string(n_occ_a) +
            " against n_mo=" + std::to_string(n_mo));

    // HOMO is the last occupied; LUMO is the first virtual (or HOMO itself
    // if there are no virtuals, which only happens at minimal basis on H-like
    // systems -- flag that explicitly).
    homo_eV_ = static_cast<double>(eps(n_occ_a - 1)) * kHartree2eV;
    lumo_eV_ = (n_occ_a < n_mo)
                   ? static_cast<double>(eps(n_occ_a)) * kHartree2eV
                   : homo_eV_;
    double gap = lumo_eV_ - homo_eV_;

    if (os_.is_open()) {
        os_ << ctx.step << '\t'
            << homo_eV_ << '\t'
            << lumo_eV_ << '\t'
            << gap << '\n';
        os_.flush();
    }
    if (also_stdout_)
        std::printf("%6d  %12.6f  %12.6f  %12.6f\n",
                    ctx.step, homo_eV_, lumo_eV_, gap);
    ++n_rows_;
}

void OrbitalInfoObserver::destroy(const cpu::DestroyContext& /*ctx*/) {
    if (os_.is_open()) os_.close();
}

} // namespace scf
} // namespace cpu
} // namespace ATV
