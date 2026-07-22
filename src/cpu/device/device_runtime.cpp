/// @file cpu/device/device_runtime.cpp
/// @brief cpu::device::DeviceRuntime implementation (moved out of the former
///        CPU execution engine by backend-mirror-05-dissolve).

#include "atomiverse/cpu/device/device_runtime.h"

#include "atomiverse/cpu/component/context.h"

#include <cstdio>
#include <cstdlib>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {
namespace device {

void DeviceRuntime::set_threads(Context& ctx, int n) {
    ctx.n_threads = (n < 1) ? 0 : n;
}

void DeviceRuntime::configure_threads(Context& ctx) {
#ifdef _OPENMP
    if (ctx.n_threads > 0)
        omp_set_num_threads(ctx.n_threads);
    else
        ctx.n_threads = omp_get_max_threads();
#else
    ctx.n_threads = 1;
#endif

    // Warn once per process if the user will pay NUMA / cache penalties for
    // letting the OS migrate threads. Load-bearing on 2-socket / NUMA boxes.
    static bool s_omp_bind_warned = false;
    if (!s_omp_bind_warned && ctx.n_threads > 1) {
        const char* bind = std::getenv("OMP_PROC_BIND");
        const bool unset = (bind == nullptr) || std::string(bind).empty();
        const bool unbinding = unset || std::string(bind) == "false"
                               || std::string(bind) == "FALSE";
        if (unbinding) {
            int my_rank = 0;
#ifdef ATV_HAS_MPI
            int mpi_init_check = 0;
            MPI_Initialized(&mpi_init_check);
            if (mpi_init_check)
                MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
            if (my_rank == 0) {
                std::fprintf(
                    stderr,
                    "[atv] warning: OMP_PROC_BIND is unset and "
                    "n_threads=%d > 1; consider 'export OMP_PROC_BIND=close "
                    "OMP_PLACES=cores' for stable hot-path throughput "
                    "(see benchmarks/md/pe/validation/reference/SOURCES.md "
                    "'Thread binding').\n",
                    ctx.n_threads);
                std::fflush(stderr);
            }
            s_omp_bind_warned = true;
        }
    }
}

void DeviceRuntime::init_timing(Context& ctx) {
    // Per-component timing is ON BY DEFAULT (near-zero overhead) so every driver
    // can emit the built-in end-of-run time+memory summary. Opt OUT with
    // ATV_TIMER=0/off. Re-apply the remembered preference so a pre-run
    // enable_timing(false) sticks; ATV_TIMER overrides.
    timing_.set_enabled(timing_pref_);
    if (const char* tv = std::getenv("ATV_TIMER")) {
        const std::string v = tv;
        if (v == "0" || v == "false" || v == "FALSE" || v == "off")
            timing_.set_enabled(false);
        else if (v == "1" || v == "true" || v == "TRUE" || v == "on")
            timing_.set_enabled(true);
    }
    ctx.timer = &timing_;
}

} // namespace device
} // namespace cpu
} // namespace ATV
