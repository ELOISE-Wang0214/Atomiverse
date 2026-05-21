/// @file cuda/scf/driver.cu
/// @brief SCF ZeroDensity infrastructure component.

#include "atomiverse/cuda/scf/scf.cuh"

#include "atomiverse/cuda/types.cuh"
#include <cuda_runtime.h>

namespace ATV {

using cuda::prec_t;
using cuda::real_t;
using cuda::work_t;

namespace scf {

void ZeroDensity::launch(const LaunchContext& ctx) {
    // Only zero on the first iteration -- subsequent iterations
    // get P from DensityBuilder.
    if (ctx.step != 0)
        return;

    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& scf = ctx.frame.scf;
    if (scf.n > 0) {
        int nn = scf.n * scf.n;
        cudaMemsetAsync(scf.P, 0, nn * sizeof(prec_t), stream);
    }
}

} // namespace scf
} // namespace ATV
