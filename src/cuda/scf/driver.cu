/// @file cuda/scf/driver.cu
/// @brief SCF ZeroDensity infrastructure component.

#include "atomiverse/cuda/scf/scf.cuh"

#include <cuda_runtime.h>

namespace ATV {
namespace scf {

void ZeroDensity::launch(const LaunchContext& ctx) {
    // Only zero on the first iteration -- subsequent iterations
    // get P from DensityBuilder.
    if (ctx.step != 0) return;

    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& scf = ctx.frame.scf;
    if (scf.n > 0) {
        Int nn = scf.n * scf.n;
        cudaMemsetAsync(scf.P, 0, nn * sizeof(Real), stream);
    }
}

} // namespace scf
} // namespace ATV
