/// @file regressions/scf/cuda-periodic-kpoint-hf/reg_cuda-periodic-kpoint-hf.cu
/// @brief Device KRHF on the same H₂ / STO-3G 1×1×2 cell as the CPU pin.

#include <cmath>
#include <cstdio>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cuda/scf/krhf.cuh"

int main() {
    using namespace ATV;
    const double L_A = 4.0 * units::kBohrInAngstrom;
    const double z_A = 1.4 * units::kBohrInAngstrom;
    cpu::Frame frame(2);
    frame.span(kAtom_type)[0] = frame.span(kAtom_type)[1] = 1;
    frame.span(kAtom_x)[0] = frame.span(kAtom_y)[0] = frame.span(kAtom_z)[0] =
        0.0;
    frame.span(kAtom_x)[1] = frame.span(kAtom_y)[1] = 0.0;
    frame.span(kAtom_z)[1] = z_A;
    frame.box = cpu::Box::cubic(L_A);

    cuda::scf::KRHF krhf{cpu::scf::BasisSet::STO3G, {1, 1, 2}, 0, 1e-10, 40};
    krhf.run(frame);
    if (!krhf.converged()) {
        std::fprintf(stderr, "FAIL: device KRHF did not converge\n");
        return 1;
    }
    constexpr double kEpinHa = -1.3541438798330374;
    const double diff = std::fabs(krhf.energy() - kEpinHa);
    std::printf("E_KRHF_cuda got=%.16f |diff|=%.3e iters=%d\n", krhf.energy(),
                diff, krhf.iterations());
    if (diff > 1e-6) {
        std::fprintf(stderr, "FAIL: device KRHF energy off pin\n");
        return 1;
    }
    std::printf("PASS\n");
    return 0;
}
