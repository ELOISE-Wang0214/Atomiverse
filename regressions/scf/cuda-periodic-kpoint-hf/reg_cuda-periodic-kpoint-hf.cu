/// @file regressions/scf/cuda-periodic-kpoint-hf/reg_cuda-periodic-kpoint-hf.cu
/// @brief Device KRHF on the H₂ / STO-3G 4 Bohr cell (1×1×2 and 2×1×1).

#include <cmath>
#include <cstdio>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cuda/scf/krhf.cuh"

namespace {

ATV::cpu::Frame h2_cell() {
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
    return frame;
}

int gate(const char* name, double got, double want, bool conv) {
    const double diff = std::fabs(got - want);
    std::printf("%s got=%.16f |diff|=%.3e conv=%d\n", name, got, diff,
                conv ? 1 : 0);
    if (!conv || diff > 1e-6) {
        std::fprintf(stderr, "FAIL: %s\n", name);
        return 1;
    }
    return 0;
}

} // namespace

int main() {
    using namespace ATV;
    int rc = 0;
    {
        cpu::Frame frame = h2_cell();
        cuda::scf::KRHF krhf{
            cpu::scf::BasisSet::STO3G, {1, 1, 2}, 0, 1e-10, 40};
        krhf.run(frame);
        rc |= gate("KRHF cuda 1x1x2", krhf.energy(), -1.3541438798330374,
                   krhf.converged());
    }
    {
        cpu::Frame frame = h2_cell();
        cuda::scf::KRHF krhf{
            cpu::scf::BasisSet::STO3G, {2, 1, 1}, 0, 1e-10, 40};
        krhf.run(frame);
        rc |= gate("KRHF cuda 2x1x1", krhf.energy(), -1.4471748215213558,
                   krhf.converged());
    }
    if (rc == 0)
        std::printf("PASS\n");
    return rc == 0 ? 0 : 1;
}
