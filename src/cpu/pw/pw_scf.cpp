/// @file cpu/pw/pw_scf.cpp
/// @brief pw_scf() — ≤30-line shim constructing a transient PWSolver.
#include "atomiverse/cpu/pw/pw_scf.h"

namespace ATV {
namespace cpu {
namespace pw {

PWResult pw_scf(const cpu::Frame& frame, const std::string& upf_path,
                double ecut_ry, int k_mesh, double tol, int max_iter) {
    PWConfig config;
    config.upf_path = upf_path;
    config.ecut_ry = ecut_ry;
    config.k_mesh = k_mesh;
    config.tol = tol;
    config.max_iter = max_iter;

    PWSolver solver{std::move(config)};
    solver.build(frame);
    return solver.run(frame);
}

} // namespace pw
} // namespace cpu
} // namespace ATV
