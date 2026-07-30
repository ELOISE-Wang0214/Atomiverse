/// @file cuda/scf/scf_evaluator.cu
/// @brief Synthetic-path cuda::scf::SCF only (no integral-pack ctors).

#include "atomiverse/cuda/scf/scf_evaluator.cuh"

namespace ATV {
namespace cuda {
namespace scf {

SCF::SCF(int n_occ, double alpha, double tol, int max_iter)
    : solver_(std::make_unique<SCFSolver>(n_occ, alpha, tol, max_iter)) {}

void SCF::build(const BuildContext& ctx) {
    solver_->build(ctx.runtime.stream, ctx.frame);
}

void SCF::launch(const LaunchContext& ctx) {
    solver_->run(ctx.runtime.stream, ctx.frame);
}

void SCF::destroy(const DestroyContext& ctx) {
    solver_->destroy(ctx.runtime.stream, ctx.frame);
}

} // namespace scf
} // namespace cuda
} // namespace ATV
