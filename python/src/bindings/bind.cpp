/// @file python/bind.cpp
/// @brief Top-level nanobind module definition for _atomiverse.
///
/// Atomiverse binds the driver / component API only — never a Frame data
/// type. Simulation state crosses the boundary as a `molrs.Frame` (from the
/// separate `molrs` extension); each `run()` borrows it via a PyCapsule.

#include <nanobind/nanobind.h>

namespace nb = nanobind;

// Forward declarations — each bind_*.cpp provides one init function.
void bind_md(nb::module_& m);
void bind_pw(nb::module_& m);

NB_MODULE(_atomiverse, m) {
    m.doc() = "Atomiverse: molecular simulation compute backend";

    bind_md(m);
    bind_pw(m);
    // SCF and geometry-optimization bindings are future work (their own
    // spec): they will bind the Minimize driver + LBFGS/CG Optimizer
    // components and the SCFEvaluator/build_integrals/rhf_solve family in
    // the correct driver/component shape (arch-refactor-08 Out of scope).
}
