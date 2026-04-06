/// @file python/bind.cpp
/// @brief Top-level nanobind module definition for _atomiverse.

#include <nanobind/nanobind.h>

namespace nb = nanobind;

// Forward declarations — each bind_*.cpp provides one init function.
void bind_frame(nb::module_& m);
void bind_md(nb::module_& m);
void bind_scf(nb::module_& m);
void bind_optimizer(nb::module_& m);

NB_MODULE(_atomiverse, m) {
    m.doc() = "Atomiverse: molecular simulation compute backend";

    bind_frame(m);
    bind_md(m);
    bind_scf(m);
    bind_optimizer(m);
}
