/// @file python/bind_pw.cpp
/// @brief Plane-wave DFT bindings: pw::pw_scf()

#include "molrs_frame_bridge.h"

#include <nanobind/stl/string.h>

#ifdef ATV_HAS_PW
#include "atomiverse/cpu/materialize.h"
#include "atomiverse/cpu/pw/pw_scf.h"
#endif

void bind_pw(nb::module_& m) {

#ifdef ATV_HAS_PW
    using namespace ATV;
    using namespace ATV::cpu::pw;

    nb::class_<PWResult>(m, "PWResult",
                         "Result of a plane-wave DFT calculation.")
        .def_prop_ro(
            "energy", [](const PWResult& r) { return r.energy; },
            "Total energy [Ry].")
        .def_prop_ro(
            "converged", [](const PWResult& r) { return r.converged; },
            "Whether SCF converged within max_iter.")
        .def_prop_ro(
            "iterations", [](const PWResult& r) { return r.iterations; },
            "Number of SCF iterations taken.")
        .def("__repr__", [](const PWResult& r) {
            std::string conv = r.converged ? "converged" : "NOT converged";
            return "PWResult(" + conv + ", E=" + std::to_string(r.energy)
                   + " Ry, iters=" + std::to_string(r.iterations) + ")";
        });

    m.def(
        "pw_scf",
        [](nb::object frame_obj, const std::string& upf_path, double ecut_ry,
           int k_mesh, double tol, int max_iter) {
            auto fref = borrow_molrs_frame(frame_obj);
            PWResult result;
            {
                nb::gil_scoped_release release;
                cpu::Frame work = cpu::materialize(*fref);
                result = pw_scf(work, upf_path, ecut_ry, k_mesh, tol, max_iter);
            }
            return result;
        },
        nb::arg("frame"), nb::arg("upf_path") = "",
        nb::arg("ecut_ry") = double(30.0), nb::arg("k_mesh") = 4,
        nb::arg("tol") = double(1e-6), nb::arg("max_iter") = 50,
        "Run plane-wave DFT total energy calculation for a periodic system.\n\n"
        "The Frame provides atom positions, types (Z), and simulation cell.\n"
        "Returns total energy in Rydberg.\n\n"
        "Parameters\n"
        "----------\n"
        "frame : Frame\n"
        "    Input structure with positions, atom types, and cell.\n"
        "upf_path : str\n"
        "    Path to UPF pseudopotential file (empty = toy local potential).\n"
        "ecut_ry : float\n"
        "    Plane-wave energy cutoff [Ry]. Default 30.0.\n"
        "k_mesh : int\n"
        "    Monkhorst-Pack mesh size per dimension. Default 4.\n"
        "tol : float\n"
        "    SCF convergence tolerance. Default 1e-6.\n"
        "max_iter : int\n"
        "    Maximum SCF iterations. Default 50.\n"
        "\n"
        "Returns\n"
        "-------\n"
        "PWResult\n"
        "    Result with ``.energy`` (Ry), ``.converged``, and "
        "``.iterations``.");
#endif
}
