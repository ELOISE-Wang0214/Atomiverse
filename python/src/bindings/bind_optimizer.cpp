/// @file python/bind_optimizer.cpp
/// @brief L-BFGS geometry optimizer bindings.

#include "cpu_frame_wrapper.h"

#include <nanobind/stl/string.h>

#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/observer/thermo_component.h"
#include "atomiverse/facade/stage.h"

using namespace ATV;

void bind_optimizer(nb::module_& m) {

    // ── minimizer::Phase enum ──

    nb::enum_<minimizer::Phase>(m, "MinimizerPhase",
        "Phase sequence for geometry optimization.")
        .value("search_direction",  minimizer::Phase::search_direction)
        .value("line_search",       minimizer::Phase::line_search)
        .value("force",             minimizer::Phase::force)
        .value("convergence_check", minimizer::Phase::convergence_check)
        .value("end_of_step",       minimizer::Phase::end_of_step)
    ;

    // ── LBFGSStep component ──

    nb::class_<cpu::minimizer::LBFGSStep>(m, "LBFGSStep",
        "L-BFGS step component (search direction + convergence check).")
        .def(nb::init<int, double, double, int>(),
             nb::arg("memory") = 5,
             nb::arg("alpha0") = double(0.01),
             nb::arg("ftol") = double(0.05),
             nb::arg("max_ls") = 10,
             "Construct L-BFGS step.\n"
             "memory: correction pair history size.\n"
             "alpha0: initial step size [Angstrom^2/eV].\n"
             "ftol: force convergence tolerance [eV/Angstrom].\n"
             "max_ls: maximum line search iterations.")
        .def("converged", &cpu::minimizer::LBFGSStep::converged,
             "Whether optimization has converged.")
        .def("max_force", &cpu::minimizer::LBFGSStep::max_force,
             "Current maximum force magnitude [eV/Angstrom].")
    ;

    // ── SCFEvaluator component ──

    // SCFEvaluator is not directly bound since it contains unique_ptr
    // (non-copyable). Instead, it is constructed inline in add_scf_evaluator.

    // ── LBFGS minimizer driver ──

    using CpuLBFGS = cpu::minimizer::LBFGS;
    using MinPhaseBuilder = PhaseBuilder<minimizer::Phase, cpu::device_tag>;

    nb::class_<CpuLBFGS>(m, "LBFGS",
        "L-BFGS geometry optimizer driver (CPU backend).")
        .def(nb::init<>())
        .def("add_lbfgs_step", [](CpuLBFGS& self, cpu::minimizer::LBFGSStep comp) {
            self.add(std::move(comp));
        }, nb::arg("comp"),
           "Add an LBFGSStep component with default phase bindings.")
        .def("on", [](CpuLBFGS& self, minimizer::Phase phase) -> MinPhaseBuilder {
            return self.on(phase);
        }, nb::arg("phase"), nb::rv_policy::move,
           "Return a phase builder for explicit component binding.")
        .def("build", &CpuLBFGS::build,
             "Compile the execution plan.")
        .def("run", [](CpuLBFGS& self, nb::object frame_obj, int n_steps) {
            auto& cf = nb::cast<CpuFrame&>(frame_obj);
            { nb::gil_scoped_release release;
              self.run(cf.frame, n_steps); }
        }, nb::arg("frame"), nb::arg("n_steps"),
           "Run optimization. GIL is released during computation.")
    ;

    // ── MinimizerPhaseBuilder ──

    nb::class_<MinPhaseBuilder>(m, "MinimizerPhaseBuilder",
        "Fluent builder for explicit minimizer phase binding.")
        .def("add_scf_evaluator", [](MinPhaseBuilder& self,
                                      cpu::scf::SCFConfig config,
                                      int charge,
                                      int order) -> MinPhaseBuilder& {
            return self.add(cpu::scf::SCFEvaluator{config, charge}, order);
        }, nb::arg("config"), nb::arg("charge") = 0, nb::arg("order") = 0,
           nb::rv_policy::reference,
           "Add an SCFEvaluator to this phase.")
        .def("add_thermo", [](MinPhaseBuilder& self,
                               cpu::Thermo comp,
                               int order) -> MinPhaseBuilder& {
            return self.add(std::move(comp), order);
        }, nb::arg("comp"), nb::arg("order") = 0, nb::rv_policy::reference,
           "Add a Thermo observer to this phase.")
    ;
}
