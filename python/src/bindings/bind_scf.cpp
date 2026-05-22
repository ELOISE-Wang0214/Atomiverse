/// @file python/bind_scf.cpp
/// @brief Unified SCF engine bindings: integrals, scf_solve, scf_gradient.
///
/// Breaking change from v1: replaces per-method solvers (rhf_solve, uhf_solve,
/// uks_solve) with the unified scf_solve() + convenience config factories.

#include "molrs_frame_bridge.h"

#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cstring>

#include "atomiverse/cpu/materialize.h"
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;
using namespace ATV::cpu::scf;

// ── xtensor → numpy helpers ────────────────────────────────────────────────

/// Copy a 1D xtensor to a numpy array.
static nb::ndarray<nb::numpy, double, nb::ndim<1>>
xt1d_to_numpy(const xt::xtensor<double, 1>& t) {
    auto n = static_cast<size_t>(t.shape(0));
    double* buf = new double[n];
    std::memcpy(buf, t.data(), n * sizeof(double));
    nb::capsule owner(
        buf, [](void* p) noexcept { delete[] static_cast<double*>(p); });
    return nb::ndarray<nb::numpy, double, nb::ndim<1>>(buf, {n}, owner);
}

/// Copy a 2D xtensor to a numpy array (row-major).
static nb::ndarray<nb::numpy, double, nb::ndim<2>>
xt2d_to_numpy(const xt::xtensor<double, 2>& t) {
    auto s0 = static_cast<size_t>(t.shape(0));
    auto s1 = static_cast<size_t>(t.shape(1));
    double* buf = new double[s0 * s1];
    std::memcpy(buf, t.data(), s0 * s1 * sizeof(double));
    nb::capsule owner(
        buf, [](void* p) noexcept { delete[] static_cast<double*>(p); });
    return nb::ndarray<nb::numpy, double, nb::ndim<2>>(buf, {s0, s1}, owner);
}

/// Copy a 4D xtensor to a numpy array.
static nb::ndarray<nb::numpy, double, nb::ndim<4>>
xt4d_to_numpy(const xt::xtensor<double, 4>& t) {
    size_t s[4];
    for (int i = 0; i < 4; ++i)
        s[i] = static_cast<size_t>(t.shape(i));
    size_t total = s[0] * s[1] * s[2] * s[3];
    double* buf = new double[total];
    std::memcpy(buf, t.data(), total * sizeof(double));
    nb::capsule owner(
        buf, [](void* p) noexcept { delete[] static_cast<double*>(p); });
    return nb::ndarray<nb::numpy, double, nb::ndim<4>>(
        buf, {s[0], s[1], s[2], s[3]}, owner);
}

void bind_scf(nb::module_& m) {

    // ── Enums ──

    nb::enum_<SpinType>(m, "SpinType", "Spin treatment for SCF.")
        .value("Restricted", SpinType::Restricted)
        .value("Unrestricted", SpinType::Unrestricted);

    nb::enum_<MeanFieldType>(m, "MeanFieldType", "Mean-field method type.")
        .value("HF", MeanFieldType::HF)
        .value("KS", MeanFieldType::KS);

    nb::enum_<XCType>(m, "XCType", "Exchange-correlation functional type.")
        .value("LSDA", XCType::LSDA)
        .value("PBE", XCType::PBE);

    // ── SCFConfig ──

    nb::class_<SCFConfig>(m, "SCFConfig",
                          "Configuration for the unified SCF solver.")
        .def(nb::init<>())
        .def_rw("spin", &SCFConfig::spin)
        .def_rw("method", &SCFConfig::method)
        .def_rw("xc_type", &SCFConfig::xc_type)
        .def_rw("n_alpha", &SCFConfig::n_alpha)
        .def_rw("n_beta", &SCFConfig::n_beta)
        .def_rw("tol", &SCFConfig::tol)
        .def_rw("max_iter", &SCFConfig::max_iter)
        .def_rw("use_diis", &SCFConfig::use_diis)
        .def_rw("diis_max_vectors", &SCFConfig::diis_max_vectors)
        .def_rw("diis_start_iter", &SCFConfig::diis_start_iter)
        .def_rw("damping_factor", &SCFConfig::damping_factor)
        .def_rw("n_radial", &SCFConfig::n_radial)
        .def_rw("n_angular", &SCFConfig::n_angular)
        .def("__repr__", [](const SCFConfig& c) {
            std::string spin_str = (c.spin == SpinType::Restricted)
                                       ? "Restricted"
                                       : "Unrestricted";
            std::string method_str =
                (c.method == MeanFieldType::HF) ? "HF" : "KS";
            return "SCFConfig(spin=" + spin_str + ", method=" + method_str
                   + ", n_alpha=" + std::to_string(c.n_alpha)
                   + ", n_beta=" + std::to_string(c.n_beta)
                   + ", tol=" + std::to_string(c.tol) + ")";
        });

    // ── SCFResult ──

    nb::class_<SCFResult>(m, "SCFResult",
                          "Unified SCF result with full energy decomposition.")
        // Core
        .def_prop_ro("converged",
                     [](const SCFResult& r) { return r.converged; })
        .def_prop_ro("iterations",
                     [](const SCFResult& r) { return r.iterations; })
        .def_prop_ro("energy", [](const SCFResult& r) { return r.energy; })
        .def_prop_ro("residual", [](const SCFResult& r) { return r.residual; })
        // Energy decomposition
        .def_prop_ro("one_electron_energy",
                     [](const SCFResult& r) { return r.one_electron_energy; })
        .def_prop_ro("coulomb_energy",
                     [](const SCFResult& r) { return r.coulomb_energy; })
        .def_prop_ro("exchange_energy",
                     [](const SCFResult& r) { return r.exchange_energy; })
        .def_prop_ro("xc_energy",
                     [](const SCFResult& r) { return r.xc_energy; })
        .def_prop_ro("nuclear_repulsion",
                     [](const SCFResult& r) { return r.nuclear_repulsion; })
        // Spin info
        .def_prop_ro("n_alpha", [](const SCFResult& r) { return r.n_alpha; })
        .def_prop_ro("n_beta", [](const SCFResult& r) { return r.n_beta; })
        .def_prop_ro("spin_contamination",
                     [](const SCFResult& r) { return r.spin_contamination; })
        // Orbitals and densities (as numpy arrays)
        .def("get_eigenvalues_alpha",
             [](const SCFResult& r) {
                 return xt1d_to_numpy(r.eigenvalues_alpha);
             })
        .def("get_eigenvalues_beta",
             [](const SCFResult& r) {
                 return xt1d_to_numpy(r.eigenvalues_beta);
             })
        .def("get_P_alpha",
             [](const SCFResult& r) { return xt2d_to_numpy(r.P_alpha); })
        .def("get_P_beta",
             [](const SCFResult& r) { return xt2d_to_numpy(r.P_beta); })
        .def("get_C_alpha",
             [](const SCFResult& r) { return xt2d_to_numpy(r.C_alpha); })
        .def("get_C_beta",
             [](const SCFResult& r) { return xt2d_to_numpy(r.C_beta); })
        // Grid quality (KS only)
        .def_prop_ro("integrated_n_alpha",
                     [](const SCFResult& r) { return r.integrated_n_alpha; })
        .def_prop_ro("integrated_n_beta",
                     [](const SCFResult& r) { return r.integrated_n_beta; })
        // Diagnostics
        .def_prop_ro("trace",
                     [](const SCFResult& r) {
                         nb::list trace;
                         for (const auto& rec : r.trace) {
                             nb::dict d;
                             d["iter"] = rec.iter;
                             d["energy"] = rec.energy;
                             d["residual"] = rec.residual;
                             d["diis_active"] = rec.diis_active;
                             trace.append(d);
                         }
                         return trace;
                     })
        .def("__repr__", [](const SCFResult& r) {
            std::string conv = r.converged ? "converged" : "NOT converged";
            return "SCFResult(" + conv + ", E=" + std::to_string(r.energy)
                   + " Ha, iters=" + std::to_string(r.iterations) + ")";
        });

    // ── GradientResult ──

    nb::class_<GradientResult>(m, "GradientResult",
                               "Result of an analytic gradient computation.")
        .def(
            "get_forces",
            [](const GradientResult& r) { return xt2d_to_numpy(r.forces); },
            "Forces [n_atoms, 3] in eV/Angstrom.")
        .def_prop_ro(
            "energy", [](const GradientResult& r) { return r.energy; },
            "Total energy in eV.")
        .def("__repr__", [](const GradientResult& r) {
            auto n = r.forces.shape(0);
            return "GradientResult(n_atoms=" + std::to_string(n)
                   + ", E=" + std::to_string(r.energy) + " eV)";
        });

    // ── IntegralPack ──

    nb::class_<IntegralPack>(
        m, "IntegralPack",
        "Molecular integrals (S, T, V, Hcore, ERI, V_nuc) for HF/DFT.")
        .def("get_S",
             [](const IntegralPack& ip) { return xt2d_to_numpy(ip.S); })
        .def("get_T",
             [](const IntegralPack& ip) { return xt2d_to_numpy(ip.T); })
        .def("get_V",
             [](const IntegralPack& ip) { return xt2d_to_numpy(ip.V); })
        .def("get_Hcore",
             [](const IntegralPack& ip) { return xt2d_to_numpy(ip.Hcore); })
        .def("get_ERI",
             [](const IntegralPack& ip) { return xt4d_to_numpy(ip.ERI); })
        .def_prop_ro("V_nuc", [](const IntegralPack& ip) { return ip.V_nuc; })
        .def_prop_ro("n_basis",
                     [](const IntegralPack& ip) { return ip.n_basis; })
        .def_prop_ro("n_electrons",
                     [](const IntegralPack& ip) { return ip.n_electrons; })
        .def_prop_ro("n_occ",
                     [](const IntegralPack& ip) { return ip.n_occ(); });

    // ── build_integrals_sto3g ──

    m.def(
        "build_integrals_sto3g",
        [](nb::object frame_obj, int charge) {
            auto fref = borrow_molrs_frame(frame_obj);
            IntegralPack ip;
            {
                nb::gil_scoped_release release;
                cpu::Frame work = cpu::materialize(*fref);
                ip = build_integrals(work, BasisSet::STO3G, charge, false);
            }
            return ip;
        },
        nb::arg("frame"), nb::arg("charge") = 0,
        "Compute STO-3G molecular integrals from a molrs.Frame.\n"
        "atoms.type must be atomic number Z (1=H, 2=He, 3=Li).");

    // ── scf_solve ──

    m.def(
        "scf_solve",
        [](const IntegralPack& ip, const SCFConfig& config) {
            SCFResult result;
            {
                nb::gil_scoped_release release;
                result = scf_solve(ip, config);
            }
            return result;
        },
        nb::arg("integrals"), nb::arg("config"),
        "Run the unified SCF solver. Method/spin dispatched via SCFConfig.");

    // ── scf_gradient ──

    m.def(
        "scf_gradient",
        [](nb::object frame_obj, const IntegralPack& ip,
           const SCFResult& result, const SCFConfig& config) {
            auto fref = borrow_molrs_frame(frame_obj);
            GradientResult grad;
            {
                nb::gil_scoped_release release;
                cpu::Frame work = cpu::materialize(*fref);
                grad = scf_gradient(work, ip, result, config);
            }
            return grad;
        },
        nb::arg("frame"), nb::arg("integrals"), nb::arg("result"),
        nb::arg("config"),
        "Compute analytic nuclear gradient (forces in eV/Angstrom).");

    // ── Convenience config factories ──

    m.def(
        "rhf_config",
        [](const IntegralPack& ip, double tol) { return rhf_config(ip, tol); },
        nb::arg("integrals"), nb::arg("tol") = double(1e-8),
        "Create an SCFConfig for restricted Hartree-Fock.");

    m.def(
        "uhf_config",
        [](const IntegralPack& ip, int n_alpha, int n_beta, double tol) {
            return uhf_config(ip, n_alpha, n_beta, tol);
        },
        nb::arg("integrals"), nb::arg("n_alpha"), nb::arg("n_beta"),
        nb::arg("tol") = double(1e-8),
        "Create an SCFConfig for unrestricted Hartree-Fock.");

    m.def(
        "uks_config",
        [](const IntegralPack& ip, int n_alpha, int n_beta, XCType xc,
           double tol) { return uks_config(ip, n_alpha, n_beta, xc, tol); },
        nb::arg("integrals"), nb::arg("n_alpha"), nb::arg("n_beta"),
        nb::arg("xc") = XCType::LSDA, nb::arg("tol") = double(1e-8),
        "Create an SCFConfig for unrestricted Kohn-Sham DFT.");
}
