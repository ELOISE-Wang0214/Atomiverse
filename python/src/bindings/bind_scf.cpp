/// @file python/bind_scf.cpp
/// @brief SCF solver bindings: integrals, RHF, UHF, UKS.

#include "cpu_frame_wrapper.h"

#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cpu/scf/uhf_solve.h"
#include "atomiverse/cpu/scf/uks_solve.h"

using namespace ATV::cpu::scf;

/// Copy a 2D xtensor to a numpy array (row-major).
static nb::ndarray<nb::numpy, Real, nb::ndim<2>>
xt2d_to_numpy(const xt::xtensor<Real, 2>& t) {
    auto s0 = static_cast<size_t>(t.shape(0));
    auto s1 = static_cast<size_t>(t.shape(1));
    Real* buf = new Real[s0 * s1];
    std::memcpy(buf, t.data(), s0 * s1 * sizeof(Real));
    nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<Real*>(p); });
    return nb::ndarray<nb::numpy, Real, nb::ndim<2>>(buf, {s0, s1}, owner);
}

/// Copy a 4D xtensor to a numpy array.
static nb::ndarray<nb::numpy, Real, nb::ndim<4>>
xt4d_to_numpy(const xt::xtensor<Real, 4>& t) {
    size_t s[4];
    for (int i = 0; i < 4; ++i) s[i] = static_cast<size_t>(t.shape(i));
    size_t total = s[0] * s[1] * s[2] * s[3];
    Real* buf = new Real[total];
    std::memcpy(buf, t.data(), total * sizeof(Real));
    nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<Real*>(p); });
    return nb::ndarray<nb::numpy, Real, nb::ndim<4>>(buf, {s[0], s[1], s[2], s[3]}, owner);
}

/// Copy a 1D xtensor to a numpy array.
static nb::ndarray<nb::numpy, Real, nb::ndim<1>>
xt1d_to_numpy(const xt::xtensor<Real, 1>& t) {
    auto n = static_cast<size_t>(t.shape(0));
    Real* buf = new Real[n];
    std::memcpy(buf, t.data(), n * sizeof(Real));
    nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<Real*>(p); });
    return nb::ndarray<nb::numpy, Real, nb::ndim<1>>(buf, {n}, owner);
}


void bind_scf(nb::module_& m) {

    // ── IntegralPack ──
    nb::class_<IntegralPack>(m, "IntegralPack",
        "Molecular integrals (S, T, V, Hcore, ERI, V_nuc) for HF/DFT.")
        .def("get_S", [](const IntegralPack& ip) { return xt2d_to_numpy(ip.S); })
        .def("get_T", [](const IntegralPack& ip) { return xt2d_to_numpy(ip.T); })
        .def("get_V", [](const IntegralPack& ip) { return xt2d_to_numpy(ip.V); })
        .def("get_Hcore", [](const IntegralPack& ip) { return xt2d_to_numpy(ip.Hcore); })
        .def("get_ERI", [](const IntegralPack& ip) { return xt4d_to_numpy(ip.ERI); })
        .def_prop_ro("V_nuc", [](const IntegralPack& ip) { return ip.V_nuc; })
        .def_prop_ro("n_basis", [](const IntegralPack& ip) { return ip.n_basis; })
        .def_prop_ro("n_electrons", [](const IntegralPack& ip) { return ip.n_electrons; })
        .def_prop_ro("n_occ", [](const IntegralPack& ip) { return ip.n_occ(); })
    ;

    // ── build_integrals_sto3g ──
    m.def("build_integrals_sto3g",
        [](nb::object frame_obj, int charge) {
            auto& cf = nb::cast<CpuFrame&>(frame_obj);
            IntegralPack ip;
            { nb::gil_scoped_release release;
              ip = build_integrals_sto3g(cf.frame, charge); }
            return ip;
        },
        nb::arg("frame"), nb::arg("charge") = 0,
        "Compute STO-3G molecular integrals from a CpuFrame.\n"
        "atoms.type must be atomic number Z (1=H, 2=He, 3=Li)."
    );

    // ── RHFResult ──
    nb::class_<RHFResult>(m, "RHFResult", "Result of an RHF calculation.")
        .def_prop_ro("converged", [](const RHFResult& r) { return r.converged; })
        .def_prop_ro("iterations", [](const RHFResult& r) { return r.iterations; })
        .def_prop_ro("energy", [](const RHFResult& r) { return r.energy; })
        .def_prop_ro("electronic_energy", [](const RHFResult& r) { return r.electronic_energy; })
        .def_prop_ro("residual", [](const RHFResult& r) { return r.residual; })
        .def("get_eigenvalues", [](const RHFResult& r) { return xt1d_to_numpy(r.eigenvalues); })
        .def("get_P", [](const RHFResult& r) { return xt2d_to_numpy(r.P); })
        .def("get_C", [](const RHFResult& r) { return xt2d_to_numpy(r.C); })
        .def_prop_ro("trace", [](const RHFResult& r) {
            nb::list trace;
            for (auto& rec : r.trace) {
                nb::dict d;
                d["iter"] = rec.iter;
                d["energy"] = rec.energy;
                d["residual"] = rec.residual;
                trace.append(d);
            }
            return trace;
        })
    ;

    // ── rhf_solve ──
    m.def("rhf_solve",
        [](const IntegralPack& ip, Real tol, int max_iter) {
            RHFResult result;
            { nb::gil_scoped_release release;
              result = rhf_solve(ip, tol, max_iter); }
            return result;
        },
        nb::arg("integrals"), nb::arg("tol") = Real(1e-8), nb::arg("max_iter") = 100,
        "Run Restricted Hartree-Fock on CPU (xtensor-blas)."
    );

    // ── UHFResult ──
    nb::class_<UHFResult>(m, "UHFResult", "Result of a UHF calculation.")
        .def_prop_ro("converged", [](const UHFResult& r) { return r.converged; })
        .def_prop_ro("iterations", [](const UHFResult& r) { return r.iterations; })
        .def_prop_ro("energy", [](const UHFResult& r) { return r.energy; })
        .def_prop_ro("electronic_energy", [](const UHFResult& r) { return r.electronic_energy; })
        .def_prop_ro("residual", [](const UHFResult& r) { return r.residual; })
        .def_prop_ro("n_alpha", [](const UHFResult& r) { return r.n_alpha; })
        .def_prop_ro("n_beta", [](const UHFResult& r) { return r.n_beta; })
        .def_prop_ro("spin_contamination", [](const UHFResult& r) { return r.spin_contamination; })
        .def("get_eigenvalues_alpha", [](const UHFResult& r) { return xt1d_to_numpy(r.eigenvalues_alpha); })
        .def("get_eigenvalues_beta", [](const UHFResult& r) { return xt1d_to_numpy(r.eigenvalues_beta); })
        .def("get_P_alpha", [](const UHFResult& r) { return xt2d_to_numpy(r.P_alpha); })
        .def("get_P_beta", [](const UHFResult& r) { return xt2d_to_numpy(r.P_beta); })
        .def_prop_ro("trace", [](const UHFResult& r) {
            nb::list trace;
            for (auto& rec : r.trace) {
                nb::dict d;
                d["iter"] = rec.iter;
                d["energy"] = rec.energy;
                d["residual"] = rec.residual;
                trace.append(d);
            }
            return trace;
        })
    ;

    // ── uhf_solve ──
    m.def("uhf_solve",
        [](const IntegralPack& ip, int n_alpha, int n_beta, Real tol, int max_iter) {
            UHFResult result;
            { nb::gil_scoped_release release;
              result = uhf_solve(ip, n_alpha, n_beta, tol, max_iter); }
            return result;
        },
        nb::arg("integrals"), nb::arg("n_alpha"), nb::arg("n_beta"),
        nb::arg("tol") = Real(1e-8), nb::arg("max_iter") = 100,
        "Run Unrestricted Hartree-Fock on CPU."
    );

    // ── XCType enum ──
    nb::enum_<XCType>(m, "XCType", "Exchange-correlation functional type.")
        .value("LSDA", XCType::LSDA)
        .value("PBE",  XCType::PBE)
    ;

    // ── UKSResult ──
    nb::class_<UKSResult>(m, "UKSResult", "Result of a UKS (DFT) calculation.")
        .def_prop_ro("converged", [](const UKSResult& r) { return r.converged; })
        .def_prop_ro("iterations", [](const UKSResult& r) { return r.iterations; })
        .def_prop_ro("energy", [](const UKSResult& r) { return r.energy; })
        .def_prop_ro("electronic_energy", [](const UKSResult& r) { return r.electronic_energy; })
        .def_prop_ro("one_electron_energy", [](const UKSResult& r) { return r.one_electron_energy; })
        .def_prop_ro("coulomb_energy", [](const UKSResult& r) { return r.coulomb_energy; })
        .def_prop_ro("xc_energy", [](const UKSResult& r) { return r.xc_energy; })
        .def_prop_ro("residual", [](const UKSResult& r) { return r.residual; })
        .def_prop_ro("spin_contamination", [](const UKSResult& r) { return r.spin_contamination; })
        .def("get_eigenvalues_alpha", [](const UKSResult& r) { return xt1d_to_numpy(r.eigenvalues_alpha); })
        .def("get_eigenvalues_beta", [](const UKSResult& r) { return xt1d_to_numpy(r.eigenvalues_beta); })
        .def("get_P_alpha", [](const UKSResult& r) { return xt2d_to_numpy(r.P_alpha); })
        .def("get_P_beta", [](const UKSResult& r) { return xt2d_to_numpy(r.P_beta); })
        .def_prop_ro("trace", [](const UKSResult& r) {
            nb::list trace;
            for (auto& rec : r.trace) {
                nb::dict d;
                d["iter"] = rec.iter;
                d["energy"] = rec.energy;
                d["residual"] = rec.residual;
                trace.append(d);
            }
            return trace;
        })
    ;

    // ── uks_solve ──
    m.def("uks_solve",
        [](const IntegralPack& ip, int n_alpha, int n_beta,
           XCType xc_type, int n_radial, int n_angular,
           Real tol, int max_iter) {
            UKSResult result;
            { nb::gil_scoped_release release;
              result = uks_solve(ip, n_alpha, n_beta, xc_type,
                                 n_radial, n_angular, tol, max_iter); }
            return result;
        },
        nb::arg("integrals"), nb::arg("n_alpha"), nb::arg("n_beta"),
        nb::arg("xc_type") = XCType::LSDA,
        nb::arg("n_radial") = 50, nb::arg("n_angular") = 110,
        nb::arg("tol") = Real(1e-8), nb::arg("max_iter") = 100,
        "Run Unrestricted Kohn-Sham DFT (LSDA or PBE) on CPU."
    );
}
