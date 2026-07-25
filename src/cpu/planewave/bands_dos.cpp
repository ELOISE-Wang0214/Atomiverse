/// @file cpu/planewave/bands_dos.cpp
/// @brief Band structure / DOS post-processing on frozen PW \(V_\mathrm{eff}\).

#include "atomiverse/cpu/planewave/bands_dos.h"

#include "atomiverse/core/log.h"
#include "atomiverse/core/vec3.h"
#include "atomiverse/cpu/planewave/internal/davidson.h"
#include "atomiverse/cpu/planewave/internal/kpoint.h"
#include "atomiverse/cpu/planewave/internal/pw_setup.h"
#include "atomiverse/cpu/planewave/internal/restricted_scf.h"
#include "atomiverse/cpu/planewave/qe_units.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ATV {
namespace cpu {
namespace pw {

namespace {

Vec3 arr_to_vec(const std::array<double, 3>& a) {
    return {static_cast<AccReal>(a[0]), static_cast<AccReal>(a[1]),
            static_cast<AccReal>(a[2])};
}

std::array<double, 3> vec_to_arr(const Vec3& v) {
    return {static_cast<double>(v[0]), static_cast<double>(v[1]),
            static_cast<double>(v[2])};
}

/// Diagonalize `n_bands` states at every k currently on `prob` (NSCF).
std::vector<std::vector<double>> nscf_evals_ry(internal::PWProblem& prob,
                                               int n_bands, double dav_tol) {
    using namespace internal;
    PWWorkspace& ws = prob.ws;
    const int n_kpts = static_cast<int>(prob.kg.kpoints.size());
    DavidsonSolver davidson(n_bands, AccReal(dav_tol), 50);
    davidson.prepare(ws);

    std::vector<std::vector<double>> evals(static_cast<std::size_t>(n_kpts));
    auto* saved_gvec = ws.gvec;
    const int n_g0 = ws.n_g;

    for (int ik = 0; ik < n_kpts; ++ik) {
        auto& gvk = prob.gv_k[static_cast<std::size_t>(ik)];
        ws.gvec = &gvk;
        ws.n_g = gvk.count();
        ws.T_G.resize({static_cast<std::size_t>(ws.n_g)});
        for (int ig = 0; ig < ws.n_g; ++ig)
            ws.T_G[ig] = gvk.ekin[ig];
        if (!prob.kb_proj_k.empty())
            ws.kb_proj = prob.kb_proj_k[static_cast<std::size_t>(ik)];
        if (!prob.atom_wfc_k.empty())
            ws.atom_wfc = prob.atom_wfc_k[static_cast<std::size_t>(ik)];

        davidson.step(ws);

        evals[static_cast<std::size_t>(ik)].resize(
            static_cast<std::size_t>(n_bands));
        for (int b = 0; b < n_bands; ++b) {
            // Davidson eigenvalues are Hartree; report Ry.
            const double e_ha = static_cast<double>(ws.eigenvalues[b]);
            evals[static_cast<std::size_t>(ik)][static_cast<std::size_t>(b)] =
                e_ha / qe_units::rydberg_to_hartree;
        }
    }

    ws.gvec = saved_gvec;
    ws.n_g = n_g0;
    return evals;
}

std::vector<Vec3> lerp_segment(const Vec3& a, const Vec3& b, int n_seg,
                               bool include_start) {
    // n_seg samples along the segment; include_start controls whether `a`
    // is emitted (false when chaining segments that share the join).
    std::vector<Vec3> pts;
    if (n_seg < 2)
        n_seg = 2;
    const int start = include_start ? 0 : 1;
    for (int i = start; i < n_seg; ++i) {
        const AccReal t = AccReal(i) / AccReal(n_seg - 1);
        pts.push_back({a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]),
                       a[2] + t * (b[2] - a[2])});
    }
    return pts;
}

} // namespace

double gap_from_spectrum(const std::vector<std::vector<double>>& evals,
                         int n_occ, double& e_vbm, double& e_cbm, int& vbm_k,
                         int& cbm_k, int& vbm_b, int& cbm_b) {
    if (evals.empty())
        throw std::invalid_argument("gap_from_spectrum: empty spectrum");
    if (n_occ < 1)
        throw std::invalid_argument("gap_from_spectrum: n_occ < 1");
    const int n_bands = static_cast<int>(evals[0].size());
    if (n_bands < n_occ + 1)
        throw std::invalid_argument("gap_from_spectrum: need at least one "
                                    "empty band (n_bands > n_occ)");

    e_vbm = -1e300;
    e_cbm = 1e300;
    vbm_k = cbm_k = 0;
    vbm_b = n_occ - 1;
    cbm_b = n_occ;

    for (std::size_t ik = 0; ik < evals.size(); ++ik) {
        if (static_cast<int>(evals[ik].size()) < n_occ + 1)
            throw std::invalid_argument(
                "gap_from_spectrum: inconsistent band count");
        const double v = evals[ik][static_cast<std::size_t>(n_occ - 1)];
        const double c = evals[ik][static_cast<std::size_t>(n_occ)];
        if (v > e_vbm) {
            e_vbm = v;
            vbm_k = static_cast<int>(ik);
            vbm_b = n_occ - 1;
        }
        if (c < e_cbm) {
            e_cbm = c;
            cbm_k = static_cast<int>(ik);
            cbm_b = n_occ;
        }
    }
    return e_cbm - e_vbm;
}

void gaussian_dos(const std::vector<double>& evals,
                  const std::vector<double>& weights, double e_min,
                  double e_max, int n_bins, double sigma,
                  std::vector<double>& energy, std::vector<double>& dos) {
    if (evals.size() != weights.size())
        throw std::invalid_argument(
            "gaussian_dos: evals/weights size mismatch");
    if (n_bins < 2)
        throw std::invalid_argument("gaussian_dos: n_bins < 2");
    if (!(sigma > 0.0))
        throw std::invalid_argument("gaussian_dos: sigma must be > 0");
    if (!(e_max > e_min))
        throw std::invalid_argument("gaussian_dos: need e_max > e_min");

    energy.resize(static_cast<std::size_t>(n_bins));
    dos.assign(static_cast<std::size_t>(n_bins), 0.0);
    const double de = (e_max - e_min) / static_cast<double>(n_bins - 1);
    const double pref = 2.0 / (std::sqrt(2.0 * 3.14159265358979323846) * sigma);
    const double inv_2s2 = 1.0 / (2.0 * sigma * sigma);

    for (int i = 0; i < n_bins; ++i)
        energy[static_cast<std::size_t>(i)] =
            e_min + de * static_cast<double>(i);

    for (std::size_t n = 0; n < evals.size(); ++n) {
        const double e = evals[n];
        const double w = weights[n];
        for (int i = 0; i < n_bins; ++i) {
            const double x = energy[static_cast<std::size_t>(i)] - e;
            dos[static_cast<std::size_t>(i)] +=
                w * pref * std::exp(-x * x * inv_2s2);
        }
    }
}

void Bands::set_diamond_path(int n_seg) {
    if (n_seg < 2)
        n_seg = 2;
    // Standard FCC primitive reciprocal path (crystal coordinates).
    const Vec3 G{0, 0, 0};
    const Vec3 X{0.5, 0.0, 0.5};
    const Vec3 W{0.5, 0.25, 0.75};
    const Vec3 L{0.5, 0.5, 0.5};
    const Vec3 K{0.375, 0.375, 0.75};

    std::vector<Vec3> path;
    auto append = [&](const Vec3& a, const Vec3& b, bool include_start) {
        auto seg = lerp_segment(a, b, n_seg, include_start);
        path.insert(path.end(), seg.begin(), seg.end());
    };
    append(G, X, true);
    append(X, W, false);
    append(W, L, false);
    append(L, G, false);
    append(G, K, false);

    kfrac.clear();
    kfrac.reserve(path.size());
    for (const auto& p : path)
        kfrac.push_back(vec_to_arr(p));
}

BandsResult Bands::compute(const PW& method, const cpu::Frame& frame) const {
    if (method.spin_polarized)
        throw std::invalid_argument("Bands::compute: spin-polarized bands not "
                                    "supported in this foothold");
    if (method.smearing != Smearing::none)
        throw std::invalid_argument("Bands::compute: metallic smearing not "
                                    "supported for bands foothold");
    if (n_bands < 2)
        throw std::invalid_argument("Bands::compute: n_bands must be >= 2");

    std::vector<std::array<double, 3>> path = kfrac;
    if (path.empty()) {
        Bands tmp;
        tmp.set_diamond_path(4);
        path = tmp.kfrac;
    }

    internal::PWProblem prob;
    const PWResult scf = internal::run_restricted_scf(frame, method, prob);
    if (n_bands < prob.n_bands)
        throw std::invalid_argument(
            "Bands::compute: n_bands must be >= occupied band count");

    std::vector<Vec3> kvec;
    kvec.reserve(path.size());
    for (const auto& k : path)
        kvec.push_back(arr_to_vec(k));
    rebuild_k_tables(prob, frame, method.upf_path, method.ecut_ry, kvec);

    BandsResult out;
    out.kfrac = path;
    out.scf_converged = scf.converged;
    out.energy_ry = scf.energy;
    out.eigenvalues_ry = nscf_evals_ry(prob, n_bands, method.tol * 0.1);

    const int n_occ = prob.n_bands; // n_e / 2
    out.gap_ry =
        gap_from_spectrum(out.eigenvalues_ry, n_occ, out.e_vbm_ry, out.e_cbm_ry,
                          out.vbm_k, out.cbm_k, out.vbm_band, out.cbm_band);

    return out;
}

DosResult Dos::compute(const PW& method, const cpu::Frame& frame) const {
    if (method.spin_polarized)
        throw std::invalid_argument(
            "Dos::compute: spin-polarized DOS not supported in this foothold");
    if (method.smearing != Smearing::none)
        throw std::invalid_argument(
            "Dos::compute: metallic smearing not supported for DOS foothold");
    if (n_bands < 2)
        throw std::invalid_argument("Dos::compute: n_bands must be >= 2");
    if (k_mesh < 1)
        throw std::invalid_argument("Dos::compute: k_mesh must be >= 1");
    if (!(sigma_ry > 0.0))
        throw std::invalid_argument("Dos::compute: sigma_ry must be > 0");
    if (n_bins < 2)
        throw std::invalid_argument("Dos::compute: n_bins must be >= 2");

    // SCF density on method.k_mesh (can be cheaper than DOS mesh).
    internal::PWProblem prob;
    const PWResult scf = internal::run_restricted_scf(frame, method, prob);
    if (n_bands < prob.n_bands)
        throw std::invalid_argument(
            "Dos::compute: n_bands must be >= occupied band count");

    // Replace k-mesh with denser MP grid for the DOS eigenvalues.
    auto kg = internal::generate_mp_mesh(k_mesh, k_mesh, k_mesh);
    rebuild_k_tables(prob, frame, method.upf_path, method.ecut_ry, kg.kpoints);

    auto evals = nscf_evals_ry(prob, n_bands, method.tol * 0.1);
    const int n_kpts = static_cast<int>(evals.size());
    const int n_occ = prob.n_bands;

    DosResult out;
    out.scf_converged = scf.converged;
    int vk = 0, ck = 0, vb = 0, cb = 0;
    out.gap_ry = gap_from_spectrum(evals, n_occ, out.e_vbm_ry, out.e_cbm_ry, vk,
                                   ck, vb, cb);

    std::vector<double> flat_e, flat_w;
    flat_e.reserve(static_cast<std::size_t>(n_kpts * n_bands));
    flat_w.reserve(static_cast<std::size_t>(n_kpts * n_bands));
    for (int ik = 0; ik < n_kpts; ++ik) {
        const double w =
            static_cast<double>(kg.weights[static_cast<std::size_t>(ik)]);
        for (int b = 0; b < n_bands; ++b) {
            flat_e.push_back(evals[static_cast<std::size_t>(ik)]
                                  [static_cast<std::size_t>(b)]);
            flat_w.push_back(w);
        }
    }

    double emin = e_min_ry;
    double emax = e_max_ry;
    if (!(emax > emin)) {
        emin = *std::min_element(flat_e.begin(), flat_e.end()) - 4.0 * sigma_ry;
        emax = *std::max_element(flat_e.begin(), flat_e.end()) + 4.0 * sigma_ry;
    }

    gaussian_dos(flat_e, flat_w, emin, emax, n_bins, sigma_ry, out.energy_ry,
                 out.dos);

    // Trapezoid integrate DOS up to VBM (closed-shell factor already in D).
    const double de =
        (out.energy_ry.back() - out.energy_ry.front())
        / static_cast<double>(static_cast<int>(out.energy_ry.size()) - 1);
    double integ = 0;
    for (std::size_t i = 0; i + 1 < out.energy_ry.size(); ++i) {
        if (out.energy_ry[i] > out.e_vbm_ry)
            break;
        const double e1 = out.energy_ry[i];
        const double e2 = out.energy_ry[i + 1];
        const double d1 = out.dos[i];
        const double d2 = out.dos[i + 1];
        if (e2 <= out.e_vbm_ry) {
            integ += 0.5 * (d1 + d2) * de;
        } else {
            // Partial last bin up to VBM.
            const double t = (out.e_vbm_ry - e1) / (e2 - e1);
            const double d_v = d1 + t * (d2 - d1);
            integ += 0.5 * (d1 + d_v) * (out.e_vbm_ry - e1);
            break;
        }
    }
    out.n_elec_to_vbm = integ;
    return out;
}

} // namespace pw
} // namespace cpu
} // namespace ATV
