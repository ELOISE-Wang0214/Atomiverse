/// @file regressions/pw/pw-uks-fe-smearing/reg_pw-uks-fe-smearing.cpp
/// @brief END-TO-END: BCC Fe collinear UKS + Fermi–Dirac metallic smearing.
///
/// See REFERENCE.md. Public API only. Requires ATV_BUILD_PW.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/planewave/pw_solver.h"

using namespace ATV;
using namespace ATV::cpu::pw;

namespace {

/// Primitive BCC cell from cubic lattice constant a [Å], one Fe at origin.
cpu::Frame make_fe_bcc_primitive(double a_ang) {
    cpu::Frame f(1);
    f.span(kAtom_type)[0] = 26;
    f.span(kAtom_mass)[0] = 55.845;
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    const double h = a_ang / 2.0;
    Mat3 cell;
    // Primitive BCC: (-h,h,h), (h,-h,h), (h,h,-h)
    cell(0, 0) = -h;
    cell(0, 1) = h;
    cell(0, 2) = h;
    cell(1, 0) = h;
    cell(1, 1) = -h;
    cell(1, 2) = h;
    cell(2, 0) = h;
    cell(2, 1) = h;
    cell(2, 2) = -h;
    f.box.set_cell(cell);
    return f;
}

std::string fe_upf() {
#ifdef ATV_SOURCE_DIR
    {
        std::string p =
            std::string(ATV_SOURCE_DIR) + "/tests/data/pw/Fe_ONCV_PBE-1.2.upf";
        if (std::FILE* f = std::fopen(p.c_str(), "r")) {
            std::fclose(f);
            return p;
        }
    }
#endif
    if (const char* env = std::getenv("ATV_SOURCE_DIR")) {
        std::string p = std::string(env) + "/tests/data/pw/Fe_ONCV_PBE-1.2.upf";
        if (std::FILE* f = std::fopen(p.c_str(), "r")) {
            std::fclose(f);
            return p;
        }
    }
    const char* cands[] = {
        "tests/data/pw/Fe_ONCV_PBE-1.2.upf",
        "../tests/data/pw/Fe_ONCV_PBE-1.2.upf",
        "../../tests/data/pw/Fe_ONCV_PBE-1.2.upf",
        "../../../tests/data/pw/Fe_ONCV_PBE-1.2.upf",
    };
    for (const char* c : cands) {
        if (std::FILE* f = std::fopen(c, "r")) {
            std::fclose(f);
            return c;
        }
    }
    return "tests/data/pw/Fe_ONCV_PBE-1.2.upf";
}

} // namespace

int main() {
    // Offline pins — see REFERENCE.md
    constexpr double want_e_ry = -205.8365363098; // free energy [Ry]
    constexpr double tol_e_ry = 0.05;             // Ry — modest ecut / k-mesh
    constexpr double want_ne = 16.0;              // SG15 Fe Z_val
    constexpr double tol_ne = 1e-4;
    constexpr double mag_min = 0.1; // ferro seed; self-consistent M ~ 1

    const std::string upf = fe_upf();
    if (std::FILE* f = std::fopen(upf.c_str(), "r")) {
        std::fclose(f);
    } else {
        std::printf("FAIL: missing UPF %s\n", upf.c_str());
        return 1;
    }

    auto frame = make_fe_bcc_primitive(/*a_ang=*/2.87);

    PW pw;
    pw.upf_path = upf;
    pw.ecut_ry = 12.0; // modest algorithm pin (≤ 60 s)
    pw.k_mesh = 2;
    pw.tol = 2e-3;
    pw.max_iter = 30;
    pw.mix_beta = 0.3;
    pw.pulay_depth = 5;
    pw.xc = cpu::scf::XCType::PBE;
    pw.spin_polarized = true;
    pw.multiplicity = 5; // seed M≈4; free under smear
    pw.smearing = Smearing::fermi_dirac;
    pw.degauss = 0.02; // Ry
    pw.extra_bands = 4;

    pw.build(frame);
    const auto& r = pw.run(frame);

    std::printf("E=%.10f Ry (%.10f Ha)  mag=%.6f  e_fermi=%.6f Ry  "
                "n_occ=%.8f  e_ent=%.6f Ry  conv=%d iters=%d\n",
                r.energy, r.energy_ha, r.magnetization, r.e_fermi, r.n_elec_occ,
                r.e_entropy, (int)r.converged, r.iterations);

    bool ok = true;
    if (!std::isfinite(r.energy) || !std::isfinite(r.energy_ha)
        || !std::isfinite(r.e_fermi) || !std::isfinite(r.magnetization)
        || !std::isfinite(r.e_entropy)) {
        std::printf("FAIL: non-finite energy/fermi/mag/entropy\n");
        ok = false;
    }
    if (!r.converged) {
        std::printf("FAIL: SCF not converged\n");
        ok = false;
    }
    if (std::fabs(r.n_elec_occ - want_ne) > tol_ne) {
        std::printf("FAIL: n_elec_occ=%.8f want=%.1f ± %g\n", r.n_elec_occ,
                    want_ne, tol_ne);
        ok = false;
    }
    if (!(r.magnetization > mag_min)) {
        std::printf("FAIL: magnetization=%.6f expected > %g (ferro seed)\n",
                    r.magnetization, mag_min);
        ok = false;
    }
    if (!(std::fabs(r.e_fermi) < 1e3)) {
        std::printf("FAIL: absurd e_fermi=%.6f Ry\n", r.e_fermi);
        ok = false;
    }
    if (std::fabs(r.energy - want_e_ry) > tol_e_ry) {
        std::printf("FAIL: E=%.10f Ry  pin=%.10f ± %g\n", r.energy, want_e_ry,
                    tol_e_ry);
        ok = false;
    }

    if (!ok)
        return 1;
    std::printf("PASS: pw-uks-fe-smearing (BCC Fe UKS + FD smear)\n");
    return 0;
}
