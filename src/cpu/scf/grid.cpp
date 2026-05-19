/// @file cpu/scf/grid.cpp
/// @brief Molecular integration grid: Lebedev angular + Euler-Maclaurin radial + Becke partition.
///
/// Lebedev angular quadrature points tabulated for 6-point and 110-point grids.
/// Higher-order grids (194, 302, 590) can be added later.
///
/// Radial quadrature uses Euler-Maclaurin mapping with Bragg-Slater atomic radii.
/// Becke partitioning uses the k(μ) = 1.5μ - 0.5μ³ step function applied 3 times.

#include "atomiverse/cpu/scf/grid.h"
#include "atomiverse/cpu/scf/lebedev_data.h"
#include "atomiverse/core/constants.h"

#include <cmath>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <string>

namespace ATV {
namespace cpu {
namespace scf {

static constexpr double kPi = 3.14159265358979323846;

// ── Bragg-Slater atomic radii (Bohr) ────────────────────────────────────

/// @brief Bragg-Slater atomic radius for element Z, returned in Bohr.
///
/// Radii for Z = 1..86 (H through Rn).  Most values are the empirical
/// crystal radii of Slater:
///   J.C. Slater, "Atomic Radii in Crystals", JCP 41, 3199 (1964).
///   DOI: 10.1063/1.1725697
///
/// Noble gases (He, Ne, Ar, Kr, Xe, Rn) use calculated radii from:
///   E. Clementi, D.L. Raimondi, W.P. Reinhardt, JCP 47, 1300 (1967).
///
/// Astatine (At, Z = 85) uses the covalent radius from:
///   B. Cordero et al., Dalton Trans. 2008, 2832.
///   DOI: 10.1039/b801115j
///
/// The table is stored in Bohr (Angstrom values * kAngstromToBohr).
/// Throws std::invalid_argument for Z < 1 or Z > 86.
///
/// @param Z Atomic number (1 <= Z <= 86).
/// @return Bragg-Slater radius in Bohr.
double bragg_slater_radius(int Z) {
    // Angstrom values -- index 0 is unused (Z starts at 1).
    //
    // Period 1:  H    He
    // Period 2:  Li   Be   B    C    N    O    F    Ne
    // Period 3:  Na   Mg   Al   Si   P    S    Cl   Ar
    // Period 4:  K    Ca   Sc   Ti   V    Cr   Mn   Fe   Co   Ni   Cu   Zn
    //            Ga   Ge   As   Se   Br   Kr
    // Period 5:  Rb   Sr   Y    Zr   Nb   Mo   Tc   Ru   Rh   Pd   Ag   Cd
    //            In   Sn   Sb   Te   I    Xe
    // Period 6:  Cs   Ba   La   Ce   Pr   Nd   Pm   Sm   Eu   Gd   Tb   Dy
    //            Ho   Er   Tm   Yb   Lu   Hf   Ta   W    Re   Os   Ir   Pt
    //            Au   Hg   Tl   Pb   Bi   Po   At   Rn
    static constexpr double kA2B = kAngstromToBohr;
    static constexpr int kMaxZ = 86;
    // clang-format off
    static constexpr double radii[kMaxZ + 1] = {
        0.00 * kA2B,  //  0  placeholder (unused)
        // ── Period 1 ─────────────────────────────────────────────
        0.25 * kA2B,  //  1  H                 [Slater 1964]
        0.31 * kA2B,  //  2  He                [Clementi 1967]
        // ── Period 2 ─────────────────────────────────────────────
        1.45 * kA2B,  //  3  Li                [Slater 1964]
        1.05 * kA2B,  //  4  Be                [Slater 1964]
        0.85 * kA2B,  //  5  B                 [Slater 1964]
        0.70 * kA2B,  //  6  C                 [Slater 1964]
        0.65 * kA2B,  //  7  N                 [Slater 1964]
        0.60 * kA2B,  //  8  O                 [Slater 1964]
        0.50 * kA2B,  //  9  F                 [Slater 1964]
        0.38 * kA2B,  // 10  Ne                [Clementi 1967]
        // ── Period 3 ─────────────────────────────────────────────
        1.80 * kA2B,  // 11  Na                [Slater 1964]
        1.50 * kA2B,  // 12  Mg                [Slater 1964]
        1.25 * kA2B,  // 13  Al                [Slater 1964]
        1.10 * kA2B,  // 14  Si                [Slater 1964]
        1.00 * kA2B,  // 15  P                 [Slater 1964]
        1.00 * kA2B,  // 16  S                 [Slater 1964]
        1.00 * kA2B,  // 17  Cl                [Slater 1964]
        0.71 * kA2B,  // 18  Ar                [Clementi 1967]
        // ── Period 4 ─────────────────────────────────────────────
        2.20 * kA2B,  // 19  K                 [Slater 1964]
        1.80 * kA2B,  // 20  Ca                [Slater 1964]
        1.60 * kA2B,  // 21  Sc                [Slater 1964]
        1.40 * kA2B,  // 22  Ti                [Slater 1964]
        1.35 * kA2B,  // 23  V                 [Slater 1964]
        1.40 * kA2B,  // 24  Cr                [Slater 1964]
        1.40 * kA2B,  // 25  Mn                [Slater 1964]
        1.40 * kA2B,  // 26  Fe                [Slater 1964]
        1.35 * kA2B,  // 27  Co                [Slater 1964]
        1.35 * kA2B,  // 28  Ni                [Slater 1964]
        1.35 * kA2B,  // 29  Cu                [Slater 1964]
        1.35 * kA2B,  // 30  Zn                [Slater 1964]
        1.30 * kA2B,  // 31  Ga                [Slater 1964]
        1.25 * kA2B,  // 32  Ge                [Slater 1964]
        1.15 * kA2B,  // 33  As                [Slater 1964]
        1.15 * kA2B,  // 34  Se                [Slater 1964]
        1.15 * kA2B,  // 35  Br                [Slater 1964]
        0.88 * kA2B,  // 36  Kr                [Clementi 1967]
        // ── Period 5 ─────────────────────────────────────────────
        2.35 * kA2B,  // 37  Rb                [Slater 1964]
        2.00 * kA2B,  // 38  Sr                [Slater 1964]
        1.80 * kA2B,  // 39  Y                 [Slater 1964]
        1.55 * kA2B,  // 40  Zr                [Slater 1964]
        1.45 * kA2B,  // 41  Nb                [Slater 1964]
        1.45 * kA2B,  // 42  Mo                [Slater 1964]
        1.35 * kA2B,  // 43  Tc                [Slater 1964]
        1.30 * kA2B,  // 44  Ru                [Slater 1964]
        1.35 * kA2B,  // 45  Rh                [Slater 1964]
        1.40 * kA2B,  // 46  Pd                [Slater 1964]
        1.60 * kA2B,  // 47  Ag                [Slater 1964]
        1.55 * kA2B,  // 48  Cd                [Slater 1964]
        1.55 * kA2B,  // 49  In                [Slater 1964]
        1.45 * kA2B,  // 50  Sn                [Slater 1964]
        1.45 * kA2B,  // 51  Sb                [Slater 1964]
        1.40 * kA2B,  // 52  Te                [Slater 1964]
        1.40 * kA2B,  // 53  I                 [Slater 1964]
        1.08 * kA2B,  // 54  Xe                [Clementi 1967]
        // ── Period 6 ─────────────────────────────────────────────
        2.60 * kA2B,  // 55  Cs                [Slater 1964]
        2.15 * kA2B,  // 56  Ba                [Slater 1964]
        1.95 * kA2B,  // 57  La                [Slater 1964]
        1.85 * kA2B,  // 58  Ce                [Slater 1964]
        1.85 * kA2B,  // 59  Pr                [Slater 1964]
        1.85 * kA2B,  // 60  Nd                [Slater 1964]
        1.85 * kA2B,  // 61  Pm                [Slater 1964]
        1.85 * kA2B,  // 62  Sm                [Slater 1964]
        1.85 * kA2B,  // 63  Eu                [Slater 1964]
        1.80 * kA2B,  // 64  Gd                [Slater 1964]
        1.75 * kA2B,  // 65  Tb                [Slater 1964]
        1.75 * kA2B,  // 66  Dy                [Slater 1964]
        1.75 * kA2B,  // 67  Ho                [Slater 1964]
        1.75 * kA2B,  // 68  Er                [Slater 1964]
        1.75 * kA2B,  // 69  Tm                [Slater 1964]
        1.75 * kA2B,  // 70  Yb                [Slater 1964]
        1.75 * kA2B,  // 71  Lu                [Slater 1964]
        1.55 * kA2B,  // 72  Hf                [Slater 1964]
        1.45 * kA2B,  // 73  Ta                [Slater 1964]
        1.35 * kA2B,  // 74  W                 [Slater 1964]
        1.35 * kA2B,  // 75  Re                [Slater 1964]
        1.30 * kA2B,  // 76  Os                [Slater 1964]
        1.35 * kA2B,  // 77  Ir                [Slater 1964]
        1.35 * kA2B,  // 78  Pt                [Slater 1964]
        1.35 * kA2B,  // 79  Au                [Slater 1964]
        1.50 * kA2B,  // 80  Hg                [Slater 1964]
        1.90 * kA2B,  // 81  Tl                [Slater 1964]
        1.80 * kA2B,  // 82  Pb                [Slater 1964]
        1.60 * kA2B,  // 83  Bi                [Slater 1964]
        1.90 * kA2B,  // 84  Po                [Slater 1964]
        1.50 * kA2B,  // 85  At                [Cordero 2008]
        1.20 * kA2B,  // 86  Rn                [Clementi 1967]
    };
    // clang-format on
    if (Z < 1 || Z > kMaxZ) {
        throw std::invalid_argument(
            "bragg_slater_radius: unsupported element Z="
            + std::to_string(Z) + ". Supported range: Z=1..86");
    }
    return radii[Z];
}

// ── Lebedev angular quadrature ──────────────────────────────────────────
//
// All angular grids come from the tabulated Lebedev-Laikov data in
// lebedev_data.h (orders 6, 50, 110, 194, 302, 590). The earlier
// in-source octahedral-generator code and the Gauss-Legendre x uniform-phi
// product grid that was used as a placeholder for "Lebedev-110" have
// been removed -- they produced wrong weight distributions for spherical
// harmonics. See `build_molecular_grid` below for the dispatch.

// ── Radial quadrature ───────────────────────────────────────────────────

// Euler-Maclaurin / Mura-Knowles-type log mapping.
// Maps i ∈ [1, n_radial] to r ∈ [0, ∞) with Bragg-Slater scaling.
//
// Treutler-Ahlrichs mapping (widely used, good for all-electron):
//   r_i = R_bragg * (1 + cos(θ_i))^alpha / (1 - cos(θ_i))
//
// Simpler approach: Murray-Handy-Laming (Mol. Phys. 78, 997, 1993):
//   r_i = R_bragg * (x_i / (1 - x_i))^2
//   x_i = i / (n + 1)
//   w_i = 2 R_bragg * x_i^2 / (1 - x_i)^3 * r_i^2 / (n+1)
//
// But simplest is the Becke mapping:
//   x_i = cos(π(n_radial - i + 0.5) / n_radial)   Chebyshev points
//   r_i = R * (1 + x_i) / (1 - x_i)
//   w_i = 2R / (1 - x_i)^2 * π/n_radial * sin(arccos(x_i)) * r_i^2

// RadialPoint is declared in grid.h (public API).
//
// Only Mura-Knowles is used; the earlier Becke (Chebyshev-mapped) radial
// quadrature was removed -- it left an ad-hoc radius scaling that did not
// match PySCF's reference grid and contributed to the ~mHa UKS-vs-PySCF
// offset documented in the libxc-xc-fix spec.

// ── Mura-Knowles radial quadrature ─────────────────────────────────────
// Reference: Mura & Knowles, JCP 104, 9848 (1996), equation (5).
// Uses the log3 mapping: r = -R * ln(1 - x^3), x = (i+0.5)/n.
// Matches PySCF's default radial grid (dft.gen_grid.RADII_ADJUST = none).

std::vector<RadialPoint> radial_grid_mura_knowles(int n, double R_bragg) {
    std::vector<RadialPoint> pts(n);
    for (int i = 0; i < n; ++i) {
        double x = (i + 0.5) / n;
        double x3 = x * x * x;
        double r = -R_bragg * std::log(1.0 - x3);

        // dr/dx = R * 3x^2 / (1 - x^3)
        double dr_dx = R_bragg * 3.0 * x * x / (1.0 - x3);
        // Uniform weight: 1/n
        double w_quad = 1.0 / n;

        // Total weight = r^2 * dr/dx * w_quad (4pi from angular)
        pts[i].r = r;
        pts[i].w = r * r * dr_dx * w_quad;
    }
    return pts;
}

// ── Becke partitioning ──────────────────────────────────────────────────

// Becke step function: k(μ) = 1.5μ - 0.5μ³
static double becke_step(double mu) {
    return 1.5 * mu - 0.5 * mu * mu * mu;
}

// Becke partition weight for atom A at point r.
// Uses Becke's fuzzy-cell scheme with 3 iterations of the step function.
// With Bragg-Slater heteronuclear adjustment.
static void becke_partition_weights(const native::Atom* atoms, int natm,
                                    double rx, double ry, double rz,
                                    double* weights) {
    if (natm == 1) {
        weights[0] = 1.0;
        return;
    }

    // Compute distances from grid point to each atom
    std::vector<double> dist(natm);
    for (int a = 0; a < natm; ++a) {
        double dx = rx - atoms[a].pos[0];
        double dy = ry - atoms[a].pos[1];
        double dz = rz - atoms[a].pos[2];
        dist[a] = std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    // Compute inter-atomic distances
    // R_AB = |R_A - R_B|
    // μ_AB = (r_A - r_B) / R_AB  (elliptical coordinate)

    // Initialize raw Becke weights to 1
    for (int a = 0; a < natm; ++a) weights[a] = 1.0;

    for (int a = 0; a < natm; ++a) {
        for (int b = a + 1; b < natm; ++b) {
            double dx = atoms[a].pos[0] - atoms[b].pos[0];
            double dy = atoms[a].pos[1] - atoms[b].pos[1];
            double dz = atoms[a].pos[2] - atoms[b].pos[2];
            double R_AB = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (R_AB < 1e-14) continue;

            double mu = (dist[a] - dist[b]) / R_AB;

            // Heteronuclear adjustment (Becke, eq. 21)
            double chi = bragg_slater_radius(atoms[a].Z) / bragg_slater_radius(atoms[b].Z);
            double u_AB = (chi - 1.0) / (chi + 1.0);
            double a_AB = u_AB / (u_AB * u_AB - 1.0);
            // Clamp a_AB to [-0.5, 0.5]
            a_AB = std::max(-0.5, std::min(0.5, a_AB));
            mu += a_AB * (1.0 - mu * mu);

            // Apply step function 3 times
            double s = mu;
            s = becke_step(s);
            s = becke_step(s);
            s = becke_step(s);

            // s_AB ∈ [-1, 1], convert to cell function
            // P_A *= 0.5*(1 - s_AB),  P_B *= 0.5*(1 + s_AB)
            weights[a] *= 0.5 * (1.0 - s);
            weights[b] *= 0.5 * (1.0 + s);
        }
    }

    // Normalize
    double total = 0;
    for (int a = 0; a < natm; ++a) total += weights[a];
    if (total > 1e-30) {
        for (int a = 0; a < natm; ++a) weights[a] /= total;
    }
}

// ── Public API ──────────────────────────────────────────────────────────

MolecularGrid build_molecular_grid(const native::Atom* atoms, int natm,
                                   int n_radial, int n_angular) {
    MolecularGrid grid;
    grid.n_radial = n_radial;
    grid.n_angular = n_angular;

    // Get angular grid from tabulated Lebedev-Laikov data.
    // Select the largest available grid that does not exceed n_angular.
    const lebedev::Point* leb_pts = nullptr;
    int actual_n_angular = 0;
    if (n_angular >= 590) {
        leb_pts = lebedev::grid_590; actual_n_angular = lebedev::n_590;
    } else if (n_angular >= 302) {
        leb_pts = lebedev::grid_302; actual_n_angular = lebedev::n_302;
    } else if (n_angular >= 194) {
        leb_pts = lebedev::grid_194; actual_n_angular = lebedev::n_194;
    } else if (n_angular >= 110) {
        leb_pts = lebedev::grid_110; actual_n_angular = lebedev::n_110;
    } else if (n_angular >= 50) {
        leb_pts = lebedev::grid_50;  actual_n_angular = lebedev::n_50;
    } else {
        leb_pts = lebedev::grid_6;   actual_n_angular = lebedev::n_6;
    }

    // Pre-allocate
    grid.points.reserve(natm * n_radial * actual_n_angular);

    std::vector<double> becke_w(natm);

    for (int a = 0; a < natm; ++a) {
        double R_bragg = bragg_slater_radius(atoms[a].Z);
        auto rad_pts = radial_grid_mura_knowles(n_radial, R_bragg);

        for (int ir = 0; ir < n_radial; ++ir) {
            double r = rad_pts[ir].r;
            double w_rad = rad_pts[ir].w;

            for (int ia = 0; ia < actual_n_angular; ++ia) {
                double px = atoms[a].pos[0] + r * leb_pts[ia].x;
                double py = atoms[a].pos[1] + r * leb_pts[ia].y;
                double pz = atoms[a].pos[2] + r * leb_pts[ia].z;

                // Becke partition weight
                becke_partition_weights(atoms, natm, px, py, pz, becke_w.data());

                // Lebedev weights sum to 1.0, multiply by 4π for spherical integral
                double w = 4.0 * kPi * w_rad * leb_pts[ia].w * becke_w[a];

                if (w > 1e-30) {
                    grid.points.push_back({
                        static_cast<double>(px),
                        static_cast<double>(py),
                        static_cast<double>(pz),
                        static_cast<double>(w)
                    });
                }
            }
        }
    }

    return grid;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
