/// @file cpu/pw/internal/ewald.cpp
/// @brief 3D Ewald summation implementation.

#include "atomiverse/cpu/pw/internal/ewald.h"
#include "atomiverse/core/units.h"
#include <cmath>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

EwaldSummation::EwaldSummation(const std::vector<Vec3>& pos_ang,
                               const std::vector<AccReal>& charges,
                               const ::ATV::Box& cell, AccReal ecut) {
    // Convert to Bohr internally
    AccReal a2b = AccReal(1.8897259886);
    Box cell_bohr = cell;
    cell_bohr.cell = cell_bohr.cell * a2b;

    std::vector<Vec3> pos(charges.size());
    for (size_t i = 0; i < charges.size(); ++i)
        pos[i] = {pos_ang[i][0] * a2b, pos_ang[i][1] * a2b,
                  pos_ang[i][2] * a2b};

    AccReal vol = cell_bohr.volume();
    if (vol <= 0)
        return;

    // Ewald alpha: auto-select, calibrated to match QE ewlad.f90
    AccReal rcut = cell_bohr.lengths().norm() * AccReal(0.5);
    alpha_ = std::sqrt(-std::log(AccReal(2e-17))) / rcut;
    AccReal alpha = alpha_;
    AccReal alpha_sq_4 = AccReal(4) * alpha * alpha;

    int n = static_cast<int>(charges.size());

    // Real space
    e_real_ = 0;
    AccReal r_cut =
        AccReal(10); // Angstrom in Bohr units? Let's use cell-based cutoff
    // actually: r_cut such that erfc(alpha*r) < 1e-12, roughly r_cut ~ 5/alpha
    r_cut = AccReal(5) / alpha;
    auto lens = cell_bohr.lengths();
    int ncell_r = static_cast<int>(
                      std::ceil(r_cut / std::min({lens[0], lens[1], lens[2]})))
                  + 1;

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            AccReal zi = charges[i], zj = charges[j];
            AccReal fac = (i == j) ? AccReal(0.5) : AccReal(1);
            for (int hx = -ncell_r; hx <= ncell_r; ++hx)
                for (int hy = -ncell_r; hy <= ncell_r; ++hy)
                    for (int hz = -ncell_r; hz <= ncell_r; ++hz) {
                        Vec3 dr = {pos[i][0] - pos[j][0]
                                       + AccReal(hx) * cell_bohr.cell(0, 0)
                                       + AccReal(hy) * cell_bohr.cell(0, 1)
                                       + AccReal(hz) * cell_bohr.cell(0, 2),
                                   pos[i][1] - pos[j][1]
                                       + AccReal(hx) * cell_bohr.cell(1, 0)
                                       + AccReal(hy) * cell_bohr.cell(1, 1)
                                       + AccReal(hz) * cell_bohr.cell(1, 2),
                                   pos[i][2] - pos[j][2]
                                       + AccReal(hx) * cell_bohr.cell(2, 0)
                                       + AccReal(hy) * cell_bohr.cell(2, 1)
                                       + AccReal(hz) * cell_bohr.cell(2, 2)};
                        AccReal r = dr.norm();
                        if (r > 0 && r < r_cut) {
                            e_real_ += fac * zi * zj * std::erfc(alpha * r) / r;
                        }
                        if (r < AccReal(1e-12) && i == j && hx == 0 && hy == 0
                            && hz == 0) {
                            // Self term handled separately
                        }
                    }
        }
    }

    // Reciprocal space
    e_recip_ = 0;
    AccReal g_cut = std::sqrt(AccReal(8) * ecut); // 4x E_cut for Ewald recip
    auto recip = cell_bohr.recip_cell();
    Vec3 b0{recip(0, 0), recip(0, 1), recip(0, 2)};
    Vec3 b1{recip(1, 0), recip(1, 1), recip(1, 2)};
    Vec3 b2{recip(2, 0), recip(2, 1), recip(2, 2)};
    int gmax =
        static_cast<int>(g_cut / std::min({b0.norm(), b1.norm(), b2.norm()}))
        + 1;

    for (int hx = -gmax; hx <= gmax; ++hx)
        for (int hy = -gmax; hy <= gmax; ++hy)
            for (int hz = -gmax; hz <= gmax; ++hz) {
                if (hx == 0 && hy == 0 && hz == 0)
                    continue;
                Vec3 g = b0 * AccReal(hx) + b1 * AccReal(hy) + b2 * AccReal(hz);
                AccReal g2 = g.norm_squared();
                if (g2 > g_cut * g_cut)
                    continue;
                std::complex<AccReal> sg(0, 0);
                for (int i = 0; i < n; ++i) {
                    AccReal phase =
                        g[0] * pos[i][0] + g[1] * pos[i][1] + g[2] * pos[i][2];
                    sg += std::complex<AccReal>(charges[i] * std::cos(phase),
                                                -charges[i] * std::sin(phase));
                }
                AccReal sgn = std::norm(sg);
                e_recip_ += sgn * std::exp(-g2 / alpha_sq_4) / g2;
            }
    e_recip_ *= AccReal(2) * AccReal(units::pi) / vol;

    // Self energy
    e_self_ = 0;
    for (int i = 0; i < n; ++i)
        e_self_ -=
            charges[i] * charges[i] * alpha / std::sqrt(AccReal(units::pi));
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
