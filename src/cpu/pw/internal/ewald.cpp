/// @file cpu/pw/internal/ewald.cpp
/// @brief 3D Ewald summation implementation.

#include "atomiverse/cpu/pw/internal/ewald.h"
#include "atomiverse/core/constants.h"
#include <cmath>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

EwaldSummation::EwaldSummation(const std::vector<Vec3>& pos_ang,
                               const std::vector<Real>& charges,
                               const ::ATV::Box& cell, Real ecut) {
    // Convert to Bohr internally
    Real a2b = Real(1.8897259886);
    Box cell_bohr = cell;
    cell_bohr.cell = cell_bohr.cell * a2b;

    std::vector<Vec3> pos(charges.size());
    for (size_t i = 0; i < charges.size(); ++i)
        pos[i] = {pos_ang[i][0]*a2b, pos_ang[i][1]*a2b, pos_ang[i][2]*a2b};

    Real vol = cell_bohr.volume();
    if (vol <= 0) return;

    // Ewald alpha: auto-select, calibrated to match QE ewlad.f90
    Real rcut = cell_bohr.lengths().norm() * Real(0.5);
    alpha_ = std::sqrt(-std::log(Real(2e-17))) / rcut;
    Real alpha = alpha_;
    Real alpha_sq_4 = Real(4) * alpha * alpha;

    int n = static_cast<int>(charges.size());

    // Real space
    e_real_ = 0;
    Real r_cut = Real(10); // Angstrom in Bohr units? Let's use cell-based cutoff
    // actually: r_cut such that erfc(alpha*r) < 1e-12, roughly r_cut ~ 5/alpha
    r_cut = Real(5) / alpha;
    auto lens = cell_bohr.lengths();
    int ncell_r = static_cast<int>(std::ceil(r_cut / std::min({lens[0], lens[1], lens[2]}))) + 1;

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            Real zi = charges[i], zj = charges[j];
            Real fac = (i == j) ? Real(0.5) : Real(1);
            for (int hx = -ncell_r; hx <= ncell_r; ++hx)
                for (int hy = -ncell_r; hy <= ncell_r; ++hy)
                    for (int hz = -ncell_r; hz <= ncell_r; ++hz) {
                        Vec3 dr = {
                            pos[i][0] - pos[j][0] + Real(hx)*cell_bohr.cell(0,0) + Real(hy)*cell_bohr.cell(0,1) + Real(hz)*cell_bohr.cell(0,2),
                            pos[i][1] - pos[j][1] + Real(hx)*cell_bohr.cell(1,0) + Real(hy)*cell_bohr.cell(1,1) + Real(hz)*cell_bohr.cell(1,2),
                            pos[i][2] - pos[j][2] + Real(hx)*cell_bohr.cell(2,0) + Real(hy)*cell_bohr.cell(2,1) + Real(hz)*cell_bohr.cell(2,2)
                        };
                        Real r = dr.norm();
                        if (r > 0 && r < r_cut) {
                            e_real_ += fac * zi * zj * std::erfc(alpha * r) / r;
                        }
                        if (r < Real(1e-12) && i == j && hx==0 && hy==0 && hz==0) {
                            // Self term handled separately
                        }
                    }
        }
    }

    // Reciprocal space
    e_recip_ = 0;
    Real g_cut = std::sqrt(Real(8) * ecut); // 4x E_cut for Ewald recip
    auto recip = cell_bohr.recip_cell();
    Vec3 b0{recip(0,0), recip(0,1), recip(0,2)};
    Vec3 b1{recip(1,0), recip(1,1), recip(1,2)};
    Vec3 b2{recip(2,0), recip(2,1), recip(2,2)};
    int gmax = static_cast<int>(g_cut / std::min({b0.norm(), b1.norm(), b2.norm()})) + 1;

    for (int hx = -gmax; hx <= gmax; ++hx)
        for (int hy = -gmax; hy <= gmax; ++hy)
            for (int hz = -gmax; hz <= gmax; ++hz) {
                if (hx == 0 && hy == 0 && hz == 0) continue;
                Vec3 g = b0*Real(hx) + b1*Real(hy) + b2*Real(hz);
                Real g2 = g.norm_squared();
                if (g2 > g_cut * g_cut) continue;
                std::complex<Real> sg(0, 0);
                for (int i = 0; i < n; ++i) {
                    Real phase = g[0]*pos[i][0] + g[1]*pos[i][1] + g[2]*pos[i][2];
                    sg += std::complex<Real>(charges[i] * std::cos(phase),
                                             -charges[i] * std::sin(phase));
                }
                Real sgn = std::norm(sg);
                e_recip_ += sgn * std::exp(-g2 / alpha_sq_4) / g2;
            }
    e_recip_ *= Real(2) * Real(kPI) / vol;

    // Self energy
    e_self_ = 0;
    for (int i = 0; i < n; ++i)
        e_self_ -= charges[i] * charges[i] * alpha / std::sqrt(Real(kPI));
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
