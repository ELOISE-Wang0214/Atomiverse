/// @file cpu/gpw/internal/gpw_stress.cpp
/// @brief GPW cell stress via central FD of SCF energy under diagonal strain.
///
/// v1 delivers a correct orthorhombic stress oracle (same strategy as PW
/// strain-FD). Strain is applied to each box edge at fixed fractional
/// coordinates. Analytic Pulay / Nielsen–Martin decomposition is a follow-on;
/// shear (off-diagonal) needs a triclinic GPW cell and is left 0 here.

#include "atomiverse/cpu/gpw/internal/gpw_stress.h"

#include <cmath>
#include <stdexcept>

namespace ATV {
namespace cpu {
namespace gpw {
namespace internal {
namespace {

GpwDft copy_gpw(const GpwDft& m) {
    GpwDft p;
    for (int d = 0; d < 3; ++d)
        p.cell[d] = m.cell[d];
    p.cutoff_ry = m.cutoff_ry;
    p.ngrids = m.ngrids;
    p.rel_cutoff = m.rel_cutoff;
    p.n_occ = m.n_occ;
    p.n_alpha = m.n_alpha;
    p.n_beta = m.n_beta;
    p.max_iter = m.max_iter;
    p.tol = m.tol;
    p.eps_default = m.eps_default;
    p.mix = m.mix;
    p.diis_depth = m.diis_depth;
    p.lindep_threshold = m.lindep_threshold;
    p.level_shift = m.level_shift;
    p.xc = m.xc;
    p.nonperiodic = m.nonperiodic;
    p.local_pp_on_grid = m.local_pp_on_grid;
    p.d3 = m.d3;
    p.d3bj = m.d3bj;
    p.use_ot = m.use_ot;
    return p;
}

/// Scale atom positions and all basis-function centers along axis `ax` by
/// factor `s` (fixed fractional coordinates under orthorhombic strain).
std::vector<GpwAtom> scale_axis(const std::vector<GpwAtom>& atoms, int ax,
                                double s) {
    std::vector<GpwAtom> out = atoms;
    for (auto& a : out) {
        a.pos[ax] *= s;
        for (auto& bf : a.basis)
            for (auto& p : bf)
                p.center[ax] *= s;
    }
    return out;
}

double scf_energy_ha(const std::vector<GpwAtom>& atoms, const GpwDft& method) {
    GpwDft solver = copy_gpw(method);
    solver.build(atoms);
    const GpwResult& r = solver.run(atoms);
    if (!r.converged)
        throw std::runtime_error("gpw_stress: SCF did not converge");
    return r.e_total;
}

} // namespace

std::array<double, 6> gpw_stress_impl(const std::vector<GpwAtom>& atoms,
                                      const GpwDft& method, double h) {
    if (atoms.empty())
        throw std::invalid_argument("gpw_stress: empty atom list");
    if (!(h > 0.0) || !(h < 0.1))
        throw std::invalid_argument("gpw_stress: strain step h out of range");

    const double vol = method.cell[0] * method.cell[1] * method.cell[2]; // Bohr³
    if (!(vol > 0.0))
        throw std::invalid_argument("gpw_stress: non-positive cell volume");

    std::array<double, 6> stress{};
    // Diagonal Voigt: xx,yy,zz
    for (int a = 0; a < 3; ++a) {
        const double sp = 1.0 + h;
        const double sm = 1.0 - h;

        std::vector<GpwAtom> plus = scale_axis(atoms, a, sp);
        std::vector<GpwAtom> minus = scale_axis(atoms, a, sm);

        GpwDft mp = copy_gpw(method);
        GpwDft mm = copy_gpw(method);
        mp.cell[a] = method.cell[a] * sp;
        mm.cell[a] = method.cell[a] * sm;

        const double ep = scf_energy_ha(plus, mp);
        const double em = scf_energy_ha(minus, mm);
        // CP2K STRESS_TENSOR convention (matches NUMERICAL twin): 
        // σ_aa = −(1/Ω) ∂E/∂ε_aa with E in Hartree, Ω in Bohr³.
        // (Positive σ when expanding lowers energy — same sign as CP2K bar print.)
        stress[static_cast<std::size_t>(a)] = -(ep - em) / (2.0 * h) / vol;
    }
    // Off-diagonal shear: orthorhombic GPW has no triclinic cell — leave 0.
    stress[3] = stress[4] = stress[5] = 0.0;
    return stress;
}

} // namespace internal

std::array<double, 6> gpw_stress(const std::vector<GpwAtom>& atoms,
                                 const GpwDft& method, double h) {
    return internal::gpw_stress_impl(atoms, method, h);
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
