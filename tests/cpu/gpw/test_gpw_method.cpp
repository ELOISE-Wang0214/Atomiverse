/// @file tests/cpu/gpw/test_gpw_method.cpp
/// Concrete GpwDft method API: members + build/run, no GpwConfig bag.
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/collocate.h"

using ATV::cpu::gpw::GpwAtom;
using ATV::cpu::gpw::GpwDft;
using ATV::cpu::gpw::GridBasisFunction;
using ATV::cpu::gpw::GridPrimitive;

namespace {

GridBasisFunction s_basis(double alpha, const double pos[3]) {
    GridBasisFunction bf;
    GridPrimitive p;
    p.alpha = alpha;
    p.center[0] = pos[0];
    p.center[1] = pos[1];
    p.center[2] = pos[2];
    p.l[0] = p.l[1] = p.l[2] = 0;
    p.coeff = 1.0;
    bf.push_back(p);
    return bf;
}

GpwAtom hydrogen(const double pos[3]) {
    GpwAtom a;
    a.z_val = 1.0;
    a.z = 1;
    a.pos[0] = pos[0];
    a.pos[1] = pos[1];
    a.pos[2] = pos[2];
    a.basis = {s_basis(0.5, pos)};
    a.rloc = 0.2;
    a.cloc = {-4.17890044, 0.72446331};
    return a;
}

} // namespace

TEST(GpwMethod, MembersAndFiniteEnergy) {
    const double box = 8.0;
    const double rA[3] = {box / 2 - 0.7, box / 2, box / 2};
    const double rB[3] = {box / 2 + 0.7, box / 2, box / 2};
    std::vector<GpwAtom> atoms = {hydrogen(rA), hydrogen(rB)};

    GpwDft m;
    m.cell[0] = m.cell[1] = m.cell[2] = box;
    m.cutoff_ry = 120.0;
    m.n_occ = 1;
    m.max_iter = 50;
    m.tol = 1e-7;
    m.mix = 0.5;
    m.build(atoms);
    const auto& r = m.run(atoms);
    EXPECT_TRUE(std::isfinite(r.e_total));
    EXPECT_TRUE(r.converged);
    EXPECT_EQ(&r, &m.last_result());
}

TEST(GpwMethod, VariableCellRunRebuildsFromCellMember) {
    // run() derives grids / lattice sums from `cell` per call: changing the
    // cell (fixed fractional coordinates) between runs is well-defined and
    // must change the energy — the contract variable-cell relax stands on.
    const double box = 8.0;
    const double rA[3] = {box / 2 - 0.7, box / 2, box / 2};
    const double rB[3] = {box / 2 + 0.7, box / 2, box / 2};
    std::vector<GpwAtom> atoms = {hydrogen(rA), hydrogen(rB)};

    GpwDft m;
    m.cell[0] = m.cell[1] = m.cell[2] = box;
    m.cutoff_ry = 120.0;
    m.n_occ = 1;
    m.max_iter = 50;
    m.tol = 1e-7;
    m.mix = 0.5;
    m.build(atoms);
    const double e0 = m.run(atoms).e_total;

    const double s = 1.05; // strain x only (orthorhombic anisotropic)
    m.cell[0] = box * s;
    std::vector<GpwAtom> strained = atoms;
    for (auto& a : strained) {
        a.pos[0] *= s;
        for (auto& bf : a.basis)
            for (auto& p : bf)
                p.center[0] *= s;
    }
    const auto& r1 = m.run(strained);
    EXPECT_TRUE(r1.converged);
    EXPECT_TRUE(std::isfinite(r1.e_total));
    EXPECT_NE(r1.e_total, e0);
}
