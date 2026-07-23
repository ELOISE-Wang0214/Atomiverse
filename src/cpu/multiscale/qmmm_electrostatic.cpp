/// @file src/cpu/multiscale/qmmm_electrostatic.cpp
#include "atomiverse/cpu/multiscale/qmmm_electrostatic.h"

#include <cmath>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {

QmmmElectrostatic::QmmmElectrostatic(std::vector<TagInt> qm_tags,
                                     double mm_charge, double epsilon,
                                     double sigma, double cutoff)
    : mm_q_(mm_charge), eps_(epsilon), sig_(sigma), rcut_(cutoff) {
    if (qm_tags.empty())
        throw std::invalid_argument("QmmmElectrostatic: empty QM tags");
    if (!(eps_ > 0.0) || !(sig_ > 0.0) || !(cutoff > 0.0))
        throw std::invalid_argument("QmmmElectrostatic: invalid LJ");
    for (TagInt t : qm_tags)
        qm_.insert(t);
}

void QmmmElectrostatic::build(const BuildContext&) {}

void QmmmElectrostatic::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = f.nlocal();
    const auto tags = f.span(kAtom_tag);
    const auto type = f.span(kAtom_type);
    const StorageReal* x = f.span(kAtom_x);
    const StorageReal* y = f.span(kAtom_y);
    const StorageReal* z = f.span(kAtom_z);
    double& E = f.meta.get(kEnergy);

    std::vector<int> qm_rows, mm_rows;
    for (int i = 0; i < n; ++i)
        (qm_.count(tags[i]) ? qm_rows : mm_rows).push_back(i);
    if (qm_rows.empty())
        throw std::runtime_error("QmmmElectrostatic: no QM atoms");

    cpu::Frame fqm(static_cast<int>(qm_rows.size()));
    for (std::size_t k = 0; k < qm_rows.size(); ++k) {
        const int i = qm_rows[k];
        fqm.span(kAtom_x)[k] = x[i];
        fqm.span(kAtom_y)[k] = y[i];
        fqm.span(kAtom_z)[k] = z[i];
        fqm.span(kAtom_type)[k] = type[i];
        fqm.span(kAtom_tag)[k] = static_cast<TagInt>(k + 1);
        fqm.span(kAtom_mass)[k] = 1.0;
    }
    fqm.simbox = f.simbox;

    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-8;
    rhf.max_iter = 80;
    rhf.point_charges.clear();
    for (int i : mm_rows) {
        cpu::scf::RHF::PointCharge pc;
        pc.x = x[i];
        pc.y = y[i];
        pc.z = z[i];
        pc.q = mm_q_;
        rhf.point_charges.push_back(pc);
    }
    rhf.run(fqm);
    e_qm_ = fqm.meta.get(kEnergy); // eV from method writeback
    E += e_qm_;

    // LJ MM + cross (not QM–QM)
    e_mm_ = 0;
    const double rc2 = rcut_ * rcut_;
    const double sig2 = sig_ * sig_;
    const double s2c = sig2 / rc2;
    const double s6c = s2c * s2c * s2c;
    const double vshift = 4.0 * eps_ * (s6c * s6c - s6c);
    for (int i = 0; i < n; ++i) {
        const bool iqm = qm_.count(tags[i]) > 0;
        for (int j = i + 1; j < n; ++j) {
            const bool jqm = qm_.count(tags[j]) > 0;
            if (iqm && jqm)
                continue;
            double dx = x[j] - x[i], dy = y[j] - y[i], dz = z[j] - z[i];
            f.simbox.min_image(dx, dy, dz);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= rc2 || r2 == 0.0)
                continue;
            const double r2inv = 1.0 / r2;
            const double s2 = sig2 * r2inv;
            const double s6 = s2 * s2 * s2;
            e_mm_ += 4.0 * eps_ * (s6 * s6 - s6) - vshift;
        }
    }
    E += e_mm_;
}

} // namespace cpu
} // namespace ATV
