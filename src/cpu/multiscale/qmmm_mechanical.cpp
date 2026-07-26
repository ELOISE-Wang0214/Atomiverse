/// @file src/cpu/multiscale/qmmm_mechanical.cpp
#include "atomiverse/cpu/multiscale/qmmm_mechanical.h"

#include <cmath>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {

QmmmMechanical::QmmmMechanical(std::vector<TagInt> qm_tags, double epsilon,
                               double sigma, double cutoff)
    : eps_(epsilon), sig_(sigma), rcut_(cutoff) {
    if (qm_tags.empty())
        throw std::invalid_argument("QmmmMechanical: empty QM tag set");
    if (!(eps_ > 0.0) || !(sig_ > 0.0) || !(cutoff > 0.0))
        throw std::invalid_argument("QmmmMechanical: invalid LJ parameters");
    for (TagInt t : qm_tags)
        qm_.insert(t);
    dftb_.with_forces = false;
}

void QmmmMechanical::build(const BuildContext& ctx) { dftb_.build(ctx); }

void QmmmMechanical::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = f.nlocal();
    const auto tags = f.span(kAtom_tag);
    const auto type = f.span(kAtom_type);
    const StorageReal* x = f.span(kAtom_x);
    const StorageReal* y = f.span(kAtom_y);
    const StorageReal* z = f.span(kAtom_z);
    double& E = f.meta.get(kEnergy);

    // Build a temporary QM-only frame for DFTB2.
    std::vector<int> qm_rows;
    for (int i = 0; i < n; ++i)
        if (qm_.count(tags[i]))
            qm_rows.push_back(i);
    if (qm_rows.empty())
        throw std::runtime_error("QmmmMechanical: no QM atoms local");

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
    fqm.box = f.box;
    fqm.meta.get(kEnergy) = 0;
    {
        LaunchContext lqm{ctx.runtime, fqm, ctx.stage, ctx.step};
        lqm.forces_requested = false;
        dftb_.launch(lqm);
    }
    e_qm_ = fqm.meta.get(kEnergy);
    E += e_qm_;

    // LJ mechanical: all pairs except pure QM–QM (those are in QM energy).
    e_mm_ = 0;
    const double rc2 = rcut_ * rcut_;
    const double sig2 = sig_ * sig_;
    const double rc = rcut_;
    const double s2c = sig2 / rc2;
    const double s6c = s2c * s2c * s2c;
    const double vshift = 4.0 * eps_ * (s6c * s6c - s6c);
    for (int i = 0; i < n; ++i) {
        const bool iqm = qm_.count(tags[i]) > 0;
        for (int j = i + 1; j < n; ++j) {
            const bool jqm = qm_.count(tags[j]) > 0;
            if (iqm && jqm)
                continue; // mechanical: no classical QM–QM
            double dx = x[j] - x[i], dy = y[j] - y[i], dz = z[j] - z[i];
            f.box.min_image(dx, dy, dz);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= rc2 || r2 == 0.0)
                continue;
            const double r2inv = 1.0 / r2;
            const double s2 = sig2 * r2inv;
            const double s6 = s2 * s2 * s2;
            const double s12 = s6 * s6;
            e_mm_ += 4.0 * eps_ * (s12 - s6) - vshift;
            if (ctx.forces_requested) {
                const double fpair = 24.0 * eps_ * r2inv * (2.0 * s12 - s6);
                f.span(kAtom_fx)[i] -= fpair * dx;
                f.span(kAtom_fy)[i] -= fpair * dy;
                f.span(kAtom_fz)[i] -= fpair * dz;
                f.span(kAtom_fx)[j] += fpair * dx;
                f.span(kAtom_fy)[j] += fpair * dy;
                f.span(kAtom_fz)[j] += fpair * dz;
            }
        }
    }
    E += e_mm_;
}

} // namespace cpu
} // namespace ATV
