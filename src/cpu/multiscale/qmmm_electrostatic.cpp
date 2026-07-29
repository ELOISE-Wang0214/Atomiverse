/// @file src/cpu/multiscale/qmmm_electrostatic.cpp
#include "atomiverse/cpu/multiscale/qmmm_electrostatic.h"

#include <cmath>
#include <stdexcept>
#include <unordered_map>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {

QmmmElectrostatic::QmmmElectrostatic(std::vector<TagInt> qm_tags,
                                     double mm_charge, double epsilon,
                                     double sigma, double cutoff,
                                     std::vector<LinkBond> links,
                                     double link_scale)
    : links_(std::move(links)), mm_q_(mm_charge), eps_(epsilon), sig_(sigma),
      rcut_(cutoff), link_scale_(link_scale) {
    if (qm_tags.empty())
        throw std::invalid_argument("QmmmElectrostatic: empty QM tags");
    if (!(eps_ > 0.0) || !(sig_ > 0.0) || !(cutoff > 0.0))
        throw std::invalid_argument("QmmmElectrostatic: invalid LJ");
    if (!(link_scale_ > 0.0) || !(link_scale_ < 1.0))
        throw std::invalid_argument(
            "QmmmElectrostatic: link_scale must be in (0,1)");
    for (TagInt t : qm_tags)
        qm_.insert(t);
    for (const auto& lb : links_) {
        if (lb.qm_tag == lb.mm_tag)
            throw std::invalid_argument(
                "QmmmElectrostatic: link bond with identical tags");
        if (!qm_.count(lb.qm_tag))
            throw std::invalid_argument(
                "QmmmElectrostatic: link qm_tag not in QM set");
        if (qm_.count(lb.mm_tag))
            throw std::invalid_argument(
                "QmmmElectrostatic: link mm_tag must not be in QM set");
    }
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

    std::unordered_map<TagInt, int> tag_to_row;
    std::vector<int> qm_rows, mm_rows;
    for (int i = 0; i < n; ++i) {
        tag_to_row[tags[i]] = i;
        (qm_.count(tags[i]) ? qm_rows : mm_rows).push_back(i);
    }
    if (qm_rows.empty())
        throw std::runtime_error("QmmmElectrostatic: no QM atoms");

    // Resolve link bonds; MM partners are excluded from point-charge embedding.
    std::unordered_set<TagInt> link_mm;
    struct ResolvedLink {
        int qm_row;
        int mm_row;
        int fqm_h_row; // index of link H in the private QM frame
    };
    std::vector<ResolvedLink> resolved;
    resolved.reserve(links_.size());
    for (const auto& lb : links_) {
        auto iq = tag_to_row.find(lb.qm_tag);
        auto im = tag_to_row.find(lb.mm_tag);
        if (iq == tag_to_row.end() || im == tag_to_row.end())
            throw std::runtime_error(
                "QmmmElectrostatic: link bond atom missing from frame");
        link_mm.insert(lb.mm_tag);
        resolved.push_back({iq->second, im->second, -1});
    }
    last_n_links_ = static_cast<int>(resolved.size());

    const int n_qm_real = static_cast<int>(qm_rows.size());
    const int n_fqm = n_qm_real + last_n_links_;
    cpu::Frame fqm(n_fqm);
    auto qx = fqm.span(kAtom_x);
    auto qy = fqm.span(kAtom_y);
    auto qz = fqm.span(kAtom_z);
    auto qtype = fqm.span(kAtom_type);
    auto qtag = fqm.span(kAtom_tag);
    auto qmass = fqm.span(kAtom_mass);
    auto qfx = fqm.span(kAtom_fx);
    auto qfy = fqm.span(kAtom_fy);
    auto qfz = fqm.span(kAtom_fz);
    for (std::size_t k = 0; k < qm_rows.size(); ++k) {
        const int i = qm_rows[k];
        qx[k] = x[i];
        qy[k] = y[i];
        qz[k] = z[i];
        qtype[k] = type[i];
        qtag[k] = static_cast<TagInt>(k + 1);
        qmass[k] = 1.0;
        qfx[k] = qfy[k] = qfz[k] = 0;
    }
    for (std::size_t li = 0; li < resolved.size(); ++li) {
        const int hrow = n_qm_real + static_cast<int>(li);
        resolved[li].fqm_h_row = hrow;
        const int iq = resolved[li].qm_row;
        const int im = resolved[li].mm_row;
        double hx, hy, hz;
        place_link_hydrogen(x[iq], y[iq], z[iq], x[im], y[im], z[im],
                            link_scale_, hx, hy, hz);
        qx[hrow] = hx;
        qy[hrow] = hy;
        qz[hrow] = hz;
        qtype[hrow] = 1; // H
        qtag[hrow] = static_cast<TagInt>(hrow + 1);
        qmass[hrow] = 1.0;
        qfx[hrow] = qfy[hrow] = qfz[hrow] = 0;
    }
    fqm.box = f.box;

    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-8;
    rhf.max_iter = 80;
    rhf.compute_force = ctx.forces_requested;
    rhf.point_charges.clear();
    for (int i : mm_rows) {
        if (link_mm.count(tags[i]))
            continue; // link MM partner excluded from embedding charges
        cpu::EmbeddingSite pc;
        pc.x = x[i];
        pc.y = y[i];
        pc.z = z[i];
        pc.q = mm_q_;
        rhf.point_charges.push_back(pc);
    }
    rhf.run(fqm);
    e_qm_ = fqm.meta.get(kEnergy); // eV from method writeback
    E += e_qm_;

    AccReal* fx = nullptr;
    AccReal* fy = nullptr;
    AccReal* fz = nullptr;
    if (ctx.forces_requested) {
        fx = f.span(kAtom_fx);
        fy = f.span(kAtom_fy);
        fz = f.span(kAtom_fz);
    }

    // Project QM forces (including link H) back onto the real Frame.
    if (fx) {
        const AccReal* qfx = fqm.span(kAtom_fx);
        const AccReal* qfy = fqm.span(kAtom_fy);
        const AccReal* qfz = fqm.span(kAtom_fz);
        for (std::size_t k = 0; k < qm_rows.size(); ++k) {
            const int i = qm_rows[k];
            fx[i] += qfx[k];
            fy[i] += qfy[k];
            fz[i] += qfz[k];
        }
        for (const auto& rl : resolved) {
            const int h = rl.fqm_h_row;
            project_link_force(link_scale_, qfx[h], qfy[h], qfz[h],
                               fx[rl.qm_row], fy[rl.qm_row], fz[rl.qm_row],
                               fx[rl.mm_row], fy[rl.mm_row], fz[rl.mm_row]);
        }
    }

    // LJ MM + cross (not QM–QM). Force pointers resolved once above.
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
            f.box.min_image(dx, dy, dz);
            const double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= rc2 || r2 == 0.0)
                continue;
            const double r2inv = 1.0 / r2;
            const double s2 = sig2 * r2inv;
            const double s6 = s2 * s2 * s2;
            e_mm_ += 4.0 * eps_ * (s6 * s6 - s6) - vshift;
            if (fx) {
                const double fpair = 24.0 * eps_ * r2inv * (2.0 * s6 * s6 - s6);
                fx[i] -= fpair * dx;
                fy[i] -= fpair * dy;
                fz[i] -= fpair * dz;
                fx[j] += fpair * dx;
                fy[j] += fpair * dy;
                fz[j] += fpair * dz;
            }
        }
    }
    E += e_mm_;
}

} // namespace cpu
} // namespace ATV
