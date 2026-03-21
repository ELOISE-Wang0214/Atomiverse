/// @file cpu/observer/density.cpp
/// @brief Electron density grid evaluation: STO-3G basis function assembly + rho(r).

#include "atomiverse/cpu/observer/density.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "atomiverse/core/constants.h"
#include "atomiverse/core/gto.h"

namespace ATV {
namespace cpu {

namespace {

// ── STO-3G basis data (H, He, Li) ────────────────────────────────────────────

struct STO3GShell {
    int    n_prim;
    double alpha[3];
    double coeff[3];
};

struct STO3GBasis {
    int       n_shells;
    STO3GShell shells[2];
};

// Index by atomic number Z (entry 0 unused).
// Each STO3GBasis has shells[2]; unused second shell is zero-initialized.
static const STO3GBasis kSto3GBasis[] = {
    // Z=0: unused
    {0, {}},
    // Z=1: H — 1 s-shell
    {1, { {3, {3.42525091, 0.62391373, 0.16885540},
              {0.15432897, 0.53532814, 0.44463454}},
          {} }},
    // Z=2: He — 1 s-shell
    {1, { {3, {6.36242139, 1.15892300, 0.31364979},
              {0.15432897, 0.53532814, 0.44463454}},
          {} }},
    // Z=3: Li — 2 s-shells (1s + 2s)
    {2, { {3, {16.1195750,  2.93620070, 0.79465370},
              { 0.15432897,  0.53532814,  0.44463454}},
          {3, { 0.63628970,  0.14786010,  0.04808870},
              {-0.09996723,  0.39951283,  0.70011547}} }},
};

static constexpr int kMaxSupportedZ = 3;

// ── Shell instance at a specific center ──────────────────────────────────────

struct ShellInstance {
    double center_bohr[3];
    double alpha[3];
    double coeff[3];
};

// Evaluate a contracted s-type GTO at (x, y, z) in Bohr.
double eval_shell_s(const ShellInstance& sh,
                    double x_bohr, double y_bohr, double z_bohr)
{
    const double dx = x_bohr - sh.center_bohr[0];
    const double dy = y_bohr - sh.center_bohr[1];
    const double dz = z_bohr - sh.center_bohr[2];
    const double r2 = dx*dx + dy*dy + dz*dz;
    double val = 0.0;
    for (int p = 0; p < 3; ++p)
        val += sh.coeff[p] * phys_norm_s(sh.alpha[p]) * std::exp(-sh.alpha[p] * r2);
    return val;
}

// Build the ordered list of shell instances for the given geometry.
// Returns empty vector if any Z is out of range.
std::vector<ShellInstance> build_shells(
    const std::vector<Int>&  type_id,
    const std::vector<Real>& x,
    const std::vector<Real>& y,
    const std::vector<Real>& z)
{
    std::vector<ShellInstance> shells;
    for (size_t i = 0; i < type_id.size(); ++i) {
        const int Z = static_cast<int>(type_id[i]);
        if (Z <= 0 || Z > kMaxSupportedZ) return {};
        const auto& basis = kSto3GBasis[Z];
        for (int s = 0; s < basis.n_shells; ++s) {
            ShellInstance sh{};
            sh.center_bohr[0] = static_cast<double>(x[i]) * kAngstromToBohr;
            sh.center_bohr[1] = static_cast<double>(y[i]) * kAngstromToBohr;
            sh.center_bohr[2] = static_cast<double>(z[i]) * kAngstromToBohr;
            for (int p = 0; p < 3; ++p) {
                sh.alpha[p] = basis.shells[s].alpha[p];
                sh.coeff[p] = basis.shells[s].coeff[p];
            }
            shells.push_back(sh);
        }
    }
    return shells;
}

} // namespace

// ── eval_density_grid ─────────────────────────────────────────────────────────

DensityGrid eval_density_grid(
    const std::vector<Int>&  type_id,
    const std::vector<Real>& x,
    const std::vector<Real>& y,
    const std::vector<Real>& z,
    const Real* P, Int n_basis,
    int nx, int ny, int nz, Real padding)
{
    DensityGrid grid;
    if (nx <= 0 || ny <= 0 || nz <= 0 || n_basis <= 0 || type_id.empty())
        return grid;

    const auto shells = build_shells(type_id, x, y, z);
    if (shells.empty() || static_cast<Int>(shells.size()) != n_basis)
        return grid;

    // Grid bounds from atom positions + padding
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double min_z = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    double max_z = std::numeric_limits<double>::lowest();
    for (size_t i = 0; i < x.size(); ++i) {
        min_x = std::min(min_x, static_cast<double>(x[i]));
        min_y = std::min(min_y, static_cast<double>(y[i]));
        min_z = std::min(min_z, static_cast<double>(z[i]));
        max_x = std::max(max_x, static_cast<double>(x[i]));
        max_y = std::max(max_y, static_cast<double>(y[i]));
        max_z = std::max(max_z, static_cast<double>(z[i]));
    }
    const double pad = static_cast<double>(padding);
    min_x -= pad; min_y -= pad; min_z -= pad;
    max_x += pad; max_y += pad; max_z += pad;

    grid.nx = nx; grid.ny = ny; grid.nz = nz;
    grid.ox = min_x; grid.oy = min_y; grid.oz = min_z;
    grid.lx = std::max(max_x - min_x, 1e-6);
    grid.ly = std::max(max_y - min_y, 1e-6);
    grid.lz = std::max(max_z - min_z, 1e-6);

    const size_t n  = static_cast<size_t>(n_basis);
    const size_t total = static_cast<size_t>(nx) *
                         static_cast<size_t>(ny) *
                         static_cast<size_t>(nz);
    grid.rho.assign(total, 0.0);

    for (int ix = 0; ix < nx; ++ix) {
        const double rx = min_x + (static_cast<double>(ix) / nx) * grid.lx;
        const double rx_bohr = rx * kAngstromToBohr;
        for (int iy = 0; iy < ny; ++iy) {
            const double ry = min_y + (static_cast<double>(iy) / ny) * grid.ly;
            const double ry_bohr = ry * kAngstromToBohr;
            for (int iz = 0; iz < nz; ++iz) {
                const double rz = min_z + (static_cast<double>(iz) / nz) * grid.lz;
                const double rz_bohr = rz * kAngstromToBohr;

                // Evaluate all basis functions at this point
                double phi[32] = {};  // max n_basis for H/He/Li molecules
                for (size_t mu = 0; mu < n; ++mu)
                    phi[mu] = eval_shell_s(shells[mu], rx_bohr, ry_bohr, rz_bohr);

                // rho(r) = sum_{mu,nu} P_{mu,nu} phi_mu phi_nu
                double val = 0.0;
                for (size_t mu = 0; mu < n; ++mu)
                    for (size_t nu = 0; nu < n; ++nu)
                        val += static_cast<double>(P[mu * n + nu]) * phi[mu] * phi[nu];

                grid.rho[(static_cast<size_t>(ix) * ny + iy) * nz + iz] = val;
            }
        }
    }
    return grid;
}

} // namespace cpu
} // namespace ATV
