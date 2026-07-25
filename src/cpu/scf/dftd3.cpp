/// @file cpu/scf/dftd3.cpp
/// @brief Molecular D3BJ::energy — Frame (Å) → PBE-D3(BJ) energy [Hartree].

#include "atomiverse/cpu/scf/dftd3.h"

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/gpw/internal/dftd3.h"
#include "atomiverse/cpu/gpw/internal/dftd3_params.h"

namespace ATV {
namespace cpu {
namespace scf {

double D3BJ::energy(const Frame& frame) const {
    const int n = static_cast<int>(frame.nlocal());
    if (n <= 0)
        return 0.0;

    const auto type = frame.span(kAtom_type);
    const auto x = frame.span(kAtom_x);
    const auto y = frame.span(kAtom_y);
    const auto z = frame.span(kAtom_z);

    std::vector<int> Z(static_cast<std::size_t>(n));
    std::vector<std::array<double, 3>> pos(static_cast<std::size_t>(n));
    const double a2b = units::kAngstromInBohr;
    for (int i = 0; i < n; ++i) {
        const int zi = static_cast<int>(type[i]);
        if (!gpw::D3Params::supports_element(zi))
            throw std::runtime_error(
                "D3BJ: element Z=" + std::to_string(zi)
                + " is not in the in-tree D3 table (no fallback)");
        Z[static_cast<std::size_t>(i)] = zi;
        pos[static_cast<std::size_t>(i)] = {x[i] * a2b, y[i] * a2b, z[i] * a2b};
    }

    const double cell[3] = {0.0, 0.0, 0.0};
    const gpw::D3Params params{};
    return gpw::dftd3_pbe_bj(Z, pos, cell, /*periodic=*/false, params,
                             rcut_bohr, s6, s8, a1, a2)
        .energy;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
