/// @file cpu/recorder/cube.cpp
/// @brief Gaussian cube-format I/O.

#include "atomiverse/cpu/recorder/cube.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/native_integrals.h"

namespace ATV {
namespace cpu {
namespace recorder {

CubeGrid bounding_box_grid(const scf::IntegralPack& integrals, double padding,
                           double spacing) {
    if (integrals.atoms.empty())
        throw std::runtime_error("bounding_box_grid: empty integral pack");
    if (spacing <= 0.0)
        throw std::invalid_argument("bounding_box_grid: spacing must be > 0");

    double xmin = integrals.atoms[0].pos[0];
    double xmax = xmin;
    double ymin = integrals.atoms[0].pos[1];
    double ymax = ymin;
    double zmin = integrals.atoms[0].pos[2];
    double zmax = zmin;
    for (const auto& a : integrals.atoms) {
        xmin = std::min(xmin, a.pos[0]);
        xmax = std::max(xmax, a.pos[0]);
        ymin = std::min(ymin, a.pos[1]);
        ymax = std::max(ymax, a.pos[1]);
        zmin = std::min(zmin, a.pos[2]);
        zmax = std::max(zmax, a.pos[2]);
    }
    CubeGrid g{};
    g.origin = {xmin - padding, ymin - padding, zmin - padding};
    double lx = (xmax - xmin) + 2 * padding;
    double ly = (ymax - ymin) + 2 * padding;
    double lz = (zmax - zmin) + 2 * padding;
    g.dims[0] = std::max(2, static_cast<int>(std::ceil(lx / spacing)) + 1);
    g.dims[1] = std::max(2, static_cast<int>(std::ceil(ly / spacing)) + 1);
    g.dims[2] = std::max(2, static_cast<int>(std::ceil(lz / spacing)) + 1);
    g.spacing = {spacing, spacing, spacing};
    return g;
}

void write_cube(const std::string& filename, const std::string& title1,
                const std::string& title2, const scf::IntegralPack& integrals,
                const CubeGrid& grid, const std::vector<double>& values) {
    if (static_cast<long>(values.size()) != grid.n_points())
        throw std::invalid_argument(
            "write_cube: values.size() != grid.n_points()");

    std::FILE* fp = std::fopen(filename.c_str(), "w");
    if (!fp)
        throw std::runtime_error("write_cube: cannot open " + filename);

    const int natm = static_cast<int>(integrals.atoms.size());
    std::fprintf(fp, "%s\n", title1.c_str());
    std::fprintf(fp, "%s\n", title2.c_str());

    // Header: n_atoms, origin[3].  Negative n_atoms flags "no additional scalar
    // field header"; positive is the standard cube convention used by most
    // visualizers.
    std::fprintf(fp, "%5d %12.6f %12.6f %12.6f\n", natm, grid.origin[0],
                 grid.origin[1], grid.origin[2]);
    // Grid axes: N, dx, dy, dz per line (diagonal axes in our convention).
    std::fprintf(fp, "%5d %12.6f %12.6f %12.6f\n", grid.dims[0],
                 grid.spacing[0], 0.0, 0.0);
    std::fprintf(fp, "%5d %12.6f %12.6f %12.6f\n", grid.dims[1], 0.0,
                 grid.spacing[1], 0.0);
    std::fprintf(fp, "%5d %12.6f %12.6f %12.6f\n", grid.dims[2], 0.0, 0.0,
                 grid.spacing[2]);

    // Atom block: Z, charge (placeholder = Z), x, y, z (Bohr).
    for (const auto& a : integrals.atoms) {
        std::fprintf(fp, "%5d %12.6f %12.6f %12.6f %12.6f\n", a.Z,
                     static_cast<double>(a.Z), a.pos[0], a.pos[1], a.pos[2]);
    }

    // Values: iz is the fastest-running axis. Format: 6 values per line,
    // scientific notation.
    const int nx = grid.dims[0];
    const int ny = grid.dims[1];
    const int nz = grid.dims[2];
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            int col = 0;
            for (int iz = 0; iz < nz; ++iz) {
                long idx = (static_cast<long>(ix) * ny + iy) * nz + iz;
                std::fprintf(fp, " %12.5e", values[idx]);
                ++col;
                if (col == 6) {
                    std::fputc('\n', fp);
                    col = 0;
                }
            }
            if (col != 0)
                std::fputc('\n', fp);
        }
    }
    std::fclose(fp);
}

} // namespace recorder
} // namespace cpu
} // namespace ATV
