/// @file cpu/recorder/cube.h
/// @brief Gaussian ".cube" volumetric file I/O helpers.
///
/// Shared by `DensityCubeWriter` and `OrbitalCubeWriter`. The cube format is
/// fixed: all coordinates are written in Bohr (the SCF convention), values
/// are printed 6 per line in exponential form.
///
/// Specification: http://paulbourke.net/dataformats/cube/
#pragma once

#include <array>
#include <string>
#include <vector>

namespace ATV {
namespace cpu {
namespace scf {
struct IntegralPack;
}
namespace recorder {

/// @brief Uniform Cartesian grid in Bohr (cube-file native units).
struct CubeGrid {
    std::array<double, 3> origin;  ///< Grid origin [Bohr].
    std::array<int, 3> dims;       ///< Number of samples along x,y,z.
    std::array<double, 3> spacing; ///< Sample spacing [Bohr] along each axis.

    /// @brief Grid voxel volume [Bohr^3].
    double voxel() const { return spacing[0] * spacing[1] * spacing[2]; }

    /// @brief Total number of samples.
    long n_points() const {
        return static_cast<long>(dims[0]) * dims[1] * dims[2];
    }
};

/// @brief Build a uniform grid that wraps the molecule with a padding margin.
///
/// @param integrals Integral pack (used for atom positions -- already in Bohr).
/// @param padding   Extra margin around the molecule bounding box [Bohr].
/// @param spacing   Sample spacing [Bohr] (isotropic for all three axes).
/// @return CubeGrid centered on the bounding box.
CubeGrid bounding_box_grid(const scf::IntegralPack& integrals,
                           double padding = 4.0, double spacing = 0.25);

/// @brief Write a Gaussian cube file.
///
/// @param filename Output path.
/// @param title1   First title line.
/// @param title2   Second title line.
/// @param integrals Integral pack -- source of atom positions (Bohr) and Z.
/// @param grid     Grid description (all Bohr).
/// @param values   Grid values in row-major order (ix, iy, iz) with
///                 iz the fastest axis -- cube-file convention.
/// @throws std::runtime_error if the file cannot be opened.
void write_cube(const std::string& filename, const std::string& title1,
                const std::string& title2, const scf::IntegralPack& integrals,
                const CubeGrid& grid, const std::vector<double>& values);

} // namespace recorder
} // namespace cpu
} // namespace ATV
