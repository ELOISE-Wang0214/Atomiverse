/// @file core/mat3.h
/// @brief Mat3 -- 3x3 matrix for simulation box representation.
///
/// Stores cell vectors in row-major order. Used for both orthorhombic
/// (diagonal) and triclinic simulation boxes.
#pragma once

#include "atomiverse/core/real.h"

namespace ATV {

/// @brief 3x3 matrix -- simulation box cell vectors [Angstrom].
///
/// Row-major storage: row i is the i-th cell vector.
/// For orthorhombic boxes, use Mat3::diagonal() or Mat3::cubic().
struct Mat3 {
    Real data[3][3] = {}; ///< Row-major 3x3 storage [Angstrom].

    /// @brief Element access.
    /// @param i Row index in `[0, 2]`.
    /// @param j Column index in `[0, 2]`.
    /// @return Mutable reference to `data[i][j]`.
    Real& operator()(int i, int j) { return data[i][j]; }

    /// @brief Const element access.
    /// @param i Row index in `[0, 2]`.
    /// @param j Column index in `[0, 2]`.
    /// @return Read-only reference to `data[i][j]`.
    const Real& operator()(int i, int j) const { return data[i][j]; }

    /// @brief Construct a diagonal matrix (orthorhombic box).
    /// @param a Box length along x [Angstrom].
    /// @param b Box length along y [Angstrom].
    /// @param c Box length along z [Angstrom].
    /// @return Diagonal `Mat3` with those axis lengths.
    static Mat3 diagonal(Real a, Real b, Real c) {
        Mat3 m{};
        m.data[0][0] = a; m.data[1][1] = b; m.data[2][2] = c;
        return m;
    }

    /// @brief Construct a cubic box.
    /// @param l Edge length [Angstrom].
    /// @return Cubic `Mat3` with edge length `l`.
    static Mat3 cubic(Real l) { return diagonal(l, l, l); }
};

} // namespace ATV
