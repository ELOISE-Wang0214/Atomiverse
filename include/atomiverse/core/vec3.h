/// @file core/vec3.h
/// @brief Vec3 — 3D vector value type.
#pragma once

#include <xtensor/xfixed.hpp>
#include <cmath>
#include "atomiverse/core/real.h"

namespace ATV {

class Vec3 {
    xt::xtensor_fixed<Real, xt::xshape<3>> v_{};
public:
    Vec3() = default;
    Vec3(Real x, Real y, Real z) { v_ = {x, y, z}; }

    Real& operator[](int i) { return v_(i); }
    const Real& operator[](int i) const { return v_(i); }

    Real& x() { return v_(0); }
    Real& y() { return v_(1); }
    Real& z() { return v_(2); }
    Real x() const { return v_(0); }
    Real y() const { return v_(1); }
    Real z() const { return v_(2); }

    Vec3 operator+(const Vec3& o) const { return {v_(0)+o.v_(0), v_(1)+o.v_(1), v_(2)+o.v_(2)}; }
    Vec3 operator-(const Vec3& o) const { return {v_(0)-o.v_(0), v_(1)-o.v_(1), v_(2)-o.v_(2)}; }
    Vec3 operator*(Real s) const { return {v_(0)*s, v_(1)*s, v_(2)*s}; }
    Vec3 operator/(Real s) const { return {v_(0)/s, v_(1)/s, v_(2)/s}; }

    Real dot(const Vec3& o) const {
        return v_(0)*o.v_(0) + v_(1)*o.v_(1) + v_(2)*o.v_(2);
    }
    Vec3 cross(const Vec3& o) const {
        return {
            v_(1)*o.v_(2) - v_(2)*o.v_(1),
            v_(2)*o.v_(0) - v_(0)*o.v_(2),
            v_(0)*o.v_(1) - v_(1)*o.v_(0)
        };
    }
    Real norm_squared() const { return dot(*this); }
    Real norm() const { return std::sqrt(norm_squared()); }
    Vec3 normalized() const { Real n = norm(); return {v_(0)/n, v_(1)/n, v_(2)/n}; }

    static Vec3 zeros() { return {0, 0, 0}; }
    static Vec3 unit_x() { return {1, 0, 0}; }
    static Vec3 unit_y() { return {0, 1, 0}; }
    static Vec3 unit_z() { return {0, 0, 1}; }
};

} // namespace ATV
