/// @file cpu/observer/mulliken.cpp
/// @brief Mulliken population analysis — delegates to molrs via CXX bridge.

#include "atomiverse/cpu/observer/mulliken.h"
#include "atomiverse/molrs.hpp"

#include <cstddef>

namespace ATV {
namespace cpu {

// ── Core API (raw pointer, callable from .cu via .cpp helper) ────────────────

std::vector<Real> mulliken_charge(const cpu::Frame& frame,
                                  const Real* P_flat,
                                  const Real* S_flat,
                                  int n_basis)
{
    const int n  = frame.n_atoms();
    const int n2 = n_basis * n_basis;

    // Build type_id slice for molrs (Int -> i32 conversion)
    std::vector<int32_t> type_id_i32(n);
    for (int i = 0; i < n; ++i)
        type_id_i32[i] = static_cast<int32_t>(frame.type_id[i]);

    rust::Slice<const int32_t> tid(type_id_i32.data(), type_id_i32.size());
    rust::Slice<const double>  p(reinterpret_cast<const double*>(P_flat), n2);
    rust::Slice<const double>  s(reinterpret_cast<const double*>(S_flat), n2);

    auto result = molrs::mulliken_charge(tid, p, s, n_basis);

    // Convert rust::Vec<f64> → std::vector<Real>
    std::vector<Real> out(result.size());
    for (std::size_t i = 0; i < result.size(); ++i)
        out[i] = static_cast<Real>(result[i]);
    return out;
}

std::vector<Real> mulliken_spin(const cpu::Frame& frame,
                                const Real* Pa_flat,
                                const Real* Pb_flat,
                                const Real* S_flat,
                                int n_basis)
{
    const int n  = frame.n_atoms();
    const int n2 = n_basis * n_basis;

    std::vector<int32_t> type_id_i32(n);
    for (int i = 0; i < n; ++i)
        type_id_i32[i] = static_cast<int32_t>(frame.type_id[i]);

    rust::Slice<const int32_t> tid(type_id_i32.data(), type_id_i32.size());
    rust::Slice<const double>  pa(reinterpret_cast<const double*>(Pa_flat), n2);
    rust::Slice<const double>  pb(reinterpret_cast<const double*>(Pb_flat), n2);
    rust::Slice<const double>  s(reinterpret_cast<const double*>(S_flat), n2);

    auto result = molrs::mulliken_spin(tid, pa, pb, s, n_basis);

    std::vector<Real> out(result.size());
    for (std::size_t i = 0; i < result.size(); ++i)
        out[i] = static_cast<Real>(result[i]);
    return out;
}

// ── xtensor overloads (for CPU-side callers like examples) ───────────────────

std::vector<Real> mulliken_charge(const cpu::Frame& frame,
                                  const xt::xtensor<Real, 2>& P,
                                  const xt::xtensor<Real, 2>& S)
{
    const int n_basis = static_cast<int>(P.shape(0));
    return mulliken_charge(frame, P.data(), S.data(), n_basis);
}

std::vector<Real> mulliken_spin(const cpu::Frame& frame,
                                const xt::xtensor<Real, 2>& P_alpha,
                                const xt::xtensor<Real, 2>& P_beta,
                                const xt::xtensor<Real, 2>& S)
{
    const int n_basis = static_cast<int>(P_alpha.shape(0));
    return mulliken_spin(frame, P_alpha.data(), P_beta.data(), S.data(), n_basis);
}

} // namespace cpu
} // namespace ATV
