/// @file cpu/electrodes/charge_aspc.cpp
/// @brief Kolafa ASPC for electrode charge vectors.
///
/// Lagrange coefficients for order K (predict at point K+1 given history
/// 0..K, newest = K):
///   omega_k = product_{j=0..K, j!=k} ((K+1) - j) / (k - j)
///
/// Reference: Kolafa, J. Comput. Chem. 25, 335 (2004), Eqs. 7–11.

#include "atomiverse/cpu/electrodes/charge_aspc.h"

#include <stdexcept>

namespace ATV {
namespace cpu {
namespace electrodes {

ChargeASPC::ChargeASPC(int order) : K_(order) {
    if (K_ < 0)
        throw std::invalid_argument("ChargeASPC: order must be ≥ 0");
    compute_coefficients_();
}

void ChargeASPC::compute_coefficients_() {
    omega_.resize(static_cast<std::size_t>(K_ + 1));
    for (int k = 0; k <= K_; ++k) {
        double num = 1.0, den = 1.0;
        for (int j = 0; j <= K_; ++j) {
            if (j == k)
                continue;
            num *= static_cast<double>((K_ + 1) - j);
            den *= static_cast<double>(k - j);
        }
        omega_[static_cast<std::size_t>(k)] = num / den;
    }
}

void ChargeASPC::record(const std::vector<double>& Q_scf) {
    if (Q_scf.empty())
        throw std::invalid_argument("ChargeASPC::record: empty charge vector");
    if (dim_ < 0)
        dim_ = static_cast<int>(Q_scf.size());
    else if (static_cast<int>(Q_scf.size()) != dim_)
        throw std::invalid_argument(
            "ChargeASPC::record: charge vector length changed");

    history_.push_back(Q_scf);
    while (static_cast<int>(history_.size()) > K_ + 2)
        history_.pop_front();
}

bool ChargeASPC::predict(std::vector<double>& Q_predicted) const {
    if (static_cast<int>(history_.size()) < K_ + 1)
        return false;

    const int offset = static_cast<int>(history_.size()) - (K_ + 1);
    const std::size_t n = history_.back().size();
    Q_predicted.assign(n, 0.0);
    for (int k = 0; k <= K_; ++k) {
        const double w = omega_[static_cast<std::size_t>(k)];
        const auto& hk = history_[static_cast<std::size_t>(offset + k)];
        for (std::size_t i = 0; i < n; ++i)
            Q_predicted[i] += w * hk[i];
    }
    return true;
}

void ChargeASPC::correct(const std::vector<double>& Q_scf,
                         const std::vector<double>& Q_predicted,
                         std::vector<double>& Q_corrected) const {
    if (Q_scf.size() != Q_predicted.size())
        throw std::invalid_argument("ChargeASPC::correct: size mismatch");
    const double w_scf =
        static_cast<double>(K_ + 2) / static_cast<double>(2 * K_ + 3);
    const double w_pred =
        static_cast<double>(K_ + 1) / static_cast<double>(2 * K_ + 3);
    Q_corrected.resize(Q_scf.size());
    for (std::size_t i = 0; i < Q_scf.size(); ++i)
        Q_corrected[i] = w_scf * Q_scf[i] + w_pred * Q_predicted[i];
}

void ChargeASPC::reset() {
    history_.clear();
    dim_ = -1;
}

double ChargeASPC::coefficient(int k) const {
    if (k < 0 || k > K_)
        throw std::out_of_range("ChargeASPC::coefficient");
    return omega_[static_cast<std::size_t>(k)];
}

} // namespace electrodes
} // namespace cpu
} // namespace ATV
