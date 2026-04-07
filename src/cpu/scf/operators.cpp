/// @file cpu/scf/operators.cpp
/// @brief Two-electron operator builders: J (Coulomb), K (exchange), G (fused).

#include "atomiverse/cpu/scf/operators.h"

namespace ATV {
namespace cpu {
namespace scf {

void build_coulomb(const RealMatrix& P_total,
                   const xt::xtensor<double, 4>& ERI, int n,
                   RealMatrix& J) {
    J = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) {
            double val = 0;
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l)
                    val += P_total(k, l) * ERI(i, j, k, l);
            J(i, j) = val;
            J(j, i) = val;
        }
}

void build_exchange(const RealMatrix& P_spin,
                    const xt::xtensor<double, 4>& ERI, int n,
                    RealMatrix& K) {
    K = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) {
            double val = 0;
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l)
                    val += P_spin(k, l) * ERI(i, k, j, l);
            K(i, j) = val;
            K(j, i) = val;
        }
}

void build_coulomb_exchange_fused(const RealMatrix& P_total,
                                  const xt::xtensor<double, 4>& ERI, int n,
                                  RealMatrix& G) {
    G = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double val = 0;
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l)
                    val += P_total(k, l) * (ERI(i, j, k, l) - double(0.5) * ERI(i, k, j, l));
            G(i, j) = val;
        }
}

} // namespace scf
} // namespace cpu
} // namespace ATV
