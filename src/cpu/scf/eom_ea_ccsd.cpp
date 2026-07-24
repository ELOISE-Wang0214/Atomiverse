/// @file cpu/scf/eom_ea_ccsd.cpp
/// @brief Spin-orbital EOM-EA-CCSD electron-attachment energies.
///
/// Intermediates follow PySCF `cc.gintermediates` (Gauss & Stanton 1995 Table
/// III, EOM section (b)); the right-hand σ equations are Nooijen & Bartlett,
/// J. Chem. Phys. 102, 3629 (1995) Eqs. (30)–(31) as implemented in
/// `eom_gccsd.eaccsd_matvec`. Packed vector basis: r1[a] + r2[j,a>b].

#include "atomiverse/cpu/scf/eom_ea_ccsd.h"
#include "atomiverse/cpu/scf/internal/eom_ea_ccsd_pack.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/scf/internal/ccsd_pack.h"
#include "atomiverse/cpu/scf/rhf.h"

namespace ATV {
namespace cpu {
namespace scf {
namespace {

using T2 = xt::xtensor<double, 4>;
using T1 = xt::xtensor<double, 2>;
using R1 = xt::xtensor<double, 1>;  // virtual singles (nv)
using R2 = xt::xtensor<double, 3>;  // 2p1h (no, nv, nv)

/// τ_ij^ab = t2 + ½ fac P(ij)P(ab) t1 t1 (PySCF gintermediates.make_tau).
T2 make_tau(const T2& t2, const T1& t1, double fac) {
    const std::size_t no = t1.shape(0);
    const std::size_t nv = t1.shape(1);
    T2 tau = t2;
    const double c = 0.5 * fac;
    for (std::size_t i = 0; i < no; ++i)
        for (std::size_t j = 0; j < no; ++j)
            for (std::size_t a = 0; a < nv; ++a)
                for (std::size_t b = 0; b < nv; ++b) {
                    const double t1t1_ijab = c * t1(i, a) * t1(j, b);
                    const double t1t1_jiab = c * t1(j, a) * t1(i, b);
                    const double a_ij = t1t1_ijab - t1t1_jiab;
                    const double t1t1_ijba = c * t1(i, b) * t1(j, a);
                    const double t1t1_jiba = c * t1(j, b) * t1(i, a);
                    const double a_ij_ba = t1t1_ijba - t1t1_jiba;
                    tau(i, j, a, b) += a_ij - a_ij_ba;
                }
    return tau;
}

struct EomImds {
    T1 Foo, Fvv, Fov;
    T2 Wovvo, Woovv, Woooo, Wooov, Wovoo, Wvovv, Wvvvv, Wvvvo;
    T2 t2;
    T1 t1;
    int no = 0, nv = 0;
};

/// Build EOM-EE intermediates from a converged CCSD solution (gintermediates).
EomImds build_imds(const CCSDSolution& sol) {
    const int no = sol.no;
    const int nv = sol.nv;
    const std::size_t No = static_cast<std::size_t>(no);
    const std::size_t Nv = static_cast<std::size_t>(nv);
    const T1& t1 = sol.t1;
    const T2& t2 = sol.t2;
    const T2& oooo = sol.oooo;
    const T2& ooov = sol.ooov;
    const T2& oovv = sol.oovv;
    const T2& ovov = sol.ovov;
    const T2& ovvv = sol.ovvv;
    const T2& vvvv = sol.vvvv;
    const T1& foo = sol.foo;
    const T1& fvv = sol.fvv;
    const T1& fov = sol.fov;

    const T2 tau = make_tau(t2, t1, 1.0);
    const T2 tau_h = make_tau(t2, t1, 0.5);

    // --- cc_F* (section a) -------------------------------------------------
    // Fae = fvv - 0.5 fov·t1 + t1·ovvv^T - 0.5 tau_h·oovv
    T1 ccFvv = fvv;
    for (int a = 0; a < nv; ++a)
        for (int e = 0; e < nv; ++e) {
            double s = ccFvv(a, e);
            for (int m = 0; m < no; ++m)
                s -= 0.5 * fov(m, e) * t1(m, a);
            for (int m = 0; m < no; ++m)
                for (int f = 0; f < nv; ++f)
                    s += t1(m, f) * ovvv(m, a, f, e); // amef = ovvv(m,a,f,e)
            for (int m = 0; m < no; ++m)
                for (int n = 0; n < no; ++n)
                    for (int f = 0; f < nv; ++f)
                        s -= 0.5 * tau_h(m, n, a, f) * oovv(m, n, e, f);
            ccFvv(a, e) = s;
        }
    T1 ccFoo = foo;
    for (int m = 0; m < no; ++m)
        for (int i = 0; i < no; ++i) {
            double s = ccFoo(m, i);
            for (int e = 0; e < nv; ++e)
                s += 0.5 * fov(m, e) * t1(i, e);
            for (int n = 0; n < no; ++n)
                for (int e = 0; e < nv; ++e)
                    s += t1(n, e) * ooov(m, n, i, e);
            for (int n = 0; n < no; ++n)
                for (int e = 0; e < nv; ++e)
                    for (int f = 0; f < nv; ++f)
                        s += 0.5 * tau_h(i, n, e, f) * oovv(m, n, e, f);
            ccFoo(m, i) = s;
        }
    T1 ccFov = fov;
    for (int m = 0; m < no; ++m)
        for (int e = 0; e < nv; ++e) {
            double s = ccFov(m, e);
            for (int n = 0; n < no; ++n)
                for (int f = 0; f < nv; ++f)
                    s += t1(n, f) * oovv(m, n, e, f);
            ccFov(m, e) = s;
        }

    // EOM F = cc_F with extra ½ t1 contractions
    T1 Fvv = ccFvv;
    for (int a = 0; a < nv; ++a)
        for (int e = 0; e < nv; ++e) {
            double s = Fvv(a, e);
            for (int m = 0; m < no; ++m)
                s -= 0.5 * t1(m, a) * ccFov(m, e);
            Fvv(a, e) = s;
        }
    T1 Foo = ccFoo;
    for (int m = 0; m < no; ++m)
        for (int i = 0; i < no; ++i) {
            double s = Foo(m, i);
            for (int e = 0; e < nv; ++e)
                s += 0.5 * t1(i, e) * ccFov(m, e);
            Foo(m, i) = s;
        }
    T1 Fov = ccFov;

    // cc_Woooo
    T2 ccWoooo = oooo;
    for (int m = 0; m < no; ++m)
        for (int n = 0; n < no; ++n)
            for (int i = 0; i < no; ++i)
                for (int j = 0; j < no; ++j) {
                    double s = ccWoooo(m, n, i, j);
                    for (int e = 0; e < nv; ++e) {
                        s += t1(j, e) * ooov(m, n, i, e);
                        s -= t1(i, e) * ooov(m, n, j, e);
                    }
                    for (int e = 0; e < nv; ++e)
                        for (int f = 0; f < nv; ++f)
                            s += 0.25 * tau(i, j, e, f) * oovv(m, n, e, f);
                    ccWoooo(m, n, i, j) = s;
                }
    // Woooo = cc_Woooo + 0.25 tau oovv
    T2 Woooo = ccWoooo;
    for (int m = 0; m < no; ++m)
        for (int n = 0; n < no; ++n)
            for (int i = 0; i < no; ++i)
                for (int j = 0; j < no; ++j) {
                    double s = Woooo(m, n, i, j);
                    for (int e = 0; e < nv; ++e)
                        for (int f = 0; f < nv; ++f)
                            s += 0.25 * tau(i, j, e, f) * oovv(m, n, e, f);
                    Woooo(m, n, i, j) = s;
                }

    // cc_Wvvvv
    T2 ccWvvvv = vvvv;
    for (int a = 0; a < nv; ++a)
        for (int b = 0; b < nv; ++b)
            for (int e = 0; e < nv; ++e)
                for (int f = 0; f < nv; ++f) {
                    double s = ccWvvvv(a, b, e, f);
                    double tmp_baef = 0.0, tmp_abef = 0.0;
                    for (int m = 0; m < no; ++m) {
                        tmp_baef += t1(m, b) * ovvv(m, a, f, e);
                        tmp_abef += t1(m, a) * ovvv(m, b, f, e);
                    }
                    s += -tmp_baef + tmp_abef;
                    for (int m = 0; m < no; ++m)
                        for (int n = 0; n < no; ++n)
                            s += 0.25 * tau(m, n, a, b) * oovv(m, n, e, f);
                    ccWvvvv(a, b, e, f) = s;
                }
    // Wvvvv = cc_Wvvvv + 0.25 tau oovv
    T2 Wvvvv = ccWvvvv;
    for (int a = 0; a < nv; ++a)
        for (int b = 0; b < nv; ++b)
            for (int e = 0; e < nv; ++e)
                for (int f = 0; f < nv; ++f) {
                    double s = Wvvvv(a, b, e, f);
                    for (int m = 0; m < no; ++m)
                        for (int n = 0; n < no; ++n)
                            s += 0.25 * tau(m, n, a, b) * oovv(m, n, e, f);
                    Wvvvv(a, b, e, f) = s;
                }

    // cc_Wovvo: Wmbej = ovvo + t1·ovvv - t1·oovo - 0.5 t2 oovv - t1 t1 oovv
    // ovvo(m,b,e,j) = -ovov(m,b,j,e)
    T2 ccWovvo = xt::zeros<double>({No, Nv, Nv, No});
    for (int m = 0; m < no; ++m)
        for (int b = 0; b < nv; ++b)
            for (int e = 0; e < nv; ++e)
                for (int j = 0; j < no; ++j) {
                    double s = -ovov(m, b, j, e);
                    for (int f = 0; f < nv; ++f)
                        s += t1(j, f) * ovvv(m, b, e, f);
                    for (int n = 0; n < no; ++n)
                        s -= t1(n, b) * (-ooov(m, n, j, e)); // oovo = -ooov
                    for (int n = 0; n < no; ++n)
                        for (int f = 0; f < nv; ++f)
                            s -= 0.5 * t2(j, n, f, b) * oovv(m, n, e, f);
                    for (int n = 0; n < no; ++n)
                        for (int f = 0; f < nv; ++f)
                            s -= t1(j, f) * t1(n, b) * oovv(m, n, e, f);
                    ccWovvo(m, b, e, j) = s;
                }
    // Wovvo = cc_Wovvo - 0.5 t2 oovv
    T2 Wovvo = ccWovvo;
    for (int m = 0; m < no; ++m)
        for (int b = 0; b < nv; ++b)
            for (int e = 0; e < nv; ++e)
                for (int j = 0; j < no; ++j) {
                    double s = Wovvo(m, b, e, j);
                    for (int n = 0; n < no; ++n)
                        for (int f = 0; f < nv; ++f)
                            s -= 0.5 * t2(j, n, f, b) * oovv(m, n, e, f);
                    Wovvo(m, b, e, j) = s;
                }

    // Wooov = ooov + t1·oovv
    T2 Wooov = ooov;
    for (int m = 0; m < no; ++m)
        for (int n = 0; n < no; ++n)
            for (int i = 0; i < no; ++i)
                for (int e = 0; e < nv; ++e) {
                    double s = Wooov(m, n, i, e);
                    for (int f = 0; f < nv; ++f)
                        s += t1(i, f) * oovv(m, n, f, e); // mnfe with f,e
                    Wooov(m, n, i, e) = s;
                }

    // Wvovv: Wamef = -t1·oovv - ovvv.T(1,0,2,3)
    // Wamef[a,m,e,f] = -Σ_n t1(n,a) oovv(n,m,e,f) - ovvv(m,a,e,f)
    T2 Wvovv = xt::zeros<double>({Nv, No, Nv, Nv});
    for (int a = 0; a < nv; ++a)
        for (int m = 0; m < no; ++m)
            for (int e = 0; e < nv; ++e)
                for (int f = 0; f < nv; ++f) {
                    double s = -ovvv(m, a, e, f);
                    for (int n = 0; n < no; ++n)
                        s -= t1(n, a) * oovv(n, m, e, f);
                    Wvovv(a, m, e, f) = s;
                }

    // Wovoo (Wmbij) — full formula from gintermediates.Wovoo
    T2 Wovoo = xt::zeros<double>({No, Nv, No, No});
    for (int m = 0; m < no; ++m)
        for (int b = 0; b < nv; ++b)
            for (int i = 0; i < no; ++i)
                for (int j = 0; j < no; ++j) {
                    // eris.ooov.conj().T(2,3,0,1) → ooov(i,j,m,b) for real
                    double s = ooov(i, j, m, b);
                    // -Fov · t2
                    for (int e = 0; e < nv; ++e)
                        s -= Fov(m, e) * t2(i, j, b, e); // ijbe with e→ last
                    // -t1 · Woooo
                    for (int n = 0; n < no; ++n)
                        s -= t1(n, b) * Woooo(m, n, i, j);
                    // 0.5 ovvv · tau
                    for (int e = 0; e < nv; ++e)
                        for (int f = 0; f < nv; ++f)
                            s += 0.5 * ovvv(m, b, e, f) * tau(i, j, e, f);
                    // tmp1 = ooov · t2 : Σ_ne ooov(m,n,i,e) t2(j,n,b,e)  P(ij)
                    double t1a = 0.0, t1b = 0.0;
                    for (int n = 0; n < no; ++n)
                        for (int e = 0; e < nv; ++e) {
                            t1a += ooov(m, n, i, e) * t2(j, n, b, e);
                            t1b += ooov(m, n, j, e) * t2(i, n, b, e);
                        }
                    s += t1a - t1b;
                    // tmp2 = t1 · ovvo − t1 t2 oovv : P(ij)
                    // ovvo(m,b,e,j) = -ovov(m,b,j,e)
                    double u1 = 0.0, u2 = 0.0;
                    for (int e = 0; e < nv; ++e) {
                        u1 += t1(i, e) * (-ovov(m, b, j, e));
                        u2 += t1(j, e) * (-ovov(m, b, i, e));
                    }
                    for (int n = 0; n < no; ++n)
                        for (int e = 0; e < nv; ++e)
                            for (int f = 0; f < nv; ++f) {
                                u1 -= t1(i, e) * t2(n, j, b, f) *
                                      oovv(m, n, e, f);
                                u2 -= t1(j, e) * t2(n, i, b, f) *
                                      oovv(m, n, e, f);
                            }
                    s += u1 - u2;
                    Wovoo(m, b, i, j) = s;
                }

    // Wvvvo (Wabei)
    T2 Wvvvo = xt::zeros<double>({Nv, Nv, Nv, No});
    for (int a = 0; a < nv; ++a)
        for (int b = 0; b < nv; ++b)
            for (int e = 0; e < nv; ++e)
                for (int i = 0; i < no; ++i) {
                    // eris_vvvo = -ovvv.T(2,3,1,0) → -ovvv(i,e,b,a) wait:
                    // ovvv(m,a,e,f) → vvvo[e,f,a,m] = -ovvv(m,a,e,f)?
                    // eris_vvvo = -ovvv.transpose(2,3,1,0)
                    // vvvo[e,f,b,m] = -ovvv(m,b,e,f)
                    // For Wabei[a,b,e,i]: eris_vvvo[a,b,e,i] = -ovvv(i,e,b,a)
                    // transpose(2,3,1,0) of ovvv(m,a',e',f'):
                    //   out[e',f',a',m] = ovvv(m,a',e',f')
                    // so eris_vvvo[a,b,e,i] = -ovvv(i, e? wait indices)
                    // out[p,q,r,s] = ovvv(s,r,p,q) with p=e',q=f',r=a',s=m
                    // eris_vvvo[a,b,e,i] = -ovvv(i,e,a,b)  if (p,q,r,s)=(a,b,e,i)
                    // p=a,q=b,r=e,s=i → ovvv(i,e,a,b)
                    double s = -ovvv(i, e, a, b);

                    // 0.5 oovo · tau : oovo(m,n,e,i) = -ooov(m,n,i,e)
                    for (int m = 0; m < no; ++m)
                        for (int n = 0; n < no; ++n)
                            s += 0.5 * (-ooov(m, n, i, e)) * tau(m, n, a, b);
                    // -Fov · t2
                    for (int m = 0; m < no; ++m)
                        s -= Fov(m, e) * t2(m, i, a, b);
                    // -P(ab) ovvv · t2 : tmp1 = Σ_mf ovvv(m,b,e,f) t2(m,i,a,f)
                    double v1 = 0.0, v2 = 0.0;
                    for (int m = 0; m < no; ++m)
                        for (int f = 0; f < nv; ++f) {
                            v1 += ovvv(m, b, e, f) * t2(m, i, a, f);
                            v2 += ovvv(m, a, e, f) * t2(m, i, b, f);
                        }
                    s -= v1 - v2;
                    // -P(ab) (t1·ovvo − t1 t2 oovv)
                    // ovvo(m,b,e,i) = -ovov(m,b,i,e)
                    double w1 = 0.0, w2 = 0.0;
                    for (int m = 0; m < no; ++m) {
                        w1 += t1(m, a) * (-ovov(m, b, i, e));
                        w2 += t1(m, b) * (-ovov(m, a, i, e));
                    }
                    for (int m = 0; m < no; ++m)
                        for (int n = 0; n < no; ++n)
                            for (int f = 0; f < nv; ++f) {
                                w1 -= t1(m, a) * t2(n, i, b, f) *
                                      oovv(m, n, e, f);
                                w2 -= t1(m, b) * t2(n, i, a, f) *
                                      oovv(m, n, e, f);
                            }
                    s -= w1 - w2;
                    // + Wvvvv · t1
                    for (int f = 0; f < nv; ++f)
                        s += Wvvvv(a, b, e, f) * t1(i, f);
                    Wvvvo(a, b, e, i) = s;
                }

    EomImds im;
    im.Foo = std::move(Foo);
    im.Fvv = std::move(Fvv);
    im.Fov = std::move(Fov);
    im.Wovvo = std::move(Wovvo);
    im.Woovv = oovv;
    im.Woooo = std::move(Woooo);
    im.Wooov = std::move(Wooov);
    im.Wovoo = std::move(Wovoo);
    im.Wvovv = std::move(Wvovv);
    im.Wvvvv = std::move(Wvvvv);
    im.Wvvvo = std::move(Wvvvo);
    im.t1 = t1;
    im.t2 = t2;
    im.no = no;
    im.nv = nv;
    return im;
}


/// Packed EA vector size: nv + no * C(nv, 2).
int vector_size_ea(int no, int nv) {
    return nv + no * (nv * (nv - 1) / 2);
}

void unpack_ea(const std::vector<double>& v, int no, int nv, R1& r1, R2& r2) {
    const std::size_t No = static_cast<std::size_t>(no);
    const std::size_t Nv = static_cast<std::size_t>(nv);
    r1 = xt::zeros<double>({Nv});
    r2 = xt::zeros<double>({No, Nv, Nv});
    for (int a = 0; a < nv; ++a)
        r1(a) = v[static_cast<std::size_t>(a)];
    int k = nv;
    // r2[:, a>b] packed as tril_indices(nvir, -1): for each j, a, b<a
    for (int j = 0; j < no; ++j)
        for (int a = 0; a < nv; ++a)
            for (int b = 0; b < a; ++b) {
                const double x = v[static_cast<std::size_t>(k++)];
                r2(j, a, b) = x;
                r2(j, b, a) = -x;
            }
}

std::vector<double> pack_ea(const R1& r1, const R2& r2, int no, int nv) {
    std::vector<double> v(static_cast<std::size_t>(vector_size_ea(no, nv)), 0.0);
    for (int a = 0; a < nv; ++a)
        v[static_cast<std::size_t>(a)] = r1(a);
    int k = nv;
    for (int j = 0; j < no; ++j)
        for (int a = 0; a < nv; ++a)
            for (int b = 0; b < a; ++b)
                v[static_cast<std::size_t>(k++)] = r2(j, a, b);
    return v;
}

/// Antisymmetrize R2 in virtuals in place: ½ P(ab).
void antisym_r2(R2& r2, int no, int nv) {
    R2 s = xt::zeros<double>({static_cast<std::size_t>(no),
                              static_cast<std::size_t>(nv),
                              static_cast<std::size_t>(nv)});
    for (int j = 0; j < no; ++j)
        for (int a = 0; a < nv; ++a)
            for (int b = 0; b < nv; ++b)
                s(j, a, b) = 0.5 * (r2(j, a, b) - r2(j, b, a));
    r2 = s;
}

/// Right-hand EOM-EA-CCSD σ = H̄ r (Nooijen/Bartlett + PySCF eom_gccsd).
std::vector<double> matvec_ea(const EomImds& im, const std::vector<double>& vec) {
    const int no = im.no;
    const int nv = im.nv;
    R1 r1;
    R2 r2;
    unpack_ea(vec, no, nv, r1, r2);

    // Eq. (30) Hr1[a]
    R1 Hr1 = xt::zeros<double>({static_cast<std::size_t>(nv)});
    for (int a = 0; a < nv; ++a) {
        double s = 0.0;
        for (int c = 0; c < nv; ++c)
            s += im.Fvv(a, c) * r1(c);
        for (int l = 0; l < no; ++l)
            for (int d = 0; d < nv; ++d)
                s += im.Fov(l, d) * r2(l, a, d);
        for (int l = 0; l < no; ++l)
            for (int c = 0; c < nv; ++c)
                for (int d = 0; d < nv; ++d)
                    s += 0.5 * im.Wvovv(a, l, c, d) * r2(l, c, d);
        Hr1(a) = s;
    }

    // Eq. (31) Hr2[j,a,b]
    R2 Hr2 = xt::zeros<double>({static_cast<std::size_t>(no),
                                static_cast<std::size_t>(nv),
                                static_cast<std::size_t>(nv)});
    for (int j = 0; j < no; ++j)
        for (int a = 0; a < nv; ++a)
            for (int b = 0; b < nv; ++b) {
                double s = 0.0;
                for (int c = 0; c < nv; ++c)
                    s += im.Wvvvo(a, b, c, j) * r1(c);
                for (int l = 0; l < no; ++l)
                    s -= im.Foo(l, j) * r2(l, a, b);
                for (int c = 0; c < nv; ++c)
                    for (int d = 0; d < nv; ++d)
                        s += 0.5 * im.Wvvvv(a, b, c, d) * r2(j, c, d);
                double t = 0.0;
                for (int k = 0; k < no; ++k)
                    for (int l = 0; l < no; ++l)
                        for (int c = 0; c < nv; ++c)
                            for (int d = 0; d < nv; ++d)
                                t += im.Woovv(k, l, c, d) * r2(l, c, d) *
                                     im.t2(k, j, a, b);
                s -= 0.5 * t;
                Hr2(j, a, b) = s;
            }

    // P(ab) Fvv · r2 : tmp1[j,a,b] = Σ_c Fvv[a,c] r2[j,c,b]
    {
        R2 tmp1 = xt::zeros<double>({static_cast<std::size_t>(no),
                                     static_cast<std::size_t>(nv),
                                     static_cast<std::size_t>(nv)});
        for (int j = 0; j < no; ++j)
            for (int a = 0; a < nv; ++a)
                for (int b = 0; b < nv; ++b) {
                    double s = 0.0;
                    for (int c = 0; c < nv; ++c)
                        s += im.Fvv(a, c) * r2(j, c, b);
                    tmp1(j, a, b) = s;
                }
        for (int j = 0; j < no; ++j)
            for (int a = 0; a < nv; ++a)
                for (int b = 0; b < nv; ++b)
                    Hr2(j, a, b) += tmp1(j, a, b) - tmp1(j, b, a);
    }

    // P(ab) Wovvo · r2 : tmp2[j,a,b] = Σ_ld Wovvo[l,b,d,j] r2[l,a,d]
    {
        R2 tmp2 = xt::zeros<double>({static_cast<std::size_t>(no),
                                     static_cast<std::size_t>(nv),
                                     static_cast<std::size_t>(nv)});
        for (int j = 0; j < no; ++j)
            for (int a = 0; a < nv; ++a)
                for (int b = 0; b < nv; ++b) {
                    double s = 0.0;
                    for (int l = 0; l < no; ++l)
                        for (int d = 0; d < nv; ++d)
                            s += im.Wovvo(l, b, d, j) * r2(l, a, d);
                    tmp2(j, a, b) = s;
                }
        for (int j = 0; j < no; ++j)
            for (int a = 0; a < nv; ++a)
                for (int b = 0; b < nv; ++b)
                    Hr2(j, a, b) += tmp2(j, a, b) - tmp2(j, b, a);
    }

    antisym_r2(Hr2, no, nv);
    return pack_ea(Hr1, Hr2, no, nv);
}

/// Dense non-symmetric eigenproblem via the shared `linalg::eig_nonsymmetric`.
struct EigPair {
    double wr;
    double wi;
};

std::vector<EigPair> dgeev_vals(const std::vector<std::vector<double>>& H) {
    const int n = static_cast<int>(H.size());
    if (n == 0)
        return {};
    // RealMatrix is column-major, so A(i,j) = H[i][j] reproduces the old
    // hand-rolled column-major buffer bit-for-bit under dgeev.
    linalg::RealMatrix A = xt::zeros<double>(
        {static_cast<std::size_t>(n), static_cast<std::size_t>(n)});
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            A(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) =
                H[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    linalg::RealVector wr, wi;
    linalg::eig_nonsymmetric(A, wr, wi, nullptr, "eom_ea_ccsd");
    std::vector<EigPair> out(static_cast<std::size_t>(n));
    for (int j = 0; j < n; ++j) {
        out[static_cast<std::size_t>(j)].wr = wr(static_cast<std::size_t>(j));
        out[static_cast<std::size_t>(j)].wi = wi(static_cast<std::size_t>(j));
    }
    return out;
}

} // namespace

EOMEACCSDResult eom_ea_ccsd_energies(const IntegralPack& ints,
                                     const SlaterDeterminant& det,
                                     double e_scf_ha, int n_states, double tol,
                                     int max_iter) {
    if (n_states < 1)
        throw std::invalid_argument(
            "eom_ea_ccsd: n_states must be ≥ 1 (dimensionless).");

    CCSDSolution sol = ccsd_solve(ints, det, e_scf_ha, tol, max_iter);
    EOMEACCSDResult out;
    out.e_corr_ha = sol.result.e_corr;
    out.e_ccsd_ha = sol.result.e_tot;
    out.ccsd_iterations = sol.result.iterations;
    out.converged = sol.result.converged;

    if (sol.nv <= 0 || sol.no <= 0) {
        throw std::invalid_argument(
            "eom_ea_ccsd: empty occupied or virtual spin-orbital space.");
    }

    const int dim = vector_size_ea(sol.no, sol.nv);
    if (dim < 1)
        throw std::invalid_argument(
            "eom_ea_ccsd: empty EA space (need at least one virtual).");
    if (n_states > dim)
        throw std::invalid_argument(
            "eom_ea_ccsd: n_states=" + std::to_string(n_states) +
            " exceeds EOM-EA dimension " + std::to_string(dim) + ".");

    const EomImds im = build_imds(sol);

    std::vector<std::vector<double>> H(
        static_cast<std::size_t>(dim),
        std::vector<double>(static_cast<std::size_t>(dim), 0.0));
    for (int j = 0; j < dim; ++j) {
        std::vector<double> e(static_cast<std::size_t>(dim), 0.0);
        e[static_cast<std::size_t>(j)] = 1.0;
        const std::vector<double> sig = matvec_ea(im, e);
        for (int i = 0; i < dim; ++i)
            H[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                sig[static_cast<std::size_t>(i)];
    }

    auto eigs = dgeev_vals(H);

    // Real eigenvalues, ascending; de-duplicate spin degeneracy (α/β attach).
    std::vector<double> roots;
    for (const auto& ep : eigs) {
        if (std::abs(ep.wi) > 1e-8)
            continue;
        roots.push_back(ep.wr);
    }
    std::sort(roots.begin(), roots.end());
    std::vector<double> uniq;
    for (double w : roots) {
        if (!uniq.empty() && std::abs(uniq.back() - w) < 1e-8)
            continue;
        uniq.push_back(w);
    }
    if (static_cast<int>(uniq.size()) < n_states)
        throw std::runtime_error(
            "eom_ea_ccsd: found only " + std::to_string(uniq.size()) +
            " unique real EA roots; requested n_states=" +
            std::to_string(n_states));

    out.energy_ev.reserve(static_cast<std::size_t>(n_states));
    for (int k = 0; k < n_states; ++k)
        out.energy_ev.push_back(uniq[static_cast<std::size_t>(k)] *
                                units::kHartreeInEv);
    out.converged = true;
    return out;
}

EOMEACCSDResult RHF::eom_ea_ccsd(Frame& frame, int n_states) {
    const bool prev = compute_eri;
    compute_eri = true;
    run(frame);
    compute_eri = prev;
    if (!this->converged())
        throw std::invalid_argument(
            "RHF::eom_ea_ccsd: SCF reference did not converge — refusing to "
            "compute EOM-EA-CCSD on a non-stationary determinant (no "
            "fallback).");
    return eom_ea_ccsd_energies(integral_workspace_(), determinant(), energy(),
                                n_states);
}

} // namespace scf
} // namespace cpu
} // namespace ATV
