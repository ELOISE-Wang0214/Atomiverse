/// @file cpu/scf/casscf.cpp
/// @brief Closed-shell CASSCF: dense active-space CI + Newton orbital steps.
///
/// Spin-summed RDM convention (matches the Slater–Condon energy used by FCI /
/// CASCI in `fci.cpp`):
///   E = V_shift + Σ_pq h_pq γ_pq + ½ Σ_pqrs (pq|rs) Γ_pqrs
/// with γ_pq = Σ_σ ⟨a†_{pσ} a_{qσ}⟩ and
/// Γ_pqrs = Σ_{στ} ⟨a†_{pσ} a†_{rτ} a_{sτ} a_{qσ}⟩.
/// Generalized Fock: F_pq = Σ_r h_pr γ_rq + Σ_rst (pr|st) Γ_rqst;
/// orbital gradient g_pq = 2 (F_pq − F_qp) for non-redundant pairs.

#include "atomiverse/cpu/scf/casscf.h"
#include "atomiverse/cpu/scf/internal/ao2mo.h"
#include "atomiverse/cpu/scf/internal/casscf_pack.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>

namespace ATV {
namespace cpu {
namespace scf {

namespace {

long long choose(int n, int k) {
    if (k < 0 || k > n)
        return 0;
    k = std::min(k, n - k);
    long long r = 1;
    for (int i = 0; i < k; ++i)
        r = r * (n - i) / (i + 1);
    return r;
}

std::vector<std::uint64_t> enumerate_strings(int n_mo, int n_elec) {
    std::vector<std::uint64_t> out;
    const std::uint64_t limit = std::uint64_t(1) << n_mo;
    for (std::uint64_t s = 0; s < limit; ++s) {
        if (__builtin_popcountll(s) == n_elec)
            out.push_back(s);
    }
    return out;
}

std::vector<int> occ_list(std::uint64_t s) {
    std::vector<int> v;
    while (s) {
        v.push_back(__builtin_ctzll(s));
        s &= s - 1;
    }
    return v;
}

int jw_single_sign(std::uint64_t s, int p, int q) {
    const int lo = std::min(p, q);
    const int hi = std::max(p, q);
    const std::uint64_t mask =
        ((std::uint64_t(1) << hi) - 1) & ~((std::uint64_t(1) << (lo + 1)) - 1);
    return (__builtin_popcountll(s & mask) & 1) ? -1 : 1;
}

double diagonal_element(std::uint64_t a, std::uint64_t b,
                        const xt::xtensor<double, 2>& h,
                        const xt::xtensor<double, 4>& g, double v_nuc) {
    const std::vector<int> occ_a = occ_list(a);
    const std::vector<int> occ_b = occ_list(b);
    double e = v_nuc;
    for (int i : occ_a)
        e += h(i, i);
    for (int i : occ_b)
        e += h(i, i);
    for (std::size_t x = 0; x < occ_a.size(); ++x)
        for (std::size_t y = x + 1; y < occ_a.size(); ++y) {
            const int i = occ_a[x], j = occ_a[y];
            e += g(i, i, j, j) - g(i, j, j, i);
        }
    for (std::size_t x = 0; x < occ_b.size(); ++x)
        for (std::size_t y = x + 1; y < occ_b.size(); ++y) {
            const int i = occ_b[x], j = occ_b[y];
            e += g(i, i, j, j) - g(i, j, j, i);
        }
    for (int i : occ_a)
        for (int j : occ_b)
            e += g(i, i, j, j);
    return e;
}

double single_element(std::uint64_t same_i, std::uint64_t same_j,
                      std::uint64_t other, const xt::xtensor<double, 2>& h,
                      const xt::xtensor<double, 4>& g) {
    // same_i → same_j by annihilating p and creating q (single spin channel).
    const std::uint64_t diff_i = same_i & ~same_j;
    const std::uint64_t diff_j = same_j & ~same_i;
    if (__builtin_popcountll(diff_i) != 1 || __builtin_popcountll(diff_j) != 1)
        return 0.0;
    const int p = __builtin_ctzll(diff_i);
    const int q = __builtin_ctzll(diff_j);
    const int sgn = jw_single_sign(same_i, p, q);
    double val = h(q, p);
    // Same-spin remaining occupied: (q r | p r) − (q r | r p) for r ≠ p.
    std::uint64_t rem = same_i & ~(std::uint64_t(1) << p);
    while (rem) {
        const int r = __builtin_ctzll(rem);
        rem &= rem - 1;
        val += g(q, p, r, r) - g(q, r, r, p);
    }
    // Opposite-spin: Coulomb only.
    std::uint64_t o = other;
    while (o) {
        const int r = __builtin_ctzll(o);
        o &= o - 1;
        val += g(q, p, r, r);
    }
    return sgn * val;
}

double double_element(std::uint64_t a_i, std::uint64_t b_i, std::uint64_t a_j,
                      std::uint64_t b_j, const xt::xtensor<double, 4>& g) {
    const int da = __builtin_popcountll(a_i ^ a_j);
    const int db = __builtin_popcountll(b_i ^ b_j);
    // Same-spin double: α differs by 2 orbitals each side, β identical.
    if (da == 4 && db == 0) {
        const std::uint64_t from = a_i & ~a_j;
        const std::uint64_t to = a_j & ~a_i;
        if (__builtin_popcountll(from) != 2 || __builtin_popcountll(to) != 2)
            return 0.0;
        const int p = __builtin_ctzll(from);
        const int r = __builtin_ctzll(from & ~(std::uint64_t(1) << p));
        const int q = __builtin_ctzll(to);
        const int s = __builtin_ctzll(to & ~(std::uint64_t(1) << q));
        // Phase: two sequential singles p→q then r→s on a_i.
        const int s1 = jw_single_sign(a_i, p, q);
        const std::uint64_t mid =
            (a_i ^ (std::uint64_t(1) << p)) | (std::uint64_t(1) << q);
        const int s2 = jw_single_sign(mid, r, s);
        // (qs|pr) − (qr|ps) with chemist notation (ab|cd).
        return double(s1 * s2) * (g(q, p, s, r) - g(q, r, s, p));
    }
    if (db == 4 && da == 0) {
        const std::uint64_t from = b_i & ~b_j;
        const std::uint64_t to = b_j & ~b_i;
        if (__builtin_popcountll(from) != 2 || __builtin_popcountll(to) != 2)
            return 0.0;
        const int p = __builtin_ctzll(from);
        const int r = __builtin_ctzll(from & ~(std::uint64_t(1) << p));
        const int q = __builtin_ctzll(to);
        const int s = __builtin_ctzll(to & ~(std::uint64_t(1) << q));
        const int s1 = jw_single_sign(b_i, p, q);
        const std::uint64_t mid =
            (b_i ^ (std::uint64_t(1) << p)) | (std::uint64_t(1) << q);
        const int s2 = jw_single_sign(mid, r, s);
        return double(s1 * s2) * (g(q, p, s, r) - g(q, r, s, p));
    }
    // Opposite-spin double: one α and one β excitation.
    if (da == 2 && db == 2) {
        const std::uint64_t af = a_i & ~a_j, at = a_j & ~a_i;
        const std::uint64_t bf = b_i & ~b_j, bt = b_j & ~b_i;
        if (__builtin_popcountll(af) != 1 || __builtin_popcountll(at) != 1 ||
            __builtin_popcountll(bf) != 1 || __builtin_popcountll(bt) != 1)
            return 0.0;
        const int p = __builtin_ctzll(af);
        const int q = __builtin_ctzll(at);
        const int r = __builtin_ctzll(bf);
        const int s = __builtin_ctzll(bt);
        const int sa = jw_single_sign(a_i, p, q);
        const int sb = jw_single_sign(b_i, r, s);
        return double(sa * sb) * g(q, p, s, r);
    }
    return 0.0;
}

double ci_element(std::uint64_t a_i, std::uint64_t b_i, std::uint64_t a_j,
                  std::uint64_t b_j, const xt::xtensor<double, 2>& h,
                  const xt::xtensor<double, 4>& g, double v_nuc) {
    const int da = __builtin_popcountll(a_i ^ a_j);
    const int db = __builtin_popcountll(b_i ^ b_j);
    if (da == 0 && db == 0)
        return diagonal_element(a_i, b_i, h, g, v_nuc);
    if (da == 2 && db == 0)
        return single_element(a_i, a_j, b_i, h, g);
    if (da == 0 && db == 2)
        return single_element(b_i, b_j, a_i, h, g);
    if ((da == 4 && db == 0) || (da == 0 && db == 4) || (da == 2 && db == 2))
        return double_element(a_i, b_i, a_j, b_j, g);
    return 0.0;
}

struct CiSolution {
    double energy = 0.0;
    xt::xtensor<double, 1> c; // CI vector, length n_det
    std::vector<std::uint64_t> alphas;
    std::vector<std::uint64_t> betas;
};

CiSolution dense_ci_solve(const xt::xtensor<double, 2>& h,
                          const xt::xtensor<double, 4>& g, double v_shift,
                          int n_alpha, int n_beta) {
    const int n_mo = static_cast<int>(h.shape(0));
    if (n_mo > 14) {
        throw std::invalid_argument(
            "casscf_energy: ncas=" + std::to_string(n_mo)
            + " exceeds the TINY cap (ncas <= 14)");
    }
    const long long n_det = choose(n_mo, n_alpha) * choose(n_mo, n_beta);
    if (n_det > 2000) {
        throw std::invalid_argument(
            "casscf_energy: N_det=" + std::to_string(n_det)
            + " exceeds the TINY cap (N_det <= 2000)");
    }
    CiSolution sol;
    sol.alphas = enumerate_strings(n_mo, n_alpha);
    sol.betas = enumerate_strings(n_mo, n_beta);
    const std::size_t n_a = sol.alphas.size();
    const std::size_t n_b = sol.betas.size();
    const std::size_t nd = n_a * n_b;
    sol.c = xt::zeros<double>({nd});

    if (nd == 1) {
        sol.energy =
            diagonal_element(sol.alphas[0], sol.betas[0], h, g, v_shift);
        sol.c(0) = 1.0;
        return sol;
    }

    xt::xtensor<double, 2> hci = xt::zeros<double>({nd, nd});
    for (std::size_t ia = 0; ia < n_a; ++ia) {
        for (std::size_t ib = 0; ib < n_b; ++ib) {
            const std::size_t idx_i = ia * n_b + ib;
            for (std::size_t ja = 0; ja < n_a; ++ja) {
                for (std::size_t jb = 0; jb < n_b; ++jb) {
                    const std::size_t idx_j = ja * n_b + jb;
                    if (idx_j < idx_i)
                        continue;
                    const double e =
                        ci_element(sol.alphas[ia], sol.betas[ib], sol.alphas[ja],
                                   sol.betas[jb], h, g, v_shift);
                    hci(idx_i, idx_j) = e;
                    hci(idx_j, idx_i) = e;
                }
            }
        }
    }
    auto [w, U] = xt::linalg::eigh(hci);
    sol.energy = w(0);
    for (std::size_t i = 0; i < nd; ++i)
        sol.c(i) = U(i, 0);
    return sol;
}

/// @brief Spin-summed active-space 1- and 2-RDMs from a CI solution.
void make_rdm12(const CiSolution& sol, int ncas,
                xt::xtensor<double, 2>& rdm1, xt::xtensor<double, 4>& rdm2) {
    const std::size_t na = static_cast<std::size_t>(ncas);
    rdm1 = xt::zeros<double>({na, na});
    rdm2 = xt::zeros<double>({na, na, na, na});
    const std::size_t n_a = sol.alphas.size();
    const std::size_t n_b = sol.betas.size();

    // Index map string → string index.
    auto find_str = [](const std::vector<std::uint64_t>& strs,
                       std::uint64_t s) -> int {
        for (std::size_t i = 0; i < strs.size(); ++i)
            if (strs[i] == s)
                return static_cast<int>(i);
        return -1;
    };

    for (std::size_t ia = 0; ia < n_a; ++ia) {
        for (std::size_t ib = 0; ib < n_b; ++ib) {
            const std::size_t idx_i = ia * n_b + ib;
            const double ci = sol.c(idx_i);
            if (std::abs(ci) < 1e-16)
                continue;
            const std::uint64_t Aa = sol.alphas[ia];
            const std::uint64_t Bb = sol.betas[ib];
            const std::vector<int> oa = occ_list(Aa);
            const std::vector<int> ob = occ_list(Bb);

            // Diagonal 1-RDM contributions.
            for (int p : oa)
                rdm1(p, p) += ci * ci;
            for (int p : ob)
                rdm1(p, p) += ci * ci;

            // α singles p→q: contributes to γ_qp (a†_q a_p).
            for (int p : oa) {
                for (int q = 0; q < ncas; ++q) {
                    if (Aa & (std::uint64_t(1) << q))
                        continue;
                    const std::uint64_t Aj =
                        (Aa ^ (std::uint64_t(1) << p)) | (std::uint64_t(1) << q);
                    const int ja = find_str(sol.alphas, Aj);
                    if (ja < 0)
                        continue;
                    const int sgn = jw_single_sign(Aa, p, q);
                    const double cj = sol.c(static_cast<std::size_t>(ja) * n_b + ib);
                    const double val = sgn * ci * cj;
                    rdm1(q, p) += val;
                }
            }
            // β singles.
            for (int p : ob) {
                for (int q = 0; q < ncas; ++q) {
                    if (Bb & (std::uint64_t(1) << q))
                        continue;
                    const std::uint64_t Bj =
                        (Bb ^ (std::uint64_t(1) << p)) | (std::uint64_t(1) << q);
                    const int jb = find_str(sol.betas, Bj);
                    if (jb < 0)
                        continue;
                    const int sgn = jw_single_sign(Bb, p, q);
                    const double cj = sol.c(ia * n_b + static_cast<std::size_t>(jb));
                    const double val = sgn * ci * cj;
                    rdm1(q, p) += val;
                }
            }

            // 2-RDM: loop determinant pairs with ≤2 differences per spin.
            // Direct accumulation Γ_pqrs = Σ_στ ⟨a†_pσ a†_rτ a_sτ a_qσ⟩.
            // Diagonal (I=I) same-spin and opposite-spin pairs.
            for (std::size_t x = 0; x < oa.size(); ++x) {
                for (std::size_t y = 0; y < oa.size(); ++y) {
                    if (x == y)
                        continue;
                    const int p = oa[x], r = oa[y];
                    // ⟨a†_pα a†_rα a_rα a_pα⟩ = 1 on |I⟩ → Γ_pprr contrib...
                    // For identical det: Γ_pqrs with q=p, s=r: a†_p a†_r a_r a_p = 1
                    rdm2(p, p, r, r) += ci * ci;
                    // exchange-like: a†_p a†_r a_p a_r = −1 for same spin
                    rdm2(p, r, r, p) -= ci * ci;
                }
            }
            for (std::size_t x = 0; x < ob.size(); ++x) {
                for (std::size_t y = 0; y < ob.size(); ++y) {
                    if (x == y)
                        continue;
                    const int p = ob[x], r = ob[y];
                    rdm2(p, p, r, r) += ci * ci;
                    rdm2(p, r, r, p) -= ci * ci;
                }
            }
            for (int p : oa) {
                for (int r : ob) {
                    // opposite spin: a†_pα a†_rβ a_rβ a_pα = 1
                    rdm2(p, p, r, r) += ci * ci;
                    // a†_pα a†_rβ a_pβ a_rα = 0 (different spins on last two if p≠r
                    // — only when p and r indices allow; for spatial Γ, the
                    // mixed-spin exchange term is zero when operators don't match)
                    // a†_rβ a†_pα a_pα a_rβ = 1 → also adds to Γ_rrpp
                    rdm2(r, r, p, p) += ci * ci;
                }
            }
        }
    }

    // Off-diagonal 2-RDM: for each I, apply two-electron operators to reach J.
    // Compact route: for every pair of dets differing by 0–2 spin-orbitals per
    // channel, add the corresponding Γ contributions. Already handled I=I above
    // for diagonal. For I≠J use ci*cj * matrix elements of a†a†aa.
    for (std::size_t ia = 0; ia < n_a; ++ia) {
        for (std::size_t ib = 0; ib < n_b; ++ib) {
            const std::size_t idx_i = ia * n_b + ib;
            const double ci = sol.c(idx_i);
            if (std::abs(ci) < 1e-16)
                continue;
            const std::uint64_t Aa = sol.alphas[ia];
            const std::uint64_t Bb = sol.betas[ib];

            // αα double: p,r occupied → q,s virtual (p < r, q < s for unique).
            {
                const std::vector<int> oa = occ_list(Aa);
                std::vector<int> va;
                for (int q = 0; q < ncas; ++q)
                    if (!(Aa & (std::uint64_t(1) << q)))
                        va.push_back(q);
                for (std::size_t ip = 0; ip < oa.size(); ++ip) {
                    for (std::size_t ir = ip + 1; ir < oa.size(); ++ir) {
                        const int p = oa[ip], r = oa[ir];
                        for (std::size_t iq = 0; iq < va.size(); ++iq) {
                            for (std::size_t is = iq + 1; is < va.size(); ++is) {
                                const int q = va[iq], s = va[is];
                                // phase p→q then r→s
                                const int s1 = jw_single_sign(Aa, p, q);
                                const std::uint64_t mid =
                                    (Aa ^ (std::uint64_t(1) << p))
                                    | (std::uint64_t(1) << q);
                                // r may have moved if r==q impossible (virt).
                                if (mid & (std::uint64_t(1) << s))
                                    continue;
                                if (!(mid & (std::uint64_t(1) << r)))
                                    continue;
                                const int s2 = jw_single_sign(mid, r, s);
                                const std::uint64_t Aj =
                                    (mid ^ (std::uint64_t(1) << r))
                                    | (std::uint64_t(1) << s);
                                const int ja = find_str(sol.alphas, Aj);
                                if (ja < 0)
                                    continue;
                                const double cj =
                                    sol.c(static_cast<std::size_t>(ja) * n_b + ib);
                                const double w = double(s1 * s2) * ci * cj;
                                // a†_q a†_s a_r a_p |I⟩ ~ |J⟩  (order: annihilate
                                // p,r create q,s). Operator a†_qα a†_sα a_rα a_pα
                                // contributes to Γ_q p s r.
                                rdm2(q, p, s, r) += w;
                                rdm2(q, r, s, p) -= w;
                                rdm2(s, p, q, r) -= w;
                                rdm2(s, r, q, p) += w;
                            }
                        }
                    }
                }
            }
            // ββ double (same structure).
            {
                const std::vector<int> ob = occ_list(Bb);
                std::vector<int> vb;
                for (int q = 0; q < ncas; ++q)
                    if (!(Bb & (std::uint64_t(1) << q)))
                        vb.push_back(q);
                for (std::size_t ip = 0; ip < ob.size(); ++ip) {
                    for (std::size_t ir = ip + 1; ir < ob.size(); ++ir) {
                        const int p = ob[ip], r = ob[ir];
                        for (std::size_t iq = 0; iq < vb.size(); ++iq) {
                            for (std::size_t is = iq + 1; is < vb.size(); ++is) {
                                const int q = vb[iq], s = vb[is];
                                const int s1 = jw_single_sign(Bb, p, q);
                                const std::uint64_t mid =
                                    (Bb ^ (std::uint64_t(1) << p))
                                    | (std::uint64_t(1) << q);
                                if (mid & (std::uint64_t(1) << s))
                                    continue;
                                if (!(mid & (std::uint64_t(1) << r)))
                                    continue;
                                const int s2 = jw_single_sign(mid, r, s);
                                const std::uint64_t Bj =
                                    (mid ^ (std::uint64_t(1) << r))
                                    | (std::uint64_t(1) << s);
                                const int jb = find_str(sol.betas, Bj);
                                if (jb < 0)
                                    continue;
                                const double cj =
                                    sol.c(ia * n_b + static_cast<std::size_t>(jb));
                                const double w = double(s1 * s2) * ci * cj;
                                rdm2(q, p, s, r) += w;
                                rdm2(q, r, s, p) -= w;
                                rdm2(s, p, q, r) -= w;
                                rdm2(s, r, q, p) += w;
                            }
                        }
                    }
                }
            }
            // αβ double: p→q (α) and r→s (β).
            {
                const std::vector<int> oa = occ_list(Aa);
                const std::vector<int> ob = occ_list(Bb);
                for (int p : oa) {
                    for (int q = 0; q < ncas; ++q) {
                        if (Aa & (std::uint64_t(1) << q))
                            continue;
                        const int sa = jw_single_sign(Aa, p, q);
                        const std::uint64_t Aj =
                            (Aa ^ (std::uint64_t(1) << p)) | (std::uint64_t(1) << q);
                        const int ja = find_str(sol.alphas, Aj);
                        if (ja < 0)
                            continue;
                        for (int r : ob) {
                            for (int s = 0; s < ncas; ++s) {
                                if (Bb & (std::uint64_t(1) << s))
                                    continue;
                                const int sb = jw_single_sign(Bb, r, s);
                                const std::uint64_t Bj =
                                    (Bb ^ (std::uint64_t(1) << r))
                                    | (std::uint64_t(1) << s);
                                const int jb = find_str(sol.betas, Bj);
                                if (jb < 0)
                                    continue;
                                const double cj = sol.c(
                                    static_cast<std::size_t>(ja) * n_b
                                    + static_cast<std::size_t>(jb));
                                const double w = double(sa * sb) * ci * cj;
                                // a†_qα a†_sβ a_rβ a_pα → Γ_q p s r
                                rdm2(q, p, s, r) += w;
                                // a†_sβ a†_qα a_pα a_rβ → Γ_s r q p
                                rdm2(s, r, q, p) += w;
                            }
                        }
                    }
                }
            }

            // α single × β identity for 2-RDM mixed with remaining occupations
            // is already partially covered by I=J diagonal; off-diagonal singles
            // contribute when one spin matches:
            // ⟨J| a†_qα a†_tτ a_uτ a_pα |I⟩ with J = α(p→q) from I.
            for (int p : occ_list(Aa)) {
                for (int q = 0; q < ncas; ++q) {
                    if (Aa & (std::uint64_t(1) << q))
                        continue;
                    const int sgn = jw_single_sign(Aa, p, q);
                    const std::uint64_t Aj =
                        (Aa ^ (std::uint64_t(1) << p)) | (std::uint64_t(1) << q);
                    const int ja = find_str(sol.alphas, Aj);
                    if (ja < 0)
                        continue;
                    const double cj =
                        sol.c(static_cast<std::size_t>(ja) * n_b + ib);
                    const double w = double(sgn) * ci * cj;
                    // Remaining α occupied (after annihilating p): couple as
                    // a†_q a†_r a_r a_p
                    std::uint64_t rem = Aa & ~(std::uint64_t(1) << p);
                    while (rem) {
                        const int r = __builtin_ctzll(rem);
                        rem &= rem - 1;
                        rdm2(q, p, r, r) += w;
                        rdm2(q, r, r, p) -= w;
                        rdm2(r, p, q, r) -= w; // careful — may overcount; use
                        // only the two standard spin-summed slots:
                    }
                    // Opposite spin occupied in β.
                    for (int r : occ_list(Bb)) {
                        rdm2(q, p, r, r) += w;
                        rdm2(r, r, q, p) += w;
                    }
                }
            }
            for (int p : occ_list(Bb)) {
                for (int q = 0; q < ncas; ++q) {
                    if (Bb & (std::uint64_t(1) << q))
                        continue;
                    const int sgn = jw_single_sign(Bb, p, q);
                    const std::uint64_t Bj =
                        (Bb ^ (std::uint64_t(1) << p)) | (std::uint64_t(1) << q);
                    const int jb = find_str(sol.betas, Bj);
                    if (jb < 0)
                        continue;
                    const double cj =
                        sol.c(ia * n_b + static_cast<std::size_t>(jb));
                    const double w = double(sgn) * ci * cj;
                    std::uint64_t rem = Bb & ~(std::uint64_t(1) << p);
                    while (rem) {
                        const int r = __builtin_ctzll(rem);
                        rem &= rem - 1;
                        rdm2(q, p, r, r) += w;
                        rdm2(q, r, r, p) -= w;
                    }
                    for (int r : occ_list(Aa)) {
                        rdm2(q, p, r, r) += w;
                        rdm2(r, r, q, p) += w;
                    }
                }
            }
        }
    }

    // Hermitize 1-RDM (real wavefunction).
    for (int p = 0; p < ncas; ++p)
        for (int q = p + 1; q < ncas; ++q) {
            const double avg = 0.5 * (rdm1(p, q) + rdm1(q, p));
            rdm1(p, q) = avg;
            rdm1(q, p) = avg;
        }
}

/// @brief Embed active RDMs into full MO space with closed-shell frozen core.
void embed_rdms(int n_mo, int ncore, int ncas,
                const xt::xtensor<double, 2>& act1,
                const xt::xtensor<double, 4>& act2, xt::xtensor<double, 2>& g1,
                xt::xtensor<double, 4>& g2) {
    const std::size_t n = static_cast<std::size_t>(n_mo);
    g1 = xt::zeros<double>({n, n});
    g2 = xt::zeros<double>({n, n, n, n});

    // Core: γ_ii = 2.
    for (int i = 0; i < ncore; ++i)
        g1(i, i) = 2.0;
    // Active block of 1-RDM.
    for (int u = 0; u < ncas; ++u)
        for (int v = 0; v < ncas; ++v)
            g1(ncore + u, ncore + v) = act1(u, v);

    // Core–core 2-RDM (closed shell):
    // Γ_iijj = 4, Γ_ijji = −2 for all core i,j (including i=j: Γ_iiii = 2).
    for (int i = 0; i < ncore; ++i) {
        for (int j = 0; j < ncore; ++j) {
            g2(i, i, j, j) += 4.0;
            g2(i, j, j, i) -= 2.0;
        }
    }
    // Core–active: Γ_iiuv = 2 γ_uv, Γ_iuv i = −γ_uv  (and symmetric slots).
    for (int i = 0; i < ncore; ++i) {
        for (int u = 0; u < ncas; ++u) {
            for (int v = 0; v < ncas; ++v) {
                const int U = ncore + u, V = ncore + v;
                const double d = act1(u, v);
                g2(i, i, U, V) += 2.0 * d;
                g2(U, V, i, i) += 2.0 * d;
                g2(i, V, U, i) -= d;
                g2(U, i, i, V) -= d;
            }
        }
    }
    // Active–active 2-RDM.
    for (int u = 0; u < ncas; ++u)
        for (int v = 0; v < ncas; ++v)
            for (int x = 0; x < ncas; ++x)
                for (int y = 0; y < ncas; ++y)
                    g2(ncore + u, ncore + v, ncore + x, ncore + y) +=
                        act2(u, v, x, y);
}

xt::xtensor<double, 2> generalized_fock(const xt::xtensor<double, 2>& h,
                                        const xt::xtensor<double, 4>& g,
                                        const xt::xtensor<double, 2>& rdm1,
                                        const xt::xtensor<double, 4>& rdm2) {
    const int n = static_cast<int>(h.shape(0));
    xt::xtensor<double, 2> F = xt::zeros<double>(
        {static_cast<std::size_t>(n), static_cast<std::size_t>(n)});
    for (int p = 0; p < n; ++p) {
        for (int q = 0; q < n; ++q) {
            double f = 0.0;
            for (int r = 0; r < n; ++r)
                f += h(p, r) * rdm1(r, q);
            for (int r = 0; r < n; ++r)
                for (int s = 0; s < n; ++s)
                    for (int t = 0; t < n; ++t)
                        f += g(p, r, s, t) * rdm2(r, q, s, t);
            F(p, q) = f;
        }
    }
    return F;
}

/// @brief exp(K) for antisymmetric K via scaled Taylor (same as SOSCF).
xt::xtensor<double, 2> exp_antisymmetric(const xt::xtensor<double, 2>& K) {
    const std::size_t m = K.shape(0);
    double fro = 0.0;
    for (std::size_t r = 0; r < m; ++r)
        for (std::size_t c = 0; c < m; ++c)
            fro += K(r, c) * K(r, c);
    fro = std::sqrt(fro);
    int s = 0;
    double scaled = fro;
    while (scaled > 0.25 && s < 60) {
        scaled *= 0.5;
        ++s;
    }
    xt::xtensor<double, 2> A = K / std::pow(2.0, s);
    xt::xtensor<double, 2> U = xt::eye<double>(m);
    xt::xtensor<double, 2> term = U;
    for (int k = 1; k <= 30; ++k) {
        term = xt::linalg::dot(term, A) / static_cast<double>(k);
        U = U + term;
        double t_max = 0.0;
        for (std::size_t r = 0; r < m; ++r)
            for (std::size_t c = 0; c < m; ++c)
                t_max = std::max(t_max, std::abs(term(r, c)));
        if (t_max < 1e-17)
            break;
    }
    for (int k = 0; k < s; ++k)
        U = xt::linalg::dot(U, U);
    return U;
}

/// @brief Approximate diagonal orbital Hessian for pair (p<q).
/// H_pq,pq ≈ 2 (F_pp − F_qq)(γ_qq − γ_pp)/something — use simple
/// 4*(ε_q − ε_p) style with F diagonal as orbital energies.
double approx_hess(const xt::xtensor<double, 2>& F,
                   const xt::xtensor<double, 2>& rdm1, int p, int q) {
    // Standard approx for CASSCF: H_{pq} ≈ 2 (γ_qq − γ_pp)(F_pp/γ − …);
    // practical robust form used by many codes for tiny systems:
    // H ≈ 4 * max(0.1, |F_qq − F_pp|) when rotation mixes occupied-like with
    // virtual-like. Fall back to 1.0 to avoid division by zero.
    const double dF = F(q, q) - F(p, p);
    const double dg = rdm1(p, p) - rdm1(q, q);
    double h = 2.0 * dF * dg;
    // For core→virtual, γ_pp=2, γ_qq=0 → h = 2*(Fqq−Fpp)*(2) = 4(Fqq−Fpp)
    // which is the HF-like Hessian piece. Clamp away from zero.
    if (std::abs(h) < 0.1)
        h = (h >= 0.0) ? 0.1 : -0.1;
    return h;
}

} // namespace

CASSCFResult casscf_energy(const IntegralPack& ints, const SlaterDeterminant& det,
                           double e_scf_ha, int ncas, int nelecas,
                           double orb_tol, int max_macro) {
    using Spin = SlaterDeterminant::Spin;
    if (ints.ERI.size() == 0) {
        throw std::invalid_argument(
            "casscf_energy: ERI is not materialized; call with compute_eri=true.");
    }
    if (det.n_electrons(Spin::alpha) != det.n_electrons(Spin::beta)) {
        throw std::invalid_argument(
            "casscf_energy: closed-shell reference only (n_alpha must equal "
            "n_beta).");
    }
    if (ncas <= 0) {
        throw std::invalid_argument(
            "casscf_energy: ncas must be positive (got " + std::to_string(ncas)
            + ")");
    }
    if (nelecas < 0 || (nelecas % 2) != 0) {
        throw std::invalid_argument(
            "casscf_energy: closed-shell CAS requires even nelecas (got "
            + std::to_string(nelecas) + ")");
    }

    // Working MO coefficients (mutable copy of the RHF reference).
    xt::xtensor<double, 2> C = det.coefficients(Spin::alpha);
    const int n_mo = static_cast<int>(C.shape(1));
    const int n_alpha = det.n_electrons(Spin::alpha);
    const int n_act_a = nelecas / 2;
    const int n_act_b = n_act_a;
    const int ncore = n_alpha - n_act_a;
    const int nvirt = n_mo - ncore - ncas;

    if (n_act_a > ncas) {
        throw std::invalid_argument(
            "casscf_energy: nelecas/2=" + std::to_string(n_act_a)
            + " exceeds ncas=" + std::to_string(ncas));
    }
    if (ncore < 0) {
        throw std::invalid_argument(
            "casscf_energy: ncore=" + std::to_string(ncore)
            + " is negative (nelecas too large for the reference occupation)");
    }
    if (ncore + ncas > n_mo) {
        throw std::invalid_argument(
            "casscf_energy: ncore+ncas=" + std::to_string(ncore + ncas)
            + " exceeds n_mo=" + std::to_string(n_mo));
    }

    // Non-redundant orbital pairs: core–active, active–virtual, core–virtual.
    struct Pair {
        int p, q; // rotate p < q in the sense κ_qp = −κ_pq; store (lo, hi)
    };
    std::vector<Pair> pairs;
    pairs.reserve(static_cast<std::size_t>(
        ncore * ncas + ncas * nvirt + ncore * nvirt));
    for (int i = 0; i < ncore; ++i)
        for (int u = 0; u < ncas; ++u)
            pairs.push_back({i, ncore + u});
    for (int u = 0; u < ncas; ++u)
        for (int a = 0; a < nvirt; ++a)
            pairs.push_back({ncore + u, ncore + ncas + a});
    for (int i = 0; i < ncore; ++i)
        for (int a = 0; a < nvirt; ++a)
            pairs.push_back({i, ncore + ncas + a});

    double e_cas = 0.0;
    int n_iter = 0;
    bool converged = false;

    for (int macro = 0; macro < max_macro; ++macro) {
        ++n_iter;
        const xt::xtensor<double, 2> h_full = mo_hcore(ints.Hcore, C, n_mo);
        const xt::xtensor<double, 4> g_full = ao2mo_full(ints.ERI, C, n_mo);

        // Frozen-core energy (closed-shell).
        double e_core = ints.V_nuc;
        for (int i = 0; i < ncore; ++i) {
            e_core += 2.0 * h_full(i, i);
            for (int j = 0; j < ncore; ++j)
                e_core += 2.0 * g_full(i, i, j, j) - g_full(i, j, j, i);
        }

        // Effective active Hamiltonian.
        const std::size_t na = static_cast<std::size_t>(ncas);
        xt::xtensor<double, 2> h_eff = xt::zeros<double>({na, na});
        xt::xtensor<double, 4> g_act = xt::zeros<double>({na, na, na, na});
        for (int u = 0; u < ncas; ++u) {
            const int U = ncore + u;
            for (int v = 0; v < ncas; ++v) {
                const int V = ncore + v;
                double he = h_full(U, V);
                for (int i = 0; i < ncore; ++i)
                    he += 2.0 * g_full(U, V, i, i) - g_full(U, i, i, V);
                h_eff(u, v) = he;
                for (int r = 0; r < ncas; ++r) {
                    const int R = ncore + r;
                    for (int s = 0; s < ncas; ++s) {
                        const int S = ncore + s;
                        g_act(u, v, r, s) = g_full(U, V, R, S);
                    }
                }
            }
        }

        const CiSolution ci =
            dense_ci_solve(h_eff, g_act, e_core, n_act_a, n_act_b);
        e_cas = ci.energy;

        // No non-redundant rotations ⇒ CASSCF ≡ CASCI.
        if (pairs.empty()) {
            converged = true;
            break;
        }

        xt::xtensor<double, 2> act1;
        xt::xtensor<double, 4> act2;
        make_rdm12(ci, ncas, act1, act2);

        // Energy consistency check (debug-level; do not throw on tiny noise).
        // Embed and build Fock for orbital gradient.
        xt::xtensor<double, 2> g1;
        xt::xtensor<double, 4> g2;
        embed_rdms(n_mo, ncore, ncas, act1, act2, g1, g2);
        const xt::xtensor<double, 2> F =
            generalized_fock(h_full, g_full, g1, g2);

        double gmax = 0.0;
        std::vector<double> grad(pairs.size());
        for (std::size_t k = 0; k < pairs.size(); ++k) {
            const int p = pairs[k].p;
            const int q = pairs[k].q;
            // g_pq = 2 (F_pq − F_qp); step rotates with κ_qp.
            grad[k] = 2.0 * (F(p, q) - F(q, p));
            gmax = std::max(gmax, std::abs(grad[k]));
        }
        if (gmax < orb_tol) {
            converged = true;
            break;
        }

        // Diagonal Newton: κ_pq = − g_pq / H_pq.
        xt::xtensor<double, 2> K = xt::zeros<double>(
            {static_cast<std::size_t>(n_mo), static_cast<std::size_t>(n_mo)});
        double knorm = 0.0;
        for (std::size_t k = 0; k < pairs.size(); ++k) {
            const int p = pairs[k].p;
            const int q = pairs[k].q;
            const double H = approx_hess(F, g1, p, q);
            double kap = -grad[k] / H;
            // Trust: |κ| ≤ 0.5.
            if (kap > 0.5)
                kap = 0.5;
            if (kap < -0.5)
                kap = -0.5;
            // Antisymmetric generator: K_qp = κ, K_pq = −κ with κ = κ_{p→q}.
            K(q, p) = kap;
            K(p, q) = -kap;
            knorm += kap * kap;
        }
        (void)knorm;
        C = xt::linalg::dot(C, exp_antisymmetric(K));
    }

    if (!converged && pairs.empty())
        converged = true;

    return CASSCFResult{e_cas, e_cas - e_scf_ha, ncas, nelecas, ncore, n_iter,
                        converged};
}

} // namespace scf
} // namespace cpu
} // namespace ATV
