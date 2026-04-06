/// @file cpu/scf/sto3g_data.h
/// @brief STO-3G basis set data for Z=1..8 (H through O).
///
/// Single source of truth — shared by integrals.cpp and integral_derivatives.cpp.
/// Reference: Hehre, Stewart, Pople, J. Chem. Phys. 51, 2657 (1969).
#pragma once

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Raw STO-3G shell descriptor.
/// For sp shells, both s-type and p-type contraction coefficients are stored.
struct STO3GRawShell {
    int type;           ///< 0 = pure s, 1 = sp (shared exponents, two coeff sets)
    int n_prim;         ///< Number of primitives (always 3 for STO-3G).
    double alpha[3];    ///< Exponents.
    double coeff_s[3];  ///< s-type contraction coefficients (unnormalized).
    double coeff_p[3];  ///< p-type contraction coefficients (only for type=1).
};

/// @brief STO-3G basis for one element.
struct STO3GAtomBasis {
    int n_raw_shells;          ///< 1 for H/He, 2 for Li..O.
    STO3GRawShell raw[2];      ///< Max 2 raw shells (inner + sp).
};

/// @brief Maximum supported atomic number.
static constexpr int kMaxZ = 8;

/// @brief STO-3G basis data indexed by Z (0=unused, 1=H, ..., 8=O).
static const STO3GAtomBasis sto3g_data[] = {
    // Z=0: unused
    {0, {}},
    // Z=1 (H): 1s
    {1, {
        {0, 3, {3.42525091, 0.62391373, 0.16885540},
               {0.15432897, 0.53532814, 0.44463454}, {}},
    }},
    // Z=2 (He): 1s
    {1, {
        {0, 3, {6.36242139, 1.15892300, 0.31364979},
               {0.15432897, 0.53532814, 0.44463454}, {}},
    }},
    // Z=3 (Li): 1s + 2sp
    {2, {
        {0, 3, {16.1195750, 2.93620070, 0.79465370},
               {0.15432897, 0.53532814, 0.44463454}, {}},
        {1, 3, {0.63628970, 0.14786010, 0.04808870},
               {-0.09996723, 0.39951283, 0.70011547},
               {0.15591627, 0.60768372, 0.39195739}},
    }},
    // Z=4 (Be): 1s + 2sp
    {2, {
        {0, 3, {30.1678710, 5.4951153, 1.4871927},
               {0.15432897, 0.53532814, 0.44463454}, {}},
        {1, 3, {1.3148331, 0.3055389, 0.0993707},
               {-0.09996723, 0.39951283, 0.70011547},
               {0.15591627, 0.60768372, 0.39195739}},
    }},
    // Z=5 (B): 1s + 2sp
    {2, {
        {0, 3, {48.7911130, 8.8873622, 2.4052670},
               {0.15432897, 0.53532814, 0.44463454}, {}},
        {1, 3, {2.2369561, 0.5198205, 0.1690618},
               {-0.09996723, 0.39951283, 0.70011547},
               {0.15591627, 0.60768372, 0.39195739}},
    }},
    // Z=6 (C): 1s + 2sp
    {2, {
        {0, 3, {71.6168370, 13.0450960, 3.5305122},
               {0.15432897, 0.53532814, 0.44463454}, {}},
        {1, 3, {2.9412494, 0.6834831, 0.2222899},
               {-0.09996723, 0.39951283, 0.70011547},
               {0.15591627, 0.60768372, 0.39195739}},
    }},
    // Z=7 (N): 1s + 2sp
    {2, {
        {0, 3, {99.1061690, 18.0523120, 4.8856602},
               {0.15432897, 0.53532814, 0.44463454}, {}},
        {1, 3, {3.7804559, 0.8784966, 0.2857144},
               {-0.09996723, 0.39951283, 0.70011547},
               {0.15591627, 0.60768372, 0.39195739}},
    }},
    // Z=8 (O): 1s + 2sp
    {2, {
        {0, 3, {130.7093200, 23.8088610, 6.4436083},
               {0.15432897, 0.53532814, 0.44463454}, {}},
        {1, 3, {5.0331513, 1.1695961, 0.3803890},
               {-0.09996723, 0.39951283, 0.70011547},
               {0.15591627, 0.60768372, 0.39195739}},
    }},
};

/// @brief Count expanded basis functions for an atom.
/// s-type shell -> 1 basis function. sp shell -> 4 (1s + 3p).
/// @param Z Atomic number (1..kMaxZ).
/// @return Number of basis functions for this element in STO-3G.
inline int sto3g_n_basis(int Z) {
    int n = 0;
    const auto& basis = sto3g_data[Z];
    for (int i = 0; i < basis.n_raw_shells; ++i)
        n += (basis.raw[i].type == 0) ? 1 : 4; // pure s -> 1, sp -> 1s + 3p
    return n;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
