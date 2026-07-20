/// @file cpu/scf/scf_result.h
/// @brief Public SCF **results** and light enums — no IntegralPack, no SCFConfig.
///
/// User code: construct `RHF`/`UHF`/`RKS`/`UKS`, call `run(frame)`, read
/// `last_result()`. Config bags and integral packs are private.
#pragma once

#include <vector>
#include <xtensor/containers/xtensor.hpp>

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Per-iteration convergence record.
struct IterationRecord {
    int iter = 0;
    double energy = 0; ///< [Hartree]
    double residual = 0;
    bool diis_active = false;
};

enum class MeanFieldType { HF, KS };
enum class SpinType { Restricted, Unrestricted };
enum class InitGuess { Hcore, MINAO };

/// @brief Optional density-matrix initial guess.
struct SCFInitialGuess {
    xt::xtensor<double, 2> P_alpha;
    xt::xtensor<double, 2> P_beta;
    bool has_alpha = false;
    bool has_beta = false;
};

/// @brief Converged (or last) mean-field result.
struct SCFResult {
    bool converged = false;
    int iterations = 0;
    double energy = 0; ///< Total energy [Hartree]
    double residual = 0;

    double one_electron_energy = 0;
    double coulomb_energy = 0;
    double exchange_energy = 0;
    double xc_energy = 0;
    double nuclear_repulsion = 0;

    int n_alpha = 0;
    int n_beta = 0;
    SpinType spin = SpinType::
        Restricted; ///< Spin treatment used (restricted/unrestricted).
    double spin_contamination = 0;

    xt::xtensor<double, 1> eigenvalues_alpha;
    xt::xtensor<double, 1> eigenvalues_beta;
    xt::xtensor<double, 2> P_alpha;
    xt::xtensor<double, 2> P_beta;
    xt::xtensor<double, 2> C_alpha;
    xt::xtensor<double, 2> C_beta;

    std::vector<IterationRecord> trace;

    double integrated_n_alpha = 0;
    double integrated_n_beta = 0;
};

/// @brief Analytic / FD nuclear gradient result (public forces API).
struct GradientResult {
    xt::xtensor<double, 2> forces; ///< [n_atoms, 3] eV/Å
    double energy = 0;             ///< eV
};

} // namespace scf
} // namespace cpu
} // namespace ATV
