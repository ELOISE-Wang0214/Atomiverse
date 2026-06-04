/// @file cpu/md/step_timing.h
/// @brief MD timing categories + component-name classifier + CSV schema.
///
/// The LAMMPS-`timing breakdown`-aligned category taxonomy used by the timer
/// subsystem (`core/timer.h`) and its renderer (`cpu/md/timing_report.h`).
/// Accumulation now lives in `ATV::TimerRegistry`; this header is the single
/// source of truth for the category set, their names, the name→category
/// classifier, and the machine CSV column schema. Pure metadata; no state.
#pragma once

#include <array>
#include <string>
#include <string_view>

namespace ATV {
namespace cpu {
namespace md {

/// @brief Per-component MD timing taxonomy (mirrors LAMMPS sections).
struct StepTiming {
    enum Category { Pair, Bonded, Neigh, Comm, Reduce, Other, kNumCategories };

    static constexpr std::array<std::string_view, kNumCategories> kNames{
        "Pair", "Bonded", "Neigh", "Comm", "Reduce", "Other"};

    /// @brief Map a component `name()` to a timing category.
    ///
    /// Matches the concrete CPU MD component names: pair `lj`/`morse`;
    /// bonded `harmonic_bond`/`harmonic_angle`/`multi_harmonic_dihedral`;
    /// neighbor `neighbor_list`/cell list; communication `halo_forward`/
    /// `halo_reverse`/`atom_migration`/`rcb_balancer`; reduction
    /// `reduce_energy`. Everything else (integrators, thermostats, barostat,
    /// recorders) is `Other`.
    static Category classify(std::string_view n) {
        auto has = [&](std::string_view s) {
            return n.find(s) != std::string_view::npos;
        };
        if (has("bond") || has("angle") || has("dihedral"))
            return Bonded;
        if (has("neighbor") || has("cell"))
            return Neigh;
        if (has("halo") || has("migration") || has("rcb"))
            return Comm;
        if (has("reduce"))
            return Reduce;
        if (has("lj") || has("morse") || has("pair"))
            return Pair;
        return Other;
    }

    /// @brief CSV column schema for the machine-readable per-stage breakdown.
    /// Kept here as the single source of truth; emitted by
    /// `cpu::md::emit_timing_csv` (`cpu/md/timing_report.h`).
    static std::string csv_header() {
        return "stage,ranks,steps,total_s,nlocal_min,nlocal_max,nlocal_avg,"
               "Pair_s,Bonded_s,Neigh_s,Comm_s,Reduce_s,Other_s,"
               "Pair_pct,Bonded_pct,Neigh_pct,Comm_pct,Reduce_pct,Other_pct";
    }
};

} // namespace md
} // namespace cpu
} // namespace ATV
