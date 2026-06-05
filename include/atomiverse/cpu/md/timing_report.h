/// @file cpu/md/timing_report.h
/// @brief Cross-rank reduction + rendering of TimerRegistry data into a
/// LAMMPS-`timing breakdown`-comparable table (console) and CSV (metrics file).
///
/// Lives in `cpu/md` because it knows the 6 MD timing categories and their
/// names (`StepTiming`). The registry/data model is domain-agnostic
/// (`core/timer.h`); output is routed through the log module (`core/report.h`).
#pragma once

#include <array>
#include <cstdint>
#include <string_view>

#include "atomiverse/core/report.h"
#include "atomiverse/core/timer.h"
#include "atomiverse/cpu/md/step_timing.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief Per-section cross-rank stats (mirrors LAMMPS min/avg/max + %var).
struct SectionStat {
    double min_s = 0.0;
    double avg_s = 0.0;
    double max_s = 0.0;
    double pct_var = 0.0;   ///< 100*(max-avg)/avg -- load imbalance
    double pct_total = 0.0; ///< 100*avg/loop_avg
};

/// @brief Reduced timing report for one stage/run.
struct ReducedReport {
    std::array<SectionStat, StepTiming::kNumCategories> sections{};
    double loop_min = 0.0, loop_avg = 0.0, loop_max = 0.0;
    int nlocal_min = 0, nlocal_max = 0;
    double nlocal_avg = 0.0;
    std::uint64_t steps = 0;
    int ranks = 1;
};

/// @brief Reduce a registry's per-category totals + loop wall + nlocal across
/// ranks (MPI_COMM_WORLD when MPI is initialized; single-rank otherwise) into
/// a `ReducedReport`. Uses the batched min/max/sum allreduce idiom.
ReducedReport reduce_timing(const TimerRegistry& reg);

/// @brief Emit the human-readable per-section table (console, rank 0) through
/// the report logger.
void emit_timing_table(const ReducedReport& rep, std::string_view stage);

/// @brief Emit one CSV row (schema identical to StepTiming::csv_header()) to
/// the given Csv sink (metrics file, rank 0).
void emit_timing_csv(const ReducedReport& rep, report::Csv& csv,
                     std::string_view stage);

/// @brief Diagnostic: emit per-component (per-name) max-over-ranks wall time —
/// splits the category rollup (e.g. Comm) into its constituent components
/// (halo_forward / halo_reverse / atom_migration). Rank 0 only; intended for
/// scaling investigations. Component registration order is identical across
/// ranks, so the per-component totals reduce elementwise.
void emit_per_component(const TimerRegistry& reg, std::string_view stage);

} // namespace md
} // namespace cpu
} // namespace ATV
