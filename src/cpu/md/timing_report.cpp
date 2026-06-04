/// @file cpu/md/timing_report.cpp
/// @brief Cross-rank reduction + rendering of TimerRegistry data.

#include "atomiverse/cpu/md/timing_report.h"

#include <array>
#include <string>
#include <vector>

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

namespace {
constexpr std::size_t kN = StepTiming::kNumCategories; // 6
} // namespace

ReducedReport reduce_timing(const TimerRegistry& reg) {
    ReducedReport rep;
    rep.steps = reg.steps();

    // Local per-category totals (+ loop wall) packed into one buffer.
    const std::vector<TimerStat> rollup = reg.rollup_by_category();
    std::array<double, kN + 1> local{};
    for (std::size_t c = 0; c < kN; ++c)
        local[c] = (c < rollup.size()) ? rollup[c].total_s : 0.0;
    local[kN] = reg.loop_seconds();

    std::array<double, kN + 1> vmin = local, vmax = local, vsum = local;
    int nl = reg.nlocal() < 0 ? 0 : reg.nlocal();
    int nl_min = nl, nl_max = nl;
    long long nl_sum = nl;
    int ranks = 1;

#ifdef ATV_HAS_MPI
    int inited = 0;
    MPI_Initialized(&inited);
    if (inited) {
        MPI_Comm_size(MPI_COMM_WORLD, &ranks);
        if (ranks > 1) {
            MPI_Allreduce(local.data(), vmin.data(), kN + 1, MPI_DOUBLE,
                          MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(local.data(), vmax.data(), kN + 1, MPI_DOUBLE,
                          MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(local.data(), vsum.data(), kN + 1, MPI_DOUBLE,
                          MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&nl, &nl_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&nl, &nl_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            long long nll = nl;
            MPI_Allreduce(&nll, &nl_sum, 1, MPI_LONG_LONG, MPI_SUM,
                          MPI_COMM_WORLD);
        }
    }
#endif

    rep.ranks = ranks;
    rep.loop_min = vmin[kN];
    rep.loop_avg = vsum[kN] / ranks;
    rep.loop_max = vmax[kN];
    rep.nlocal_min = nl_min;
    rep.nlocal_max = nl_max;
    rep.nlocal_avg = static_cast<double>(nl_sum) / ranks;

    for (std::size_t c = 0; c < kN; ++c) {
        SectionStat s;
        s.min_s = vmin[c];
        s.avg_s = vsum[c] / ranks;
        s.max_s = vmax[c];
        s.pct_var = s.avg_s > 0.0 ? 100.0 * (s.max_s - s.avg_s) / s.avg_s : 0.0;
        s.pct_total =
            rep.loop_avg > 0.0 ? 100.0 * s.avg_s / rep.loop_avg : 0.0;
        rep.sections[c] = s;
    }
    return rep;
}

void emit_timing_table(const ReducedReport& rep, std::string_view stage) {
    std::string title = "Timing breakdown [" + std::string(stage) + "]  ranks="
                        + std::to_string(rep.ranks)
                        + "  steps=" + std::to_string(rep.steps);
    report::Table t(title,
                    {{"Section", 9, report::Align::Left},
                     {"min(s)"},
                     {"avg(s)"},
                     {"max(s)"},
                     {"%varavg"},
                     {"%total"}});
    for (std::size_t c = 0; c < kN; ++c) {
        const auto& s = rep.sections[c];
        t.row({std::string(StepTiming::kNames[c]), report::fmt(s.min_s),
               report::fmt(s.avg_s), report::fmt(s.max_s),
               report::fmt(s.pct_var, 1), report::fmt(s.pct_total, 1)});
    }
    t.row({"Loop", report::fmt(rep.loop_min), report::fmt(rep.loop_avg),
           report::fmt(rep.loop_max), "", "100.0"});
    t.emit();

    // Load-balance line as a second one-row table.
    report::Table lb("Load balance (owned atoms/rank)",
                     {{"nlocal_min", 12}, {"nlocal_avg", 12}, {"nlocal_max", 12}});
    lb.row({std::to_string(rep.nlocal_min), report::fmt(rep.nlocal_avg, 1),
            std::to_string(rep.nlocal_max)});
    lb.emit();
}

void emit_timing_csv(const ReducedReport& rep, report::Csv& csv,
                     std::string_view stage) {
    // Schema identical to StepTiming::csv_header():
    // stage,ranks,steps,total_s,nlocal_min,nlocal_max,nlocal_avg,<6>_s,<6>_pct
    std::string r = std::string(stage) + "," + std::to_string(rep.ranks) + ","
                    + std::to_string(rep.steps) + ","
                    + report::fmt(rep.loop_avg, 6) + ","
                    + std::to_string(rep.nlocal_min) + ","
                    + std::to_string(rep.nlocal_max) + ","
                    + report::fmt(rep.nlocal_avg, 6);
    for (std::size_t c = 0; c < kN; ++c)
        r += "," + report::fmt(rep.sections[c].avg_s, 6);
    for (std::size_t c = 0; c < kN; ++c)
        r += "," + report::fmt(rep.sections[c].pct_total, 6);
    csv.row(r);
}

} // namespace md
} // namespace cpu
} // namespace ATV
