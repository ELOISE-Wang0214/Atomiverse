/// @file tests/cuda/async/test_async_record_perf.cu
/// @brief 25× recorder-frequency wall-clock overhead gate.
///
/// Runs `polyethylene_tension --backend cuda --mode stage:2 --scale 0.005`
/// with `--print-every 25` (control) and `--print-every 1` (high-frequency)
/// five times each, takes the median wall-clock per configuration, and
/// asserts `wall_high / wall_low - 1 < 0.05`. Failure indicates a missed
/// sync somewhere on the AsyncRecordBase hot path.

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <vector>

#ifndef ATV_PE_TENSION_BIN
#error "ATV_PE_TENSION_BIN must be defined to the polyethylene_tension binary path"
#endif

namespace {

/// @brief Time one polyethylene_tension run end-to-end.
/// @param print_every Value of `--print-every`.
/// @return Wall-clock seconds.
double time_pe_tension(int print_every) {
    static int run_id = 0;
    ++run_id;
    const std::string out_dir =
        (std::filesystem::temp_directory_path()
         / ("atv_pe_perf_" + std::to_string(print_every) + "_"
            + std::to_string(run_id)))
            .string();
    std::filesystem::create_directories(out_dir);
    const std::string cmd = std::string(ATV_PE_TENSION_BIN)
                            + " --backend cuda --mode stage:2 --scale 0.005"
                            + " --seed 12345 --print-every "
                            + std::to_string(print_every)
                            + " --quiet --out-dir " + out_dir
                            + " > /dev/null 2>&1";
    auto t0 = std::chrono::steady_clock::now();
    int rc = std::system(cmd.c_str());
    auto t1 = std::chrono::steady_clock::now();
    EXPECT_EQ(rc, 0) << "polyethylene_tension exited non-zero: " << rc;
    return std::chrono::duration<double>(t1 - t0).count();
}

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

} // namespace

TEST(CudaAsyncRecordPerf, FrequencyOverheadUnder5Percent) {
    constexpr int kReps = 5;
    std::vector<double> wall_low, wall_high;
    wall_low.reserve(kReps);
    wall_high.reserve(kReps);
    // Interleave the runs so transient host load is averaged across the
    // two configurations.
    for (int r = 0; r < kReps; ++r) {
        wall_low.push_back(time_pe_tension(25));
        wall_high.push_back(time_pe_tension(1));
    }
    const double med_low = median(wall_low);
    const double med_high = median(wall_high);
    const double overhead = med_high / med_low - 1.0;
    std::printf(
        "[perf] median wall: print-every 25 = %.3fs, print-every 1 = %.3fs, "
        "overhead = %.2f%%\n",
        med_low, med_high, 100.0 * overhead);
    EXPECT_LT(overhead, 0.05)
        << "25x recorder frequency overhead is " << 100.0 * overhead
        << "% (limit 5%) — async hot path is being blocked";
}
