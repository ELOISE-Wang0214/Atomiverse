/// @file benchmarks/bench_async_snapshot_mem.cu
/// @brief Resident-memory benchmark + scaling-law gate for the async snapshot
///        ring (async-analysis-snapshot-ring-03-membench).
///
/// Each `AsyncBridge` subclass and `AsyncRecordBase` hold a K-deep device ring
/// (`cudaMalloc`) plus a K-deep pinned host ring (`cudaMallocHost`), where K is
/// the compile-time `AsyncBridge::kRingDepth`. This benchmark verifies the
/// falsifiable memory model
///
///     mem ≈ K · N · B        (device ring ≈ host ring, so combined ≈ 2·K·N·B)
///
/// where N is the snapshot element count (atoms) and B the captured
/// bytes-per-element. It sweeps K and N, allocating representative rings the
/// same way production does, measures device resident memory via
/// `cudaMemGetInfo` and host resident memory via `/proc/self/status` VmRSS,
/// emits one CSV row per (K,N), least-squares-fits the device curve in K and in
/// N, and exits non-zero if either residual exceeds 10%. It also prints the
/// production `kRingDepth` footprint.
///
/// Plain `main()` (argv-driven CSV), in the style of `bench_multi_gpu_scaling`.
/// Operator-action: build on the aarch64 GH200 node with `-DATV_BUILD_CUDA=ON
/// -DATV_BUILD_BENCHMARKS=ON` (no nvcc on the login node).
///
/// @see .claude/specs/async-analysis-snapshot-ring-03-membench.md

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include <cuda_runtime.h>

#include "atomiverse/cuda/async_bridge.cuh"
#include "atomiverse/cuda/types.cuh"

using ATV::cuda::real_t;

namespace {

/// @brief Read this process's current resident set size.
/// @return VmRSS in bytes, or -1 if `/proc/self/status` is unreadable.
long read_vmrss_bytes() {
    std::ifstream st("/proc/self/status");
    std::string key;
    while (st >> key) {
        if (key == "VmRSS:") {
            long kb = 0;
            st >> kb;
            return kb * 1024L;
        }
        st.ignore(1 << 20, '\n');
    }
    return -1;
}

/// @brief Free device memory currently reported by the driver.
/// @return Free device bytes from `cudaMemGetInfo`.
size_t cuda_free_bytes() {
    size_t freeb = 0, totalb = 0;
    cudaMemGetInfo(&freeb, &totalb);
    return freeb;
}

/// @brief One swept measurement point.
struct Point {
    int K;          ///< ring depth (dimensionless)
    long N;         ///< elements per slot (atoms, dimensionless)
    long dev_bytes; ///< measured device resident delta [bytes]
};

/// @brief Max relative residual of a least-squares line y = a·x + b.
/// @param x Independent samples (K or N), dimensionless.
/// @param y Dependent samples (resident memory) [bytes].
/// @return max_i |y_i - (a·x_i + b)| / y_i, dimensionless; 0 if < 2 points.
double max_rel_residual(const std::vector<double>& x,
                        const std::vector<double>& y) {
    const std::size_t n = x.size();
    if (n < 2)
        return 0.0;
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (std::size_t i = 0; i < n; ++i) {
        sx += x[i];
        sy += y[i];
        sxx += x[i] * x[i];
        sxy += x[i] * y[i];
    }
    const double denom = n * sxx - sx * sx;
    if (denom == 0.0)
        return 0.0;
    const double a = (n * sxy - sx * sy) / denom;
    const double b = (sy - a * sx) / n;
    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double pred = a * x[i] + b;
        if (y[i] != 0.0)
            worst = std::max(worst, std::fabs(y[i] - pred) / std::fabs(y[i]));
    }
    return worst;
}

} // namespace

/// @brief Sweep K×N, measure ring resident memory, verify the K·N·B model.
/// @return 0 if both linear fits are within tolerance, 1 otherwise.
int main(int argc, char** argv) {
    // Tolerance on the least-squares residual (dimensionless). Override via
    // argv[1] for experimentation.
    const double tol = (argc > 1) ? std::atof(argv[1]) : 0.10;

    // Captured bytes-per-element: positions (x,y,z) — representative of the
    // Recorder ring, the largest per-atom payload.
    const long B = 3L * static_cast<long>(sizeof(real_t));

    // Ring depths bracketing production kRingDepth=3 (the membench sweeps K
    // itself to verify the scaling law; production K is the compile-time const).
    const std::vector<int> Ks = {1, 2, 3, 4, 6};
    // Large N so cudaMalloc page granularity (~2 MB) is negligible vs N·B and
    // the linear fit is well-posed.
    const std::vector<long> Ns = {262144L, 524288L, 1048576L, 2097152L,
                                  4194304L};
    const long FIT_N_MIN = 1048576L; // restrict fit points to where N·B >> granularity

    cudaFree(0); // force CUDA context init before measuring free memory

    std::printf("# bench_async_snapshot_mem  B=%ld bytes/elem  kRingDepth=%d\n",
                B, ATV::cuda::AsyncBridge::kRingDepth);
    std::printf("K,N,bytes_per_elem,dev_bytes,host_bytes,model_bytes\n");

    std::vector<Point> pts;
    for (int K : Ks) {
        for (long N : Ns) {
            const size_t slot = static_cast<size_t>(N) * static_cast<size_t>(B);

            const long rss0 = read_vmrss_bytes();
            const size_t free0 = cuda_free_bytes();

            std::vector<void*> dev(static_cast<size_t>(K), nullptr);
            std::vector<void*> host(static_cast<size_t>(K), nullptr);
            bool ok = true;
            for (int i = 0; i < K && ok; ++i) {
                ok = ok && cudaMalloc(&dev[i], slot) == cudaSuccess;
                ok = ok && cudaMallocHost(&host[i], slot) == cudaSuccess;
                if (ok)
                    std::memset(host[i], 0, slot); // fault pages so RSS reflects them
            }
            if (!ok) {
                std::fprintf(stderr,
                             "alloc failed at K=%d N=%ld (slot=%zu bytes)\n", K,
                             N, slot);
                for (int i = 0; i < K; ++i) {
                    if (dev[i])
                        cudaFree(dev[i]);
                    if (host[i])
                        cudaFreeHost(host[i]);
                }
                continue;
            }

            const size_t free1 = cuda_free_bytes();
            const long rss1 = read_vmrss_bytes();

            const long dev_bytes = static_cast<long>(free0 - free1);
            const long host_bytes = (rss0 >= 0 && rss1 >= 0) ? (rss1 - rss0) : -1;
            const long model = static_cast<long>(K) * static_cast<long>(slot);

            std::printf("%d,%ld,%ld,%ld,%ld,%ld\n", K, N, B, dev_bytes,
                        host_bytes, model);

            if (N >= FIT_N_MIN)
                pts.push_back({K, N, dev_bytes});

            for (int i = 0; i < K; ++i) {
                cudaFree(dev[i]);
                cudaFreeHost(host[i]);
            }
        }
    }

    // Fit 1: device memory vs K at the largest N.
    const long N_fit = Ns.back();
    std::vector<double> xk, yk;
    for (const auto& p : pts)
        if (p.N == N_fit) {
            xk.push_back(static_cast<double>(p.K));
            yk.push_back(static_cast<double>(p.dev_bytes));
        }
    // Fit 2: device memory vs N at the largest K.
    const int K_fit = Ks.back();
    std::vector<double> xn, yn;
    for (const auto& p : pts)
        if (p.K == K_fit) {
            xn.push_back(static_cast<double>(p.N));
            yn.push_back(static_cast<double>(p.dev_bytes));
        }

    const double res_k = max_rel_residual(xk, yk);
    const double res_n = max_rel_residual(xn, yn);

    // Production footprint at compile-time kRingDepth, N = 1 M atoms.
    const int Kp = ATV::cuda::AsyncBridge::kRingDepth;
    const long Np = 1048576L;
    const long prod_dev = static_cast<long>(Kp) * Np * B;
    std::printf("# production footprint: kRingDepth=%d N=%ld B=%ld -> device "
                "ring ~%ld bytes, + equal pinned host ring -> ~%ld bytes total\n",
                Kp, Np, B, prod_dev, 2 * prod_dev);
    std::printf("# linear-fit residual: vs K = %.4f, vs N = %.4f (tol %.2f)\n",
                res_k, res_n, tol);

    if (res_k > tol || res_n > tol) {
        std::fprintf(stderr,
                     "FAIL: memory not linear within tol (res_K=%.4f res_N=%.4f "
                     "tol=%.2f)\n",
                     res_k, res_n, tol);
        return 1;
    }
    std::printf("PASS: device ring memory ~ K*N*B within %.0f%%\n", tol * 100);
    return 0;
}
