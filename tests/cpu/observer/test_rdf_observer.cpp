/// @file tests/cpu/observer/test_rdf_observer.cpp
/// @brief Tests for `RDF` async observer (Phase 1.9 port).
///
/// Verifies:
///   - registration via `sim.add(rdf)` (uses cpu::md::AsyncObserver default binding)
///   - histogram accumulates across snapshots
///   - finalize() writes the configured output file
///   - first g(r) peak position matches FCC nearest-neighbor distance
///   - frequency() respected

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/observer/rdf.h"

#include "cpu/driver_test_helpers.h"

using namespace ATV;
using atv_test::NoOpForce;

namespace {

cpu::Frame make_fcc_argon(int n_cells, double a) {
    int n = n_cells * n_cells * n_cells * 4;
    cpu::Frame frame(n);
    double basis[4][3] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.5, 0.0, 0.5},
        {0.0, 0.5, 0.5},
    };
    int idx = 0;
    for (int ix = 0; ix < n_cells; ++ix)
        for (int iy = 0; iy < n_cells; ++iy)
            for (int iz = 0; iz < n_cells; ++iz)
                for (int b = 0; b < 4; ++b) {
                    frame.atoms.x[idx] = (ix + basis[b][0]) * a;
                    frame.atoms.y[idx] = (iy + basis[b][1]) * a;
                    frame.atoms.z[idx] = (iz + basis[b][2]) * a;
                    frame.atoms.mass[idx] = 39.948;
                    frame.atoms.type[idx] = 18;
                    ++idx;
                }
    frame.box = Mat3::cubic(n_cells * a);
    return frame;
}

struct ScopedFile {
    std::string path;
    explicit ScopedFile(std::string p) : path(std::move(p)) {
        std::remove(path.c_str());
    }
    ~ScopedFile() { std::remove(path.c_str()); }
};

} // namespace

// ── histogram accumulates across snapshots ──────────────────────────────

TEST(RDFObserverAsync, AccumulatesAcrossSnapshots) {
    auto frame =
        make_fcc_argon(3, 5.26); // 108 atoms, lattice 5.26 Å (Argon-ish)
    ScopedFile out("test_rdf_out_accum.dat");

    RDF rdf{6.0, 200, /*freq=*/0, out.path};
    {
        cpu::md::MD sim;
        sim.add(cpu::md::VelocityVerlet{1e-6});
        sim.add(NoOpForce{});
        sim.add(rdf);
        sim.build();
        sim.run(
            frame,
            5); // steps 0..4 → publisher fires at 1,2,3,4 (freq=0 = every step >0)
    }

    EXPECT_GE(rdf.n_samples(), 1)
        << "RDF should have received at least one snapshot";
    long long total = 0;
    for (auto h : rdf.histogram())
        total += h;
    EXPECT_GT(total, 0) << "expected nonzero pair counts";
}

// ── finalize writes the output file with valid g(r) format ──────────────

TEST(RDFObserverAsync, FinalizeWritesOutputFile) {
    auto frame = make_fcc_argon(3, 5.26);
    ScopedFile out("test_rdf_out_finalize.dat");

    RDF rdf{6.0, 100, /*freq=*/0, out.path};
    {
        cpu::md::MD sim;
        sim.add(cpu::md::VelocityVerlet{1e-6});
        sim.add(NoOpForce{});
        sim.add(rdf);
        sim.build();
        sim.run(frame, 5);
    }

    std::ifstream in(out.path);
    ASSERT_TRUE(in.good()) << "output file '" << out.path
                           << "' missing after finalize";

    std::string header;
    std::getline(in, header);
    EXPECT_NE(header.find("g(r)"), std::string::npos)
        << "header missing 'g(r)' marker -- got: " << header;

    int n_lines = 0;
    double r, gr;
    while (in >> r >> gr)
        ++n_lines;
    EXPECT_EQ(n_lines, 100) << "should have one line per histogram bin";
}

// ── first peak in g(r) sits near FCC nearest-neighbor distance ──────────

TEST(RDFObserverAsync, FirstPeakNearFCCNearestNeighbor) {
    const double a_lat = 5.26;
    const double r_nn = a_lat / std::sqrt(2.0); // ≈ 3.72 Å for FCC

    auto frame = make_fcc_argon(4, a_lat); // 256 atoms for better statistics
    ScopedFile out("test_rdf_out_peak.dat");

    const int n_bins = 200;
    const double r_max = 6.0;
    const double bin_width = r_max / n_bins;

    RDF rdf{r_max, n_bins, /*freq=*/0, out.path};
    {
        cpu::md::MD sim;
        sim.add(cpu::md::VelocityVerlet{1e-6});
        sim.add(NoOpForce{});
        sim.add(rdf);
        sim.build();
        sim.run(frame, 3); // 2 firing snapshots
    }

    // Find first nonzero bin away from r=0.
    const auto& hist = rdf.histogram();
    int peak_bin = -1;
    long long peak_count = 0;
    // Skip the lowest few bins (no atoms within < 1 Å of each other on FCC).
    const int start_bin = static_cast<int>(2.0 / bin_width);
    for (int i = start_bin; i < n_bins; ++i) {
        if (hist[i] > peak_count) {
            peak_count = hist[i];
            peak_bin = i;
        }
    }
    ASSERT_GE(peak_bin, 0);

    const double r_peak = (peak_bin + 0.5) * bin_width;
    EXPECT_NEAR(r_peak, r_nn, 0.10)
        << "expected first g(r) peak near FCC NN distance " << r_nn
        << " Å, got " << r_peak << " Å";
}

// ── frequency() stride is respected ─────────────────────────────────────

TEST(RDFObserverAsync, FrequencyRespected) {
    auto frame = make_fcc_argon(3, 5.26);
    ScopedFile out("test_rdf_out_freq.dat");

    RDF rdf{6.0, 100, /*freq=*/3, out.path};
    {
        cpu::md::MD sim;
        sim.add(cpu::md::VelocityVerlet{1e-6});
        sim.add(NoOpForce{});
        sim.add(rdf);
        sim.build();
        sim.run(frame, 11); // step > 0 && step % 3 == 0 → 3, 6, 9 → 3 samples
    }

    EXPECT_EQ(rdf.n_samples(), 3)
        << "frequency=3 over 11 steps should fire at steps 3/6/9";
}
