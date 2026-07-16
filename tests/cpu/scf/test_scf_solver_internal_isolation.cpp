/// @file tests/cpu/scf/test_scf_solver_internal_isolation.cpp
/// @brief Phase 3 — `cpu::scf::internal::*` helpers must not leak into the
///        outer Driver/Component facade.
///
/// Two enforced invariants:
///   1. **Static-assert level**: no internal helper inherits
///      `cpu::Component`. (Compile-time guarantee that they cannot be
///      registered into any outer Plan.)
///   2. **Header-include level**: no header outside
///      `include/atomiverse/cpu/scf/internal/` includes `internal/*.h`.
///      Source files (`src/...`) and tests (`tests/...`) are allowed —
///      production-side facade headers are not. (Run-time grep gate.)
///
/// Acceptance: ac-012 (second half) in
/// `.claude/specs/cpu-driver-taxonomy.acceptance.md`.

#ifdef ATV_HAS_RHF

#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/scf/internal/density_builder.h"
#include "atomiverse/cpu/scf/internal/diis_mixer.h"
#include "atomiverse/cpu/scf/internal/dual_convergence_check.h"
#include "atomiverse/cpu/scf/internal/eigen_solver.h"
#include "atomiverse/cpu/scf/internal/hamiltonian_builder.h"

using namespace ATV;

// ── Static-assert: internal helpers do not inherit Component ────

static_assert(
    !std::is_base_of_v<cpu::Component, cpu::scf::internal::DensityBuilder>,
    "internal::DensityBuilder must not inherit cpu::Component");
static_assert(!std::is_base_of_v<cpu::Component, cpu::scf::internal::DIISMixer>,
              "internal::DIISMixer must not inherit cpu::Component");
static_assert(!std::is_base_of_v<cpu::Component,
                                 cpu::scf::internal::DualConvergenceCheck>,
              "internal::DualConvergenceCheck must not inherit cpu::Component");
static_assert(
    !std::is_base_of_v<cpu::Component, cpu::scf::internal::EigenSolver>,
    "internal::EigenSolver must not inherit cpu::Component");
static_assert(
    !std::is_base_of_v<cpu::Component, cpu::scf::internal::HamiltonianBuilder>,
    "internal::HamiltonianBuilder must not inherit cpu::Component");

TEST(SCFSolverInternalIsolation, NotComponentBaseSubclasses) {
    // The static_asserts above are the load-bearing check; this test exists
    // so ctest reports a single named pass on this invariant.
    SUCCEED() << "static_assert: no internal helper inherits cpu::Component";
}

// ── Grep gate: facade-public headers must not include internal/*.h ──

namespace {

/// Walk a header tree and, for every .h/.hpp file outside the internal/
/// subdir, scan its content for `#include "atomiverse/cpu/scf/internal/`.
void collect_offenders(const std::string& dir, std::vector<std::string>& bad) {
    // Use a tiny shell-out — the test runtime does not have <filesystem>
    // recursion turned on across all compilers we support.
    std::string cmd =
        "grep -rln '#include[ \\t]*\"atomiverse/cpu/scf/internal/' " + dir
        + " 2>/dev/null";
    FILE* p = popen(cmd.c_str(), "r");
    ASSERT_NE(p, nullptr);
    char buf[1024];
    while (std::fgets(buf, sizeof(buf), p)) {
        std::string line(buf);
        if (!line.empty() && line.back() == '\n')
            line.pop_back();
        // Skip files inside the internal/ subtree itself.
        if (line.find("/internal/") != std::string::npos)
            continue;
        bad.push_back(line);
    }
    pclose(p);
}

const char* find_repo_root() {
    // Compile-time repo root injected by tests/CMakeLists.txt — ctest runs
    // from the build directory, so a cwd probe would always skip the gate.
#ifdef ATV_REPO_ROOT
    static const std::string root = ATV_REPO_ROOT;
    std::ifstream probe(root + "/CLAUDE.md");
    if (probe.good())
        return root.c_str();
#endif
    // Fallback: running from the repo root directly.
    std::ifstream probe_cwd("CLAUDE.md");
    if (probe_cwd.good())
        return ".";
    return nullptr;
}

} // namespace

TEST(SCFSolverInternalIsolation, FacadePublicHeadersDoNotIncludeInternal) {
    const char* root = find_repo_root();
    if (!root) {
        GTEST_SKIP() << "test must run from repo root for the grep gate "
                        "(no CLAUDE.md visible at cwd)";
    }

    std::vector<std::string> bad;
    collect_offenders(std::string(root) + "/include/atomiverse", bad);
    EXPECT_TRUE(bad.empty())
        << "facade-public headers must not #include cpu/scf/internal/*.h. "
        << "Offending files:\n  " << ([&] {
               std::ostringstream s;
               for (auto& b : bad)
                   s << b << "\n  ";
               return s.str();
           })();
}

#endif // ATV_HAS_RHF
