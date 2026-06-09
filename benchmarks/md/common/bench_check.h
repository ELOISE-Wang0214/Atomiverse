#pragma once
/// @file benchmarks/md/bench_check.h
/// @brief Tiny self-contained PASS/FAIL harness for the MD
///        LAMMPS-validation benchmarks.
///
/// These benchmarks were ported out of the gtest unit-test suite
/// (`tests/cpu/md/test_*_lammps_align`, `test_lammps_*_parity`). They
/// compare Atomiverse MD components against LAMMPS reference behaviour and
/// reference ensemble averages — a validation benchmark, not a portable
/// `ctest` gate. This header gives them a harness with no gtest and no
/// google-benchmark dependency, mirroring `benchmarks/qc/bench_qc_suite`:
/// each `BENCH_CASE` prints PASS / FAIL and `main()` (via `BENCH_MAIN`)
/// returns a non-zero exit code if any case fails.
///
/// Port mapping from the gtest originals:
///   TEST(Suite, Name)  -> BENCH_CASE(Suite_Name)
///   EXPECT_NEAR        -> BM_NEAR         EXPECT_EQ    -> BM_EQ
///   EXPECT_GT/LT/GE/LE -> BM_GT/LT/GE/LE  EXPECT_TRUE  -> BM_TRUE
///   EXPECT_DOUBLE_EQ   -> BM_DOUBLE_EQ    EXPECT_FALSE -> BM_FALSE
///   ASSERT_FALSE       -> BM_ASSERT_FALSE (early-returns from the case)
///   ASSERT_TRUE        -> BM_ASSERT_TRUE  (early-returns from the case)
/// Every non-aborting BM_* check accepts a trailing `<< "diagnostic"`
/// stream exactly like its gtest counterpart.

#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace bench {

/// One registered validation case.
struct Case {
    std::string name;
    std::function<void()> fn;
};

/// Per-executable case registry (one TU == one benchmark executable).
inline std::vector<Case>& registry() {
    static std::vector<Case> r;
    return r;
}

/// Failure counter for the case currently running; reset before each case.
inline int& case_failures() {
    static int n = 0;
    return n;
}

/// File-scope helper that appends a `BENCH_CASE` to the registry at
/// static-init time, preserving source order within the file.
struct Registrar {
    Registrar(const char* name, std::function<void()> fn) {
        registry().push_back({name, std::move(fn)});
    }
};

namespace detail {

template <class T>
std::string val_str(const T& v) {
    if constexpr (std::is_enum_v<T>) {
        return std::to_string(static_cast<long long>(v));
    } else {
        std::ostringstream os;
        os << v;
        return os.str();
    }
}

inline std::string loc(const char* file, int line) {
    std::string f(file);
    const auto p = f.find_last_of("/\\");
    if (p != std::string::npos)
        f = f.substr(p + 1);
    return f + ":" + std::to_string(line);
}

/// Streamed failure reporter. Returned by value from every BM_* check.
/// On a failing check it carries the diagnostic head and, once the full
/// statement (including any trailing `<< "msg"`) is built, its destructor
/// prints the line and bumps the case failure counter. On a passing check
/// it is an inert sink.
class Reporter {
  public:
    Reporter(bool failed, std::string head)
        : failed_(failed), head_(std::move(head)) {}
    Reporter(Reporter&& o) noexcept
        : failed_(o.failed_), head_(std::move(o.head_)),
          msg_(std::move(o.msg_)) {
        o.failed_ = false; // a moved-from reporter must not report
    }
    Reporter(const Reporter&) = delete;
    Reporter& operator=(const Reporter&) = delete;
    ~Reporter() {
        if (!failed_)
            return;
        ++case_failures();
        const std::string m = msg_.str();
        std::fprintf(stderr, "    FAIL: %s%s%s\n", head_.c_str(),
                     m.empty() ? "" : " — ", m.c_str());
    }
    template <class T>
    Reporter& operator<<(const T& v) {
        if (failed_)
            msg_ << v;
        return *this;
    }

  private:
    bool failed_;
    std::string head_;
    std::ostringstream msg_;
};

inline Reporter check_true(const char* file, int line, const char* expr,
                           bool ok) {
    if (ok)
        return Reporter(false, {});
    return Reporter(true, loc(file, line) + " expected true: " + expr);
}

inline Reporter check_near(const char* file, int line, const char* ge,
                           const char* we, double got, double want,
                           double tol) {
    if (std::fabs(got - want) <= tol)
        return Reporter(false, {});
    std::ostringstream h;
    h << loc(file, line) << " NEAR(" << ge << ", " << we << "): |" << got
      << " - " << want << "| = " << std::fabs(got - want) << " > tol " << tol;
    return Reporter(true, h.str());
}

inline bool double_eq(double a, double b) {
    if (a == b)
        return true;
    const double diff = std::fabs(a - b);
    const double scale = std::fmax(std::fabs(a), std::fabs(b));
    return diff <= 4.0 * std::numeric_limits<double>::epsilon() * scale;
}

inline Reporter check_double_eq(const char* file, int line, const char* ae,
                                const char* be, double a, double b) {
    if (double_eq(a, b))
        return Reporter(false, {});
    std::ostringstream h;
    h << loc(file, line) << " DOUBLE_EQ(" << ae << ", " << be << "): " << a
      << " != " << b;
    return Reporter(true, h.str());
}

enum class Op { GT, LT, GE, LE, EQ };

inline const char* op_str(Op op) {
    switch (op) {
        case Op::GT: return ">";
        case Op::LT: return "<";
        case Op::GE: return ">=";
        case Op::LE: return "<=";
        case Op::EQ: return "==";
    }
    return "?";
}

template <class A, class B>
Reporter check_cmp(const char* file, int line, const char* ae, const char* be,
                   Op op, const A& a, const B& b) {
    bool ok = false;
    switch (op) {
        case Op::GT: ok = a > b; break;
        case Op::LT: ok = a < b; break;
        case Op::GE: ok = a >= b; break;
        case Op::LE: ok = a <= b; break;
        case Op::EQ: ok = a == b; break;
    }
    if (ok)
        return Reporter(false, {});
    std::ostringstream h;
    h << loc(file, line) << " EXPECT(" << ae << ' ' << op_str(op) << ' ' << be
      << "): " << val_str(a) << ' ' << op_str(op) << ' ' << val_str(b)
      << " is false";
    return Reporter(true, h.str());
}

/// Run every registered case, print a PASS/FAIL line each, return a
/// non-zero exit code if any case failed.
inline int run_all(const char* title) {
    std::printf("# %s\n", title);
    int failed_cases = 0;
    for (auto& c : registry()) {
        case_failures() = 0;
        c.fn();
        const bool ok = case_failures() == 0;
        std::printf("  %-58s %s\n", c.name.c_str(), ok ? "PASS" : "FAIL");
        if (!ok)
            ++failed_cases;
    }
    const std::size_t n = registry().size();
    std::printf("#\n# %zu / %zu cases passed.\n", n - failed_cases, n);
    return failed_cases == 0 ? 0 : 1;
}

} // namespace detail
} // namespace bench

#define BENCH_CASE(name)                                                     \
    static void bench_case_##name();                                         \
    static const bench::Registrar bench_reg_##name{#name,                    \
                                                   bench_case_##name};       \
    static void bench_case_##name()

#define BENCH_MAIN(title)                                                    \
    int main() { return bench::detail::run_all(title); }

#define BM_TRUE(cond)                                                        \
    bench::detail::check_true(__FILE__, __LINE__, #cond, (bool)(cond))
#define BM_FALSE(cond)                                                        \
    bench::detail::check_true(__FILE__, __LINE__, "!(" #cond ")",             \
                              !(bool)(cond))
#define BM_NEAR(got, want, tol)                                              \
    bench::detail::check_near(__FILE__, __LINE__, #got, #want, (got), (want), \
                              (tol))
#define BM_DOUBLE_EQ(a, b)                                                   \
    bench::detail::check_double_eq(__FILE__, __LINE__, #a, #b, (a), (b))
#define BM_GT(a, b)                                                          \
    bench::detail::check_cmp(__FILE__, __LINE__, #a, #b,                      \
                             bench::detail::Op::GT, (a), (b))
#define BM_LT(a, b)                                                          \
    bench::detail::check_cmp(__FILE__, __LINE__, #a, #b,                      \
                             bench::detail::Op::LT, (a), (b))
#define BM_GE(a, b)                                                          \
    bench::detail::check_cmp(__FILE__, __LINE__, #a, #b,                      \
                             bench::detail::Op::GE, (a), (b))
#define BM_LE(a, b)                                                          \
    bench::detail::check_cmp(__FILE__, __LINE__, #a, #b,                      \
                             bench::detail::Op::LE, (a), (b))
#define BM_EQ(a, b)                                                          \
    bench::detail::check_cmp(__FILE__, __LINE__, #a, #b,                      \
                             bench::detail::Op::EQ, (a), (b))
#define BM_ASSERT_FALSE(cond)                                                \
    do {                                                                     \
        if ((cond)) {                                                        \
            bench::detail::check_true(__FILE__, __LINE__, "!(" #cond ")",     \
                                      false);                                \
            return;                                                          \
        }                                                                    \
    } while (0)
#define BM_ASSERT_TRUE(cond)                                                 \
    do {                                                                     \
        if (!(cond)) {                                                       \
            bench::detail::check_true(__FILE__, __LINE__, #cond, false);      \
            return;                                                          \
        }                                                                    \
    } while (0)
