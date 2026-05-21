#ifdef ATV_HAS_RHF

#include <benchmark/benchmark.h>

#include <array>
#include <vector>

#include "atomiverse/core/constants.h"
#include "atomiverse/core/log.h"
#include "atomiverse/cpu/scf/basis_data/sto3g.h"
#include "atomiverse/cpu/scf/native_integrals.h"

#ifdef ATV_HAS_LIBCINT
extern "C" {
#include <cint.h>
#include <cint_funcs.h>
}
#undef atm
#undef bas
#endif

using namespace ATV;
using namespace ATV::cpu::scf;

namespace {

constexpr double kA2B = ATV::kAngstromToBohr;

struct AtomSpec {
    int Z;
    std::array<double, 3> center_bohr;
};

struct ShellSpec {
    int atom_of;
    int Z;
    int raw_shell;
    int l;
    int cart;
};

struct ShellHolder {
    std::array<double, 3> alpha{};
    std::array<double, 3> coeff{};
    native::Shell shell{};
};

// Map the legacy (raw_shell, l) selector onto the generic STO-3G basis_data
// table. The old sto3g_data layout fused s + sp shells; basis_data::sto3g()
// emits each l separately, so the ordering is [1s, 2s, 2p, ...].
inline int sto3g_basis_data_index(int raw_shell, int l) {
    if (raw_shell == 0 && l == 0)
        return 0; // pure-s / 1s
    if (raw_shell == 1 && l == 0)
        return 1; // 2s from sp
    if (raw_shell == 1 && l == 1)
        return 2; // 2p from sp
    return -1;
}

ShellHolder make_native_shell(const AtomSpec& atom, const ShellSpec& spec) {
    const auto& atom_shells = basis_data::sto3g()[spec.Z].shells;
    const int idx = sto3g_basis_data_index(spec.raw_shell, spec.l);
    const auto& sh = atom_shells[idx];

    ShellHolder holder;
    const int np = static_cast<int>(sh.alpha.size());
    double (*norm_at)(double) =
        (spec.l == 0) ? native::phys_norm : native::phys_norm_p;
    for (int p = 0; p < np; ++p) {
        holder.alpha[p] = sh.alpha[p];
        holder.coeff[p] = sh.coeff[p] * norm_at(sh.alpha[p]);
    }

    holder.shell.n_prim = np;
    holder.shell.alpha = holder.alpha.data();
    holder.shell.coeff = holder.coeff.data();
    holder.shell.center[0] = atom.center_bohr[0];
    holder.shell.center[1] = atom.center_bohr[1];
    holder.shell.center[2] = atom.center_bohr[2];
    holder.shell.l = spec.l;
    holder.shell.cart = spec.cart;
    return holder;
}

native::Atom make_native_atom(const AtomSpec& atom) {
    native::Atom out{};
    out.Z = atom.Z;
    out.pos[0] = atom.center_bohr[0];
    out.pos[1] = atom.center_bohr[1];
    out.pos[2] = atom.center_bohr[2];
    return out;
}

const char* one_e_case_label(int case_id) {
    switch (case_id) {
    case 0:
        return "ss_hh";
    case 1:
        return "ps_oh";
    case 2:
        return "pp_oo";
    default:
        return "unknown";
    }
}

const char* eri_case_label(int case_id) {
    switch (case_id) {
    case 0:
        return "ssss";
    case 1:
        return "psss";
    case 2:
        return "ppss";
    case 3:
        return "pppp";
    default:
        return "unknown";
    }
}

struct OneECase {
    const char* label = "";
    std::array<AtomSpec, 2> atoms{};
    std::array<ShellSpec, 2> shells{};
};

OneECase make_one_e_case(int case_id) {
    OneECase out{};
    out.label = one_e_case_label(case_id);

    switch (case_id) {
    case 0:
        out.atoms = {{
            {1, {0.0, 0.0, 0.0}},
            {1, {0.74 * kA2B, 0.0, 0.0}},
        }};
        out.shells = {{
            {0, 1, 0, 0, 0},
            {1, 1, 0, 0, 0},
        }};
        break;
    case 1:
        out.atoms = {{
            {8, {0.0, 0.0, 0.0}},
            {1, {0.96 * kA2B, 0.0, 0.0}},
        }};
        out.shells = {{
            {0, 8, 1, 1, 0},
            {1, 1, 0, 0, 0},
        }};
        break;
    default:
        out.atoms = {{
            {8, {0.0, 0.0, 0.0}},
            {8, {1.20 * kA2B, 0.0, 0.0}},
        }};
        out.shells = {{
            {0, 8, 1, 1, 0},
            {1, 8, 1, 1, 0},
        }};
        break;
    }

    return out;
}

struct ERICase {
    const char* label = "";
    std::array<AtomSpec, 4> atoms{};
    std::array<ShellSpec, 4> shells{};
};

ERICase make_eri_case(int case_id) {
    ERICase out{};
    out.label = eri_case_label(case_id);
    out.atoms = {{
        {8, {0.0 * kA2B, 0.0, 0.0}},
        {8, {1.0 * kA2B, 0.0, 0.0}},
        {1, {2.0 * kA2B, 0.0, 0.0}},
        {1, {3.0 * kA2B, 0.0, 0.0}},
    }};

    switch (case_id) {
    case 0:
        out.shells = {{
            {0, 1, 0, 0, 0},
            {1, 1, 0, 0, 0},
            {2, 1, 0, 0, 0},
            {3, 1, 0, 0, 0},
        }};
        break;
    case 1:
        out.shells = {{
            {0, 8, 1, 1, 0},
            {1, 1, 0, 0, 0},
            {2, 1, 0, 0, 0},
            {3, 1, 0, 0, 0},
        }};
        break;
    case 2:
        out.shells = {{
            {0, 8, 1, 1, 0},
            {1, 8, 1, 1, 1},
            {2, 1, 0, 0, 0},
            {3, 1, 0, 0, 0},
        }};
        break;
    default:
        out.shells = {{
            {0, 8, 1, 1, 0},
            {1, 8, 1, 1, 1},
            {2, 8, 1, 1, 0},
            {3, 8, 1, 1, 1},
        }};
        break;
    }

    return out;
}

struct NativeOneEEnv {
    std::array<ShellHolder, 2> shells{};
    std::array<native::Atom, 2> atoms{};
};

NativeOneEEnv make_native_one_e_env(int case_id) {
    const auto spec = make_one_e_case(case_id);
    NativeOneEEnv env{};
    for (int i = 0; i < 2; ++i) {
        env.shells[i] = make_native_shell(spec.atoms[spec.shells[i].atom_of],
                                          spec.shells[i]);
        env.atoms[i] = make_native_atom(spec.atoms[i]);
    }
    return env;
}

struct NativeERIEnv {
    std::array<ShellHolder, 4> shells{};
};

NativeERIEnv make_native_eri_env(int case_id) {
    const auto spec = make_eri_case(case_id);
    NativeERIEnv env{};
    for (int i = 0; i < 4; ++i)
        env.shells[i] = make_native_shell(spec.atoms[spec.shells[i].atom_of],
                                          spec.shells[i]);
    return env;
}

#ifdef ATV_HAS_LIBCINT
struct LibcintShellInfo {
    int atom_of = 0;
    int l = 0;
    int cart = 0;
    int ncart = 1;
};

struct LibcintEnv {
    std::vector<FINT> atm;
    std::vector<FINT> bas;
    std::vector<double> env;
    std::vector<LibcintShellInfo> shell_info;
    int natm = 0;
    int nbas = 0;
};

LibcintEnv make_libcint_env(const std::vector<AtomSpec>& atoms,
                            const std::vector<ShellSpec>& shells) {
    LibcintEnv out{};
    out.natm = static_cast<int>(atoms.size());
    out.nbas = static_cast<int>(shells.size());
    out.atm.assign(out.natm * ATM_SLOTS, 0);
    out.bas.assign(out.nbas * BAS_SLOTS, 0);
    out.shell_info.resize(out.nbas);

    int off = PTR_ENV_START;
    out.env.assign(off + out.natm * 3 + out.nbas * 6, 0.0);

    for (int i = 0; i < out.natm; ++i) {
        out.atm[CHARGE_OF + ATM_SLOTS * i] = atoms[i].Z;
        out.atm[PTR_COORD + ATM_SLOTS * i] = off;
        out.env[off++] = atoms[i].center_bohr[0];
        out.env[off++] = atoms[i].center_bohr[1];
        out.env[off++] = atoms[i].center_bohr[2];
    }

    for (int i = 0; i < out.nbas; ++i) {
        const auto& spec = shells[i];
        const auto& atom_shells = basis_data::sto3g()[spec.Z].shells;
        const auto& bsh =
            atom_shells[sto3g_basis_data_index(spec.raw_shell, spec.l)];
        const int np = static_cast<int>(bsh.alpha.size());
        const int ptr_exp = off;
        for (int p = 0; p < np; ++p)
            out.env[off++] = bsh.alpha[p];
        const int ptr_coeff = off;
        for (int p = 0; p < np; ++p)
            out.env[off++] = bsh.coeff[p] * CINTgto_norm(spec.l, bsh.alpha[p]);

        out.bas[ATOM_OF + BAS_SLOTS * i] = spec.atom_of;
        out.bas[ANG_OF + BAS_SLOTS * i] = spec.l;
        out.bas[NPRIM_OF + BAS_SLOTS * i] = np;
        out.bas[NCTR_OF + BAS_SLOTS * i] = 1;
        out.bas[PTR_EXP + BAS_SLOTS * i] = ptr_exp;
        out.bas[PTR_COEFF + BAS_SLOTS * i] = ptr_coeff;

        out.shell_info[i].atom_of = spec.atom_of;
        out.shell_info[i].l = spec.l;
        out.shell_info[i].cart = spec.cart;
        out.shell_info[i].ncart = (spec.l == 0) ? 1 : 3;
    }

    out.env.resize(off);
    return out;
}
#endif

struct GtoNormCase {
    const char* label;
    int l;
    double alpha;
};

GtoNormCase make_gto_norm_case(int case_id) {
    static constexpr GtoNormCase kCases[] = {
        {"s_small", 0, 0.16885540},  {"s_medium", 0, 3.42525091},
        {"s_large", 0, 130.7093200}, {"p_small", 1, 0.3803890},
        {"p_medium", 1, 1.1695961},  {"p_large", 1, 5.0331513},
    };
    return kCases[case_id];
}

} // namespace

#ifdef ATV_HAS_LIBCINT
static void BM_GtoNorm_Libcint(benchmark::State& state) {
    const auto spec = make_gto_norm_case(static_cast<int>(state.range(0)));
    state.SetLabel(spec.label);
    for (auto _ : state)
        benchmark::DoNotOptimize(CINTgto_norm(spec.l, spec.alpha));
}
BENCHMARK(BM_GtoNorm_Libcint)->DenseRange(0, 5);
#endif

static void BM_GtoNorm_Native(benchmark::State& state) {
    const auto spec = make_gto_norm_case(static_cast<int>(state.range(0)));
    state.SetLabel(spec.label);
    for (auto _ : state)
        benchmark::DoNotOptimize(native::gto_norm(spec.l, spec.alpha));
}
BENCHMARK(BM_GtoNorm_Native)->DenseRange(0, 5);

#ifdef ATV_HAS_LIBCINT
static void BM_Overlap_Libcint(benchmark::State& state) {
    const int case_id = static_cast<int>(state.range(0));
    const auto spec = make_one_e_case(case_id);
    state.SetLabel(spec.label);
    auto env = make_libcint_env({spec.atoms.begin(), spec.atoms.end()},
                                {spec.shells.begin(), spec.shells.end()});

    std::array<double, 9> buf{};
    FINT shls[2] = {0, 1};
    const int ni = env.shell_info[0].ncart;
    const int out_idx = env.shell_info[0].cart + ni * env.shell_info[1].cart;

    for (auto _ : state) {
        int1e_ovlp_cart(buf.data(), nullptr, shls, env.atm.data(), env.natm,
                        env.bas.data(), env.nbas, env.env.data(), nullptr,
                        nullptr);
        benchmark::DoNotOptimize(buf[out_idx]);
    }
}
BENCHMARK(BM_Overlap_Libcint)->DenseRange(0, 2);
#endif

static void BM_Overlap_Native(benchmark::State& state) {
    const int case_id = static_cast<int>(state.range(0));
    state.SetLabel(one_e_case_label(case_id));
    const auto env = make_native_one_e_env(case_id);
    for (auto _ : state)
        benchmark::DoNotOptimize(
            native::overlap(env.shells[0].shell, env.shells[1].shell));
}
BENCHMARK(BM_Overlap_Native)->DenseRange(0, 2);

#ifdef ATV_HAS_LIBCINT
static void BM_Kinetic_Libcint(benchmark::State& state) {
    const int case_id = static_cast<int>(state.range(0));
    const auto spec = make_one_e_case(case_id);
    state.SetLabel(spec.label);
    auto env = make_libcint_env({spec.atoms.begin(), spec.atoms.end()},
                                {spec.shells.begin(), spec.shells.end()});

    std::array<double, 9> buf{};
    FINT shls[2] = {0, 1};
    const int ni = env.shell_info[0].ncart;
    const int out_idx = env.shell_info[0].cart + ni * env.shell_info[1].cart;

    for (auto _ : state) {
        int1e_kin_cart(buf.data(), nullptr, shls, env.atm.data(), env.natm,
                       env.bas.data(), env.nbas, env.env.data(), nullptr,
                       nullptr);
        benchmark::DoNotOptimize(buf[out_idx]);
    }
}
BENCHMARK(BM_Kinetic_Libcint)->DenseRange(0, 2);
#endif

static void BM_Kinetic_Native(benchmark::State& state) {
    const int case_id = static_cast<int>(state.range(0));
    state.SetLabel(one_e_case_label(case_id));
    const auto env = make_native_one_e_env(case_id);
    for (auto _ : state)
        benchmark::DoNotOptimize(
            native::kinetic(env.shells[0].shell, env.shells[1].shell));
}
BENCHMARK(BM_Kinetic_Native)->DenseRange(0, 2);

#ifdef ATV_HAS_LIBCINT
static void BM_Nuclear_Libcint(benchmark::State& state) {
    const int case_id = static_cast<int>(state.range(0));
    const auto spec = make_one_e_case(case_id);
    state.SetLabel(spec.label);
    auto env = make_libcint_env({spec.atoms.begin(), spec.atoms.end()},
                                {spec.shells.begin(), spec.shells.end()});

    std::array<double, 9> buf{};
    FINT shls[2] = {0, 1};
    const int ni = env.shell_info[0].ncart;
    const int out_idx = env.shell_info[0].cart + ni * env.shell_info[1].cart;

    for (auto _ : state) {
        int1e_nuc_cart(buf.data(), nullptr, shls, env.atm.data(), env.natm,
                       env.bas.data(), env.nbas, env.env.data(), nullptr,
                       nullptr);
        benchmark::DoNotOptimize(buf[out_idx]);
    }
}
BENCHMARK(BM_Nuclear_Libcint)->DenseRange(0, 2);
#endif

static void BM_Nuclear_Native(benchmark::State& state) {
    const int case_id = static_cast<int>(state.range(0));
    state.SetLabel(one_e_case_label(case_id));
    const auto env = make_native_one_e_env(case_id);
    for (auto _ : state)
        benchmark::DoNotOptimize(native::nuclear(
            env.shells[0].shell, env.shells[1].shell, env.atoms.data(), 2));
}
BENCHMARK(BM_Nuclear_Native)->DenseRange(0, 2);

#ifdef ATV_HAS_LIBCINT
static void BM_ERI_Libcint(benchmark::State& state) {
    const int case_id = static_cast<int>(state.range(0));
    const auto spec = make_eri_case(case_id);
    state.SetLabel(spec.label);
    auto env = make_libcint_env({spec.atoms.begin(), spec.atoms.end()},
                                {spec.shells.begin(), spec.shells.end()});

    std::array<double, 81> buf{};
    FINT shls[4] = {0, 1, 2, 3};
    const int n0 = env.shell_info[0].ncart;
    const int n1 = env.shell_info[1].ncart;
    const int n2 = env.shell_info[2].ncart;
    const int out_idx =
        env.shell_info[0].cart
        + n0
              * (env.shell_info[1].cart
                 + n1 * (env.shell_info[2].cart + n2 * env.shell_info[3].cart));

    for (auto _ : state) {
        cint2e_cart(buf.data(), shls, env.atm.data(), env.natm, env.bas.data(),
                    env.nbas, env.env.data(), nullptr);
        benchmark::DoNotOptimize(buf[out_idx]);
    }
}
BENCHMARK(BM_ERI_Libcint)->DenseRange(0, 3);
#endif

static void BM_ERI_Native(benchmark::State& state) {
    const int case_id = static_cast<int>(state.range(0));
    state.SetLabel(eri_case_label(case_id));
    const auto env = make_native_eri_env(case_id);
    for (auto _ : state)
        benchmark::DoNotOptimize(
            native::eri(env.shells[0].shell, env.shells[1].shell,
                        env.shells[2].shell, env.shells[3].shell));
}
BENCHMARK(BM_ERI_Native)->DenseRange(0, 3);

// ═════════════════════════════════════════════════════════════════════════
// Shell-block ERI: eri_block() vs cint2e_cart() (apples-to-apples)
//
// Both compute ALL Cartesian ERIs for a shell quartet in one call.
// ═════════════════════════════════════════════════════════════════════════

namespace {

struct ShellBlockHolder {
    std::array<double, 8> alpha{};
    std::array<double, 8> coeff{};
    native::ShellBlock block{};
};

ShellBlockHolder make_native_shell_block(const AtomSpec& atom,
                                         const ShellSpec& spec) {
    const auto& atom_shells = basis_data::sto3g()[spec.Z].shells;
    const auto& bsh =
        atom_shells[sto3g_basis_data_index(spec.raw_shell, spec.l)];
    const int np = static_cast<int>(bsh.alpha.size());

    ShellBlockHolder holder;
    double (*norm_at)(double) =
        (spec.l == 0) ? native::phys_norm : native::phys_norm_p;
    for (int p = 0; p < np; ++p) {
        holder.alpha[p] = bsh.alpha[p];
        holder.coeff[p] = bsh.coeff[p] * norm_at(bsh.alpha[p]);
    }

    holder.block.n_prim = np;
    holder.block.alpha = holder.alpha.data();
    holder.block.coeff = holder.coeff.data();
    holder.block.center[0] = atom.center_bohr[0];
    holder.block.center[1] = atom.center_bohr[1];
    holder.block.center[2] = atom.center_bohr[2];
    holder.block.l = spec.l;
    holder.block.n_cart = (spec.l == 0) ? 1 : 3;
    holder.block.ao_start = 0;
    return holder;
}

struct NativeERIBlockEnv {
    std::array<ShellBlockHolder, 4> blocks{};
};

NativeERIBlockEnv make_native_eri_block_env(int case_id) {
    const auto spec = make_eri_case(case_id);
    NativeERIBlockEnv env{};
    for (int i = 0; i < 4; ++i)
        env.blocks[i] = make_native_shell_block(
            spec.atoms[spec.shells[i].atom_of], spec.shells[i]);
    return env;
}

} // namespace

/// @brief Native eri_block() — computes ALL ERIs for a shell quartet at once.
static void BM_ERIBlock_Native(benchmark::State& state) {
    const int case_id = static_cast<int>(state.range(0));
    state.SetLabel(eri_case_label(case_id));
    const auto env = make_native_eri_block_env(case_id);
    double buf[81] = {};
    for (auto _ : state) {
        std::fill(buf, buf + 81, 0.0);
        native::eri_block(env.blocks[0].block, env.blocks[1].block,
                          env.blocks[2].block, env.blocks[3].block, buf);
        benchmark::DoNotOptimize(buf[0]);
    }
}
BENCHMARK(BM_ERIBlock_Native)->DenseRange(0, 3);

#ifdef ATV_HAS_LIBCINT
/// @brief libcint cint2e_cart() — computes ALL ERIs for a shell quartet at once.
/// Direct comparison with BM_ERIBlock_Native above.
static void BM_ERIBlock_Libcint(benchmark::State& state) {
    const int case_id = static_cast<int>(state.range(0));
    const auto spec = make_eri_case(case_id);
    state.SetLabel(spec.label);
    auto env = make_libcint_env({spec.atoms.begin(), spec.atoms.end()},
                                {spec.shells.begin(), spec.shells.end()});

    double buf[81] = {};
    FINT shls[4] = {0, 1, 2, 3};

    for (auto _ : state) {
        cint2e_cart(buf, shls, env.atm.data(), env.natm, env.bas.data(),
                    env.nbas, env.env.data(), nullptr);
        benchmark::DoNotOptimize(buf[0]);
    }
}
BENCHMARK(BM_ERIBlock_Libcint)->DenseRange(0, 3);
#endif

int main(int argc, char** argv) {
    ATV::log_init(spdlog::level::warn);
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv))
        return 1;
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}

#else
int main() { return 0; }
#endif
