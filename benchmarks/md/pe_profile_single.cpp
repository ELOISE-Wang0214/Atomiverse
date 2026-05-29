/// @file benchmarks/md/pe_profile_single.cpp
/// @brief Phase 5 Track B.1 — per-component scoped-timer profile of the
///        amorphous-PE tension hot path.
///
/// Builds the same component pipeline as `polyethylene_tension.cpp`
/// but uses the base `ATV::Driver<md::Phase, cpu::device_tag>` so we
/// can wrap every hot-path component in a `TimedComponent` decorator
/// (instead of `cpu::md::MD`, whose constructor auto-adds untimed
/// infrastructure components).
///
/// Wrapped components (per acceptance criterion ac-008):
///   HaloExchangeForward, HaloExchangeReverse, NeighborList, LJ,
///   HarmonicBond, HarmonicAngle, MultiHarmonicDihedral, VelocityVerlet,
///   NoseHoover, ReduceEnergy, SnapshotPublisher.
///
/// Plus the small zero/init infrastructure components (ZeroForces,
/// ZeroEnergy, ZeroStress, ZeroKinetic, KineticEnergy, AtomMigration,
/// InitVelocity) which are also instrumented.
///
/// CSV output schema (`--out`):
///   `component,total_ms,pct_of_step,n_calls,mean_us`
///
/// Acceptance: ac-008 / ac-009 / ac-015 in
/// `.claude/specs/parallel-md-engine.acceptance.md`.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/core/device.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/backend_traits.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/comm/atom_migration.h"
#include "atomiverse/cpu/comm/halo_exchange.h"
#include "atomiverse/cpu/md/mpi/reduce.h"
#include "atomiverse/cpu/md/neighbor_list.h"
#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/potentials/harmonic_angle.h"
#include "atomiverse/cpu/potentials/harmonic_bond.h"
#include "atomiverse/cpu/potentials/multi_harmonic_dihedral.h"
#include "atomiverse/cpu/snapshot_publisher.h"
#include "atomiverse/facade/driver.h"
#include "atomiverse/facade/stage.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

using namespace ATV;
using ATV::md::Phase;

// ───────────────────────────── Profile registry ─────────────────────
namespace {

class ProfileRegistry {
  public:
    static ProfileRegistry& instance() {
        static ProfileRegistry r;
        return r;
    }

    struct Entry {
        std::atomic<uint64_t> total_ns{0};
        std::atomic<uint64_t> n_calls{0};
    };

    void record(const std::string& name, uint64_t ns) {
        auto& e = entries_[name];
        e.total_ns.fetch_add(ns, std::memory_order_relaxed);
        e.n_calls.fetch_add(1, std::memory_order_relaxed);
    }

    void dump_csv(std::ostream& os, double step_total_ms) const {
        os << "component,total_ms,pct_of_step,n_calls,mean_us\n";
        for (const auto& kv : entries_) {
            const double total_ms =
                kv.second.total_ns.load(std::memory_order_relaxed) / 1e6;
            const uint64_t n =
                kv.second.n_calls.load(std::memory_order_relaxed);
            const double pct =
                step_total_ms > 0.0 ? 100.0 * total_ms / step_total_ms : 0.0;
            const double mean_us = n > 0 ? (total_ms * 1000.0 / n) : 0.0;
            os << kv.first << "," << total_ms << "," << pct << "," << n << ","
               << mean_us << "\n";
        }
    }

    void reset() {
        for (auto& kv : entries_) {
            kv.second.total_ns.store(0, std::memory_order_relaxed);
            kv.second.n_calls.store(0, std::memory_order_relaxed);
        }
    }

  private:
    ProfileRegistry() = default;
    std::map<std::string, Entry> entries_;
};

/// @brief Wraps a `cpu::ComponentBase` with per-call wall-time recording.
///
/// Carries an explicit phase/order so it can be added via `Driver::on()`
/// without depending on the wrapped type's `default_bindings()`.
class TimedComponent : public cpu::ComponentBase {
  public:
    TimedComponent(std::unique_ptr<cpu::ComponentBase> inner, std::string label,
                   int phase, int order)
        : inner_(std::move(inner)), label_(std::move(label)), phase_(phase),
          order_(order) {}

    std::string name() const override { return label_; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{phase_, order_}};
    }

    void build(const cpu::BuildContext& ctx) override { inner_->build(ctx); }

    void launch(const cpu::LaunchContext& ctx) override {
        const auto t0 = std::chrono::steady_clock::now();
        inner_->launch(ctx);
        const auto t1 = std::chrono::steady_clock::now();
        const auto ns =
            std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)
                .count();
        ProfileRegistry::instance().record(label_, static_cast<uint64_t>(ns));
    }

    int frequency() const override { return inner_->frequency(); }

  private:
    std::unique_ptr<cpu::ComponentBase> inner_;
    std::string label_;
    int phase_;
    int order_;
};

template <typename C>
TimedComponent* make_timed(C&& comp, std::string label, int phase, int order) {
    auto inner =
        std::make_unique<std::remove_reference_t<C>>(std::forward<C>(comp));
    // Heap-allocated; lifetime is process-end (intentional leak — the
    // benchmark is one-shot and the driver requires lvalue lifetime).
    return new TimedComponent(std::move(inner), std::move(label), phase, order);
}

// ───────────────────────────── Force-field constants ────────────────
constexpr double kKcalPerMolToEv = 1.0 / 23.060547830619026;
constexpr double kBondK_kcal = 350.0;
constexpr double kBondR0 = 1.53;
constexpr double kAngleK_kcal = 60.0;
constexpr double kAngleTheta0 = 109.5 * kPI / 180.0;
constexpr double kLjEps_kcal = 0.112;
constexpr double kLjSigma = 4.01;
constexpr double kLjCutoff = 10.5;

double bond_k_eV() { return kBondK_kcal * kKcalPerMolToEv; }
double angle_k_eV() { return kAngleK_kcal * kKcalPerMolToEv; }
double lj_eps_eV() { return kLjEps_kcal * kKcalPerMolToEv; }

std::vector<double> dihedral_Cn_eV() {
    return {
        1.73 * kKcalPerMolToEv,  // C0  (MSU paper, Tschopp 2010)
        -4.49 * kKcalPerMolToEv, // C1
        0.776 * kKcalPerMolToEv, // C2
        6.99 * kKcalPerMolToEv,  // C3
        0.0,                     // C4
    };
}

cpu::Frame make_pe_box(int chains_per_axis, int n_carbons, double box_len,
                       double margin) {
    const int n_chains = chains_per_axis * chains_per_axis * chains_per_axis;
    const int n_atoms = n_chains * n_carbons;
    cpu::Frame frame(n_atoms);

    const double half_apex = (kPI - kAngleTheta0) / 2.0;
    const double dy_step = kBondR0 * std::cos(half_apex);
    const double dz_step = kBondR0 * std::sin(half_apex);

    const double cell = box_len / chains_per_axis;
    int idx = 0;
    for (int cx = 0; cx < chains_per_axis; ++cx) {
        for (int cy = 0; cy < chains_per_axis; ++cy) {
            for (int cz = 0; cz < chains_per_axis; ++cz) {
                const double cx0 = cx * cell + cell / 2.0;
                const double cz0 = cz * cell + cell / 2.0;
                const double y0 = cy * cell + margin;
                for (int i = 0; i < n_carbons; ++i) {
                    frame.atoms.x[idx] = cx0;
                    frame.atoms.y[idx] = y0 + i * dy_step;
                    frame.atoms.z[idx] = cz0 + (i % 2 == 0 ? 0.0 : dz_step);
                    frame.atoms.mass[idx] = 14.027;
                    frame.atoms.type[idx] = 0;
                    ++idx;
                }
            }
        }
    }

    const int nb_per = n_carbons - 1;
    const int na_per = n_carbons - 2;
    const int nd_per = n_carbons - 3;
    frame.bonds.atomi = xt::zeros<int>({n_chains * nb_per});
    frame.bonds.atomj = xt::zeros<int>({n_chains * nb_per});
    frame.bonds.type = xt::zeros<int>({n_chains * nb_per});
    frame.angles.atomi = xt::zeros<int>({n_chains * na_per});
    frame.angles.atomj = xt::zeros<int>({n_chains * na_per});
    frame.angles.atomk = xt::zeros<int>({n_chains * na_per});
    frame.angles.type = xt::zeros<int>({n_chains * na_per});
    frame.dihedrals.atomi = xt::zeros<int>({n_chains * nd_per});
    frame.dihedrals.atomj = xt::zeros<int>({n_chains * nd_per});
    frame.dihedrals.atomk = xt::zeros<int>({n_chains * nd_per});
    frame.dihedrals.atoml = xt::zeros<int>({n_chains * nd_per});
    frame.dihedrals.type = xt::zeros<int>({n_chains * nd_per});

    for (int c = 0; c < n_chains; ++c) {
        const int base = c * n_carbons;
        for (int k = 0; k < nb_per; ++k) {
            const int b = c * nb_per + k;
            frame.bonds.atomi[b] = base + k;
            frame.bonds.atomj[b] = base + k + 1;
        }
        for (int k = 0; k < na_per; ++k) {
            const int aa = c * na_per + k;
            frame.angles.atomi[aa] = base + k;
            frame.angles.atomj[aa] = base + k + 1;
            frame.angles.atomk[aa] = base + k + 2;
        }
        for (int k = 0; k < nd_per; ++k) {
            const int d = c * nd_per + k;
            frame.dihedrals.atomi[d] = base + k;
            frame.dihedrals.atomj[d] = base + k + 1;
            frame.dihedrals.atomk[d] = base + k + 2;
            frame.dihedrals.atoml[d] = base + k + 3;
        }
    }
    frame.box = Mat3::cubic(box_len);
    return frame;
}

struct Args {
    int steps = 1000;
    int warmup = 100;
    int chains_per_axis = 4;
    int n_atoms_per_chain = 18;
    double box_len = 60.0;
    double T_target = 300.0;
    double dt = 0.5;
    double tau_T = 100.0;
    int seed = 12345;
    std::string out_csv;
    bool quiet = false;
};

void print_help() {
    std::printf(
        "Usage: pe_profile_single [options]\n"
        "  --steps N             Steps to time (default 1000).\n"
        "  --warmup N            Warm-up steps (default 100).\n"
        "  --chains-per-axis N   NxNxN chain grid (default 4).\n"
        "  --n-atoms-per-chain N UA atoms per chain (default 18).\n"
        "  --box-len L           Box edge in Å (default 60).\n"
        "  --T K                 Setpoint (default 300).\n"
        "  --dt fs               Time step (default 0.5).\n"
        "  --tau-T fs            Thermostat τ (default 100).\n"
        "  --seed S              RNG seed (default 12345).\n"
        "  --out PATH            CSV path (default profile_<ranks>.csv).\n"
        "  --quiet               Suppress progress output.\n"
        "  --help                Show this message.\n");
}

bool parse_args(int argc, char** argv, Args& a) {
    for (int i = 1; i < argc; ++i) {
        const char* s = argv[i];
        auto take_d = [&](const char* n, double& d) {
            if (std::strcmp(s, n) != 0)
                return false;
            if (i + 1 >= argc)
                std::exit(2);
            d = std::atof(argv[++i]);
            return true;
        };
        auto take_i = [&](const char* n, int& d) {
            if (std::strcmp(s, n) != 0)
                return false;
            if (i + 1 >= argc)
                std::exit(2);
            d = std::atoi(argv[++i]);
            return true;
        };
        auto take_s = [&](const char* n, std::string& d) {
            if (std::strcmp(s, n) != 0)
                return false;
            if (i + 1 >= argc)
                std::exit(2);
            d = argv[++i];
            return true;
        };
        if (take_i("--steps", a.steps))
            continue;
        if (take_i("--warmup", a.warmup))
            continue;
        if (take_i("--chains-per-axis", a.chains_per_axis))
            continue;
        if (take_i("--n-atoms-per-chain", a.n_atoms_per_chain))
            continue;
        if (take_d("--box-len", a.box_len))
            continue;
        if (take_d("--T", a.T_target))
            continue;
        if (take_d("--dt", a.dt))
            continue;
        if (take_d("--tau-T", a.tau_T))
            continue;
        if (take_i("--seed", a.seed))
            continue;
        if (take_s("--out", a.out_csv))
            continue;
        if (std::strcmp(s, "--quiet") == 0) {
            a.quiet = true;
            continue;
        }
        if (std::strcmp(s, "--help") == 0) {
            print_help();
            std::exit(0);
        }
        std::fprintf(stderr, "unknown argument: %s\n", s);
        return false;
    }
    return true;
}

double now_seconds() {
#ifdef ATV_HAS_MPI
    return MPI_Wtime();
#else
    using clk = std::chrono::steady_clock;
    static const auto t0 = clk::now();
    return std::chrono::duration<double>(clk::now() - t0).count();
#endif
}

} // namespace

int main(int argc, char** argv) {
    Args a;
    if (!parse_args(argc, argv, a))
        return 2;

    int my_rank = 0;
    int world_size = 1;
#ifdef ATV_HAS_MPI
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

    if (a.out_csv.empty()) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "profile_%d.csv", world_size);
        a.out_csv = buf;
    }

    auto frame = make_pe_box(a.chains_per_axis, a.n_atoms_per_chain, a.box_len,
                             /*margin=*/1.0);
    const int n_atoms = frame.atoms.size();

    // Use the base Driver template directly so we control every binding.
    Driver<Phase, cpu::device_tag> sim;

    const int p_init = static_cast<int>(Phase::initial_integrate);
    const int p_pre = static_cast<int>(Phase::pre_force);
    const int p_force = static_cast<int>(Phase::force);
    const int p_post = static_cast<int>(Phase::post_force);
    const int p_final = static_cast<int>(Phase::final_integrate);
    const int p_eos = static_cast<int>(Phase::end_of_step);

    // Init helper (one-shot at step 0): bound at initial_integrate@10
    // through its default_bindings; we instrument it the same way.
    sim.on(Phase::initial_integrate)
        .add(*make_timed(
                 cpu::md::InitVelocity{a.T_target,
                                       static_cast<unsigned int>(a.seed)},
                 "InitVelocity", p_init, 10),
             10);
    sim.on(Phase::initial_integrate)
        .add(*make_timed(cpu::md::VelocityVerlet{a.dt}, "VelocityVerlet",
                         p_init, 0),
             0);
    sim.on(Phase::initial_integrate)
        .add(*make_timed(cpu::md::NoseHoover{a.T_target, a.tau_T, a.dt},
                         "NoseHoover", p_init, 100),
             100);

    sim.on(Phase::pre_force)
        .add(*make_timed(cpu::comm::AtomMigration{10}, "AtomMigration",
                         p_pre, -30),
             -30);
    sim.on(Phase::pre_force)
        .add(*make_timed(cpu::comm::HaloExchangeForward{kLjCutoff},
                         "HaloExchangeForward", p_pre, -20),
             -20);

    sim.on(Phase::force)
        .add(*make_timed(cpu::md::ZeroForces{}, "ZeroForces", p_force, -100),
             -100);
    sim.on(Phase::force)
        .add(*make_timed(cpu::md::ZeroEnergy{}, "ZeroEnergy", p_force, -99),
             -99);
    sim.on(Phase::force)
        .add(*make_timed(cpu::md::ZeroStress{}, "ZeroStress", p_force, -98),
             -98);
    sim.on(Phase::force)
        .add(*make_timed(cpu::md::NeighborList<>{kLjCutoff}, "NeighborList",
                         p_force, -10),
             -10);
    sim.on(Phase::force)
        .add(*make_timed(cpu::md::LJ{lj_eps_eV(), kLjSigma, kLjCutoff}
                             .set_special_bonds(0.0, 0.0, 0.0),
                         "LJ", p_force, 0),
             0);
    sim.on(Phase::force)
        .add(*make_timed(HarmonicBond{{bond_k_eV()}, {kBondR0}}, "HarmonicBond",
                         p_force, 1),
             1);
    sim.on(Phase::force)
        .add(*make_timed(HarmonicAngle{{angle_k_eV()}, {kAngleTheta0}},
                         "HarmonicAngle", p_force, 2),
             2);
    sim.on(Phase::force)
        .add(*make_timed(MultiHarmonicDihedral{{dihedral_Cn_eV()}},
                         "MultiHarmonicDihedral", p_force, 3),
             3);

    sim.on(Phase::post_force)
        .add(*make_timed(cpu::comm::HaloExchangeReverse{},
                         "HaloExchangeReverse", p_post, -1000),
             -1000);

    sim.on(Phase::final_integrate)
        .add(*make_timed(cpu::md::VelocityVerlet{a.dt}, "VelocityVerlet_final",
                         p_final, 0),
             0);

    sim.on(Phase::end_of_step)
        .add(*make_timed(cpu::md::ZeroKinetic{}, "ZeroKinetic", p_eos, -1000),
             -1000);
    sim.on(Phase::end_of_step)
        .add(
            *make_timed(cpu::md::KineticEnergy{}, "KineticEnergy", p_eos, -100),
            -100);
    sim.on(Phase::end_of_step)
        .add(*make_timed(cpu::md::mpi::ReduceEnergy{}, "ReduceEnergy", p_eos,
                         -50),
             -50);
    sim.on(Phase::end_of_step)
        .add(*make_timed(cpu::SnapshotPublisher{}, "SnapshotPublisher", p_eos,
                         10),
             10);

    sim.build();

    if (my_rank == 0 && !a.quiet) {
        std::fprintf(
            stderr, "pe_profile_single: ranks=%d atoms=%d warmup=%d steps=%d\n",
            world_size, n_atoms, a.warmup, a.steps);
    }

    // Warm-up: run but discard timings.
    if (a.warmup > 0)
        sim.run(frame, a.warmup);
    ProfileRegistry::instance().reset();

    const double t0 = now_seconds();
    sim.run(frame, a.steps);
    const double t1 = now_seconds();
    const double total_step_ms = 1000.0 * (t1 - t0);

    if (my_rank == 0) {
        std::ofstream csv(a.out_csv);
        if (!csv.is_open()) {
            std::fprintf(stderr, "error: cannot open %s for write\n",
                         a.out_csv.c_str());
            return 3;
        }
        ProfileRegistry::instance().dump_csv(csv, total_step_ms);
        csv.close();
        if (!a.quiet) {
            std::fprintf(stderr,
                         "  wrote %s — total_step_wall=%.2fms over %d steps "
                         "(%.4f ms/step)\n",
                         a.out_csv.c_str(), total_step_ms, a.steps,
                         total_step_ms / std::max(1, a.steps));
        }
    }

#ifdef ATV_HAS_MPI
    MPI_Finalize();
#endif
    return 0;
}
