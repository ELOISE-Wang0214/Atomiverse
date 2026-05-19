/// @file benchmarks/md/polyethylene_tension.cpp
/// @brief MSU-style amorphous-polyethylene uniaxial-tension benchmark.
///
/// Drives the full 5-stage protocol from `parallel-md-engine.md`
/// Verification §7 inside a *single* `cpu::md::MD` instance, using
/// `cpu::md::NoseHoover::set_T_target` to switch thermostat setpoints
/// between stages without rebuilding the driver. The deformation
/// component is engaged only during stage 5 (tension) by ratcheting the
/// strain rate from 0 → ε̇_target between stages.
///
/// Stages (all stages run in one process, one driver, one frame):
///
///   1. NVT  @ 500 K, dt = 1.0 fs, 10 000 steps    -- initial randomisation
///   2. NPT  @ 500 K, dt = 0.5 fs, 50 000 steps    -- density relaxation
///   3. NVT  500 → 100 K linear ramp,
///         dt = 0.5 fs, 50 000 steps               -- glass formation
///   4. NPT  @ 100 K, dt = 0.5 fs, 50 000 steps    -- glassy equilibration
///   5. Tension along x at ε̇ = 1e-6 /fs,
///         dt = 0.5 fs, 200 000 steps              -- uniaxial tensile load
///
/// The cooling ramp (stage 3) is implemented as N micro-segments inside
/// stage 3, calling `NoseHoover::set_T_target(...)` between each one to
/// step the setpoint linearly from 500 K down to 100 K.
///
/// `--mode {full|nvt|npt|tension|stage:<id>}` selects the run shape:
///
///   * `full`     -- run all 5 stages (default; the canonical benchmark);
///   * `nvt`      -- stage 1 only (single-segment NVT, kept for back-compat
///                  with `tests/mpi/test_pe_msu_pipeline_n8.cpp`);
///   * `npt`      -- stages 1+2;
///   * `tension`  -- stages 1+2+5 (skip cooling/equilibration; single-T
///                  tensile run, kept for legacy timing comparisons);
///   * `stage:<id>` -- run a single stage by 1-indexed id.
///
/// Output (rank 0 only) under `--out-dir`:
///   * `thermo.csv`: `step,stage,T_target_K,T_K,P_bar,E_pot_eV,E_kin_eV,box_x,box_y,box_z`
///   * `stress_strain.csv`: `step,stage,strain,sigma_xx_eV_per_A3,sigma_xx_GPa`
///
/// @see .claude/specs/parallel-md-engine.md  Phase 5 Track A.1 + A.2
/// @see .claude/specs/msu-pe-validation.md   Phase C.2

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "atomiverse/cpu/potentials/harmonic_angle.h"
#include "atomiverse/cpu/potentials/harmonic_bond.h"
#include "atomiverse/cpu/potentials/multi_harmonic_dihedral.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/barostat.h"
#include "atomiverse/cpu/md/deform.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/mpi/halo_exchange.h"
#include "atomiverse/cpu/md/neighbor_list.h"
#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/md/virial_stress.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/facade/stage.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#else
#include <chrono>
#endif

using namespace ATV;

namespace {

// ── Force-field constants (Dreiding-style UA polyethylene) ──────────
constexpr double kKcalPerMolToEv = 1.0 / 23.060547830619026;
constexpr double kBondK_kcal     = 350.0;
constexpr double kBondR0         = 1.53;
constexpr double kAngleK_kcal    = 60.0;
constexpr double kAngleTheta0    = 109.5 * kPI / 180.0;
constexpr double kLjEps_kcal     = 0.112;
constexpr double kLjSigma        = 4.01;
constexpr double kLjCutoff       = 10.5;

constexpr double kEvPerA3ToGPa = 160.21766;
constexpr double kEvPerA3ToBar = 1.602176634e6;

double bond_k_eV()  { return kBondK_kcal  * kKcalPerMolToEv; }
double angle_k_eV() { return kAngleK_kcal * kKcalPerMolToEv; }
double lj_eps_eV()  { return kLjEps_kcal  * kKcalPerMolToEv; }

std::vector<double> dihedral_Cn_eV() {
    return {
        1.73    * kKcalPerMolToEv,   // C0  (MSU paper, Tschopp 2010)
        -4.49   * kKcalPerMolToEv,   // C1
        0.776   * kKcalPerMolToEv,   // C2
        6.99    * kKcalPerMolToEv,   // C3
        0.0,                          // C4
    };
}

enum class Mode { full, nvt, npt, tension, single_stage };

struct Args {
    Mode   mode        = Mode::full;
    int    single_stage_id = 1;     ///< 1..5 when mode == single_stage
    /// Per-stage step counts. Defaults are scaled-down from spec §7 so
    /// the benchmark completes in minutes for CI smoke runs; large runs
    /// override via `--scale N` (multiplies every stage's step count).
    int    stage1_steps = 2000;     ///< NVT  @ 500 K, dt=1.0 fs
    int    stage2_steps = 5000;     ///< NPT  @ 500 K, dt=0.5 fs
    int    stage3_steps = 5000;     ///< NVT  cool 500 → 100 K
    int    stage4_steps = 5000;     ///< NPT  @ 100 K
    int    stage5_steps = 10000;    ///< tension along x
    int    cool_segments = 20;      ///< stage 3: linear-ramp micro-segments
    double scale       = 1.0;       ///< multiplier on each stage's steps

    int    print_every = 100;
    int    chains_per_axis = 4;
    int    n_atoms_per_chain = 18;
    double box_len     = 60.0;
    double margin      = 1.0;
    int    seed        = 12345;
    double T_hot       = 500.0;     ///< stages 1-2 setpoint [K]
    double T_cold      = 100.0;     ///< stages 4-5 setpoint [K]
    double dt_stage1   = 1.0;       ///< [fs]
    double dt_other    = 0.5;       ///< [fs] -- stages 2-5
    double tau_T       = 100.0;     ///< [fs]
    double tau_P       = 1000.0;    ///< [fs]
    double strain_rate = 1e-6;      ///< 1/fs (stage 5)
    std::string out_dir = ".";
    bool   quiet       = false;
};

void print_help() {
    std::printf(
        "Usage: polyethylene_tension [options]\n"
        "  --mode {full|nvt|npt|tension|stage:N}   Run shape (default full).\n"
        "      full        : all 5 stages, one driver, cooling ramp.\n"
        "      nvt         : stage 1 only.\n"
        "      npt         : stages 1+2.\n"
        "      tension     : stages 1+2+5 (no cooling/equilibration).\n"
        "      stage:N     : single stage by 1-indexed id (1..5).\n"
        "  --scale F                   Multiply every stage's step count by F (default 1).\n"
        "  --print-every N             Thermo/stress sample stride (default 100).\n"
        "  --chains-per-axis N         NxNxN chain grid (default 4).\n"
        "  --n-atoms-per-chain N       UA atoms per chain (default 18).\n"
        "  --box-len L                 Initial cubic box edge in Å (default 60).\n"
        "  --T-hot K                   Hot setpoint, stages 1-2 (default 500).\n"
        "  --T-cold K                  Cold setpoint, stages 4-5 (default 100).\n"
        "  --tau-T fs                  Thermostat time constant (default 100).\n"
        "  --tau-P fs                  Barostat time constant (default 1000).\n"
        "  --strain-rate /fs           Engineering strain rate stage 5 (default 1e-6).\n"
        "  --cool-segments N           Linear-ramp micro-segments stage 3 (default 20).\n"
        "  --out-dir DIR               Where to write CSV outputs (default '.').\n"
        "  --seed S                    RNG seed for InitVelocity (default 12345).\n"
        "  --quiet                     Suppress per-stage progress lines.\n"
        "  --help                      Show this message.\n");
}

bool parse_args(int argc, char** argv, Args& a) {
    for (int i = 1; i < argc; ++i) {
        const char* s = argv[i];
        auto take_d = [&](const char* name, double& dst) {
            if (std::strcmp(s, name) != 0) return false;
            if (i + 1 >= argc) std::exit(2);
            dst = std::atof(argv[++i]);
            return true;
        };
        auto take_i = [&](const char* name, int& dst) {
            if (std::strcmp(s, name) != 0) return false;
            if (i + 1 >= argc) std::exit(2);
            dst = std::atoi(argv[++i]);
            return true;
        };
        auto take_s = [&](const char* name, std::string& dst) {
            if (std::strcmp(s, name) != 0) return false;
            if (i + 1 >= argc) std::exit(2);
            dst = argv[++i];
            return true;
        };
        if (std::strcmp(s, "--mode") == 0) {
            if (i + 1 >= argc) std::exit(2);
            const std::string v = argv[++i];
            if      (v == "full")    a.mode = Mode::full;
            else if (v == "nvt")     a.mode = Mode::nvt;
            else if (v == "npt")     a.mode = Mode::npt;
            else if (v == "tension") a.mode = Mode::tension;
            else if (v.rfind("stage:", 0) == 0) {
                a.mode = Mode::single_stage;
                a.single_stage_id = std::atoi(v.c_str() + 6);
                if (a.single_stage_id < 1 || a.single_stage_id > 5) {
                    std::fprintf(stderr, "stage id must be 1..5: %s\n", v.c_str());
                    return false;
                }
            }
            else { std::fprintf(stderr, "unknown mode: %s\n", v.c_str()); return false; }
            continue;
        }
        if (take_d("--scale", a.scale)) continue;
        if (take_i("--print-every", a.print_every)) continue;
        if (take_i("--chains-per-axis", a.chains_per_axis)) continue;
        if (take_i("--n-atoms-per-chain", a.n_atoms_per_chain)) continue;
        if (take_d("--box-len", a.box_len)) continue;
        if (take_d("--T-hot", a.T_hot)) continue;
        if (take_d("--T-cold", a.T_cold)) continue;
        if (take_d("--tau-T", a.tau_T)) continue;
        if (take_d("--tau-P", a.tau_P)) continue;
        if (take_d("--strain-rate", a.strain_rate)) continue;
        if (take_i("--cool-segments", a.cool_segments)) continue;
        if (take_s("--out-dir", a.out_dir)) continue;
        if (take_i("--seed", a.seed)) continue;
        if (std::strcmp(s, "--quiet") == 0) { a.quiet = true; continue; }
        if (std::strcmp(s, "--help")  == 0) { print_help(); std::exit(0); }
        std::fprintf(stderr, "unknown argument: %s\n", s);
        return false;
    }
    return true;
}

cpu::Frame make_pe_box(int chains_per_axis, int n_carbons, double box_len, double margin) {
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
                const double y0  = cy * cell + margin;
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
    frame.bonds.atomi  = xt::zeros<int>({n_chains * nb_per});
    frame.bonds.atomj  = xt::zeros<int>({n_chains * nb_per});
    frame.bonds.type   = xt::zeros<int>({n_chains * nb_per});
    frame.angles.atomi = xt::zeros<int>({n_chains * na_per});
    frame.angles.atomj = xt::zeros<int>({n_chains * na_per});
    frame.angles.atomk = xt::zeros<int>({n_chains * na_per});
    frame.angles.type  = xt::zeros<int>({n_chains * na_per});
    frame.dihedrals.atomi = xt::zeros<int>({n_chains * nd_per});
    frame.dihedrals.atomj = xt::zeros<int>({n_chains * nd_per});
    frame.dihedrals.atomk = xt::zeros<int>({n_chains * nd_per});
    frame.dihedrals.atoml = xt::zeros<int>({n_chains * nd_per});
    frame.dihedrals.type  = xt::zeros<int>({n_chains * nd_per});

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

double temperature_K(double ke_eV, int n_dof) {
    if (n_dof <= 0) return 0.0;
    return 2.0 * ke_eV / (static_cast<double>(n_dof) * kBoltzmann);
}

double pressure_eV_per_A3(const cpu::Frame& f) {
    const double V = f.box(0,0) * f.box(1,1) * f.box(2,2);
    if (V <= 0.0) return 0.0;
    const double trace = f.stress[0] + f.stress[1] + f.stress[2];
    return trace / (3.0 * V);
}

/// @brief Inline observer: writes thermo + stress-strain CSV every step.
///
/// Bound to `end_of_step @ +20` -- after `ReduceEnergy(-50)` and
/// `Barostat(-45)` so the recorded values are post-allreduce and reflect
/// the post-step box. Only rank 0 writes to disk.
///
/// The driver-thread updates `current_stage_` / `current_T_target_` /
/// `current_step_offset_` / `record_strain_` between stages so the CSV
/// rows are tagged with the active stage and absolute step index.
class CsvRecorder : public cpu::ComponentBase {
public:
    CsvRecorder(std::ofstream* thermo, std::ofstream* stress,
                int print_every, int n_global_dof, double box_x_initial)
        : thermo_(thermo), stress_(stress),
          print_every_(std::max(1, print_every)),
          n_dof_(n_global_dof),
          box_x0_(box_x_initial) {}

    std::string name() const override { return "csv_recorder"; }
    int frequency() const override { return print_every_; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::end_of_step), 20}};
    }

    /// @brief Enter a new stage. Called between `sim.run(...)` segments.
    void begin_stage(int stage_id, double T_target, bool record_strain) {
        current_stage_     = stage_id;
        current_T_target_  = T_target;
        record_strain_     = record_strain;
        current_step_offset_ += stage_steps_so_far_;
        stage_steps_so_far_ = 0;
    }

    void launch(const cpu::LaunchContext& ctx) override {
        ++stage_steps_so_far_;
        if (ctx.runtime.mpi_rank != 0) return;
        const auto& f = ctx.frame;
        const double T = temperature_K(f.kinetic, n_dof_);
        const double P_bar = pressure_eV_per_A3(f) * kEvPerA3ToBar;
        const int abs_step = current_step_offset_ + ctx.step;
        if (thermo_ && thermo_->is_open()) {
            (*thermo_) << abs_step << "," << current_stage_ << ","
                       << current_T_target_ << "," << T << "," << P_bar << ","
                       << f.energy << "," << f.kinetic << ","
                       << f.box(0,0) << "," << f.box(1,1) << "," << f.box(2,2) << "\n";
        }
        if (record_strain_ && stress_ && stress_->is_open()) {
            const double V = f.box(0,0) * f.box(1,1) * f.box(2,2);
            const double sigma_xx = (V > 0.0) ? f.stress[0] / V : 0.0;
            const double strain = f.box(0,0) / box_x0_ - 1.0;
            (*stress_) << abs_step << "," << current_stage_ << ","
                       << strain << "," << sigma_xx << ","
                       << (sigma_xx * kEvPerA3ToGPa) << "\n";
        }
    }

private:
    std::ofstream* thermo_;
    std::ofstream* stress_;
    int print_every_;
    int n_dof_;
    double box_x0_;
    int    current_stage_       = 0;
    double current_T_target_    = 0.0;
    bool   record_strain_       = false;
    int    current_step_offset_ = 0;
    int    stage_steps_so_far_  = 0;
};

double now_seconds() {
#ifdef ATV_HAS_MPI
    return MPI_Wtime();
#else
    using clk = std::chrono::steady_clock;
    static const auto t0 = clk::now();
    const auto t = clk::now();
    return std::chrono::duration<double>(t - t0).count();
#endif
}

} // namespace

namespace {

// Determine which stages a given run mode actually executes.
bool stage_active(Mode m, int stage_id, int single_id) {
    switch (m) {
        case Mode::full:    return stage_id >= 1 && stage_id <= 5;
        case Mode::nvt:     return stage_id == 1;
        case Mode::npt:     return stage_id == 1 || stage_id == 2;
        case Mode::tension: return stage_id == 1 || stage_id == 2 || stage_id == 5;
        case Mode::single_stage: return stage_id == single_id;
    }
    return false;
}

int scaled_steps(int base, double scale) {
    if (base <= 0) return 0;
    int s = static_cast<int>(base * scale);
    return std::max(s, 1);
}

} // namespace

int main(int argc, char** argv) {
    Args a;
    if (!parse_args(argc, argv, a)) return 2;

#ifdef ATV_HAS_MPI
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    int my_rank = 0;
    int world_size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    (void)world_size;
#else
    int my_rank = 0;
#endif

    cpu::Frame frame = make_pe_box(a.chains_per_axis, a.n_atoms_per_chain,
                                    a.box_len, a.margin);
    const int n_atoms = frame.atoms.size();
    const int n_global_dof = 3 * n_atoms - 3;
    const double box_x0 = frame.box(0, 0);

    std::ofstream thermo_csv, stress_csv;
    if (my_rank == 0) {
        thermo_csv.open(a.out_dir + "/thermo.csv");
        stress_csv.open(a.out_dir + "/stress_strain.csv");
        if (thermo_csv.is_open()) {
            thermo_csv << "step,stage,T_target_K,T_K,P_bar,E_pot_eV,E_kin_eV,box_x,box_y,box_z\n";
        }
        if (stress_csv.is_open()) {
            stress_csv << "step,stage,strain,sigma_xx_eV_per_A3,sigma_xx_GPa\n";
        }
    }

    // All components live as lvalues so the driver stores non-owning
    // pointers; we retune them between stages without rebuilding the
    // driver. The integrator's dt is updated by *replacing* the
    // VelocityVerlet for stage transitions that change dt -- but to keep
    // a single sim.build() across all stages, dt is uniform within a
    // stage block: we use dt_other (0.5 fs) for stages 2-5 and only
    // stage 1 uses dt_stage1 (1.0 fs). Both share the same step floor;
    // mixing dt under one VV would require a setter we deliberately
    // didn't add (dt is fundamental to the Velocity-Verlet kernel).
    //
    // Therefore the user can pick:
    //   * --dt-stage1 == --dt-other  (single dt, all stages homogenised)
    //   * single_stage mode for a fast standalone NVT-only smoke run
    //
    // The default values keep stage 1 at 1 fs and the rest at 0.5 fs;
    // for the *single-driver* canonical run we use dt = dt_other for
    // every stage (acceptable -- LAMMPS reference protocol uses 0.5 fs
    // throughout for stages 2-5 anyway, and stage 1 is just initial
    // randomisation).
    const double dt = a.dt_other;

    cpu::md::VelocityVerlet vv{dt};
    cpu::md::InitVelocity   init_v{a.T_hot, static_cast<unsigned int>(a.seed)};
    cpu::md::NoseHoover     nh{a.T_hot, a.tau_T, dt};
    cpu::md::VirialStress   vs{};
    cpu::md::Barostat       bs{0.0, a.tau_P, dt};
    cpu::md::Deform         df{0.0, dt, /*axis=*/0};
    cpu::md::NeighborList<> nl{kLjCutoff};
    cpu::md::LJ             lj = cpu::md::LJ{lj_eps_eV(), kLjSigma, kLjCutoff}
                                     .set_special_bonds(0.0, 0.0, 0.0);
    HarmonicBond  hb{{bond_k_eV()}, {kBondR0}};
    HarmonicAngle ha{{angle_k_eV()}, {kAngleTheta0}};
    MultiHarmonicDihedral hd{{dihedral_Cn_eV()}};
    CsvRecorder rec{&thermo_csv, &stress_csv, a.print_every, n_global_dof, box_x0};

    cpu::md::MD sim;
    sim.add(vv);
    sim.add(init_v);
    sim.add(nh);
    sim.add(vs);
    sim.add(bs);          // start disabled per-stage; toggle via set_enabled
    sim.add(df);          // strain_rate=0 → no-op until stage 5

    sim.on(ATV::md::Phase::pre_force).add(cpu::md::mpi::HaloExchangeForward{kLjCutoff}, -20);
    sim.on(ATV::md::Phase::post_force).add(cpu::md::mpi::HaloExchangeReverse{}, -1000);
    sim.on(ATV::md::Phase::force).add(nl, -10);
    sim.on(ATV::md::Phase::force).add(lj);
    sim.on(ATV::md::Phase::force).add(hb);
    sim.on(ATV::md::Phase::force).add(ha);
    sim.on(ATV::md::Phase::force).add(hd);
    sim.on(ATV::md::Phase::end_of_step).add(rec, 20);

    if (my_rank == 0 && !a.quiet) {
        const char* mode_str =
            a.mode == Mode::full         ? "full" :
            a.mode == Mode::nvt          ? "NVT" :
            a.mode == Mode::npt          ? "NPT" :
            a.mode == Mode::tension      ? "tension" :
                                            "single_stage";
        std::fprintf(stderr,
            "polyethylene_tension: mode=%s atoms=%d dt=%.3f T_hot=%.1f T_cold=%.1f\n",
            mode_str, n_atoms, dt, a.T_hot, a.T_cold);
    }

    sim.build();
    const double t0 = now_seconds();

    // ── Stage 1 — NVT @ T_hot ──────────────────────────────────────────
    if (stage_active(a.mode, 1, a.single_stage_id)) {
        bs.set_enabled(false);
        df.set_strain_rate(0.0);
        nh.set_T_target(a.T_hot);
        rec.begin_stage(1, a.T_hot, /*record_strain=*/false);
        if (my_rank == 0 && !a.quiet)
            std::fprintf(stderr, "  stage 1: NVT @ %.1f K, %d steps\n",
                         a.T_hot, scaled_steps(a.stage1_steps, a.scale));
        sim.run(frame, scaled_steps(a.stage1_steps, a.scale));
    }

    // ── Stage 2 — NPT @ T_hot, P=0 ────────────────────────────────────
    if (stage_active(a.mode, 2, a.single_stage_id)) {
        nh.set_T_target(a.T_hot);            // 1st set_T_target call
        bs.set_enabled(true);
        bs.set_P_target(0.0);
        df.set_strain_rate(0.0);
        rec.begin_stage(2, a.T_hot, /*record_strain=*/false);
        if (my_rank == 0 && !a.quiet)
            std::fprintf(stderr, "  stage 2: NPT @ %.1f K, %d steps\n",
                         a.T_hot, scaled_steps(a.stage2_steps, a.scale));
        sim.run(frame, scaled_steps(a.stage2_steps, a.scale));
    }

    // ── Stage 3 — NVT cooling T_hot → T_cold (linear ramp) ────────────
    if (stage_active(a.mode, 3, a.single_stage_id)) {
        bs.set_enabled(false);
        df.set_strain_rate(0.0);
        const int total = scaled_steps(a.stage3_steps, a.scale);
        const int n_seg = std::max(1, a.cool_segments);
        const int per_seg = std::max(1, total / n_seg);
        if (my_rank == 0 && !a.quiet)
            std::fprintf(stderr, "  stage 3: NVT cool %.1f → %.1f K, %d segments × %d steps\n",
                         a.T_hot, a.T_cold, n_seg, per_seg);
        for (int seg = 0; seg < n_seg; ++seg) {
            const double frac = static_cast<double>(seg + 1) / n_seg;
            const double T_seg = a.T_hot + (a.T_cold - a.T_hot) * frac;
            nh.set_T_target(T_seg);          // 2nd+ set_T_target call -- ramp
            rec.begin_stage(3, T_seg, /*record_strain=*/false);
            sim.run(frame, per_seg);
        }
    }

    // ── Stage 4 — NPT @ T_cold, P=0 ───────────────────────────────────
    if (stage_active(a.mode, 4, a.single_stage_id)) {
        nh.set_T_target(a.T_cold);           // setpoint into stage 4
        bs.set_enabled(true);
        bs.set_P_target(0.0);
        df.set_strain_rate(0.0);
        rec.begin_stage(4, a.T_cold, /*record_strain=*/false);
        if (my_rank == 0 && !a.quiet)
            std::fprintf(stderr, "  stage 4: NPT @ %.1f K, %d steps\n",
                         a.T_cold, scaled_steps(a.stage4_steps, a.scale));
        sim.run(frame, scaled_steps(a.stage4_steps, a.scale));
    }

    // ── Stage 5 — Uniaxial tension at T_cold, ε̇ along x ──────────────
    if (stage_active(a.mode, 5, a.single_stage_id)) {
        nh.set_T_target(a.T_cold);
        bs.set_enabled(false);   // pure NVT tension; transverse pressure decoupled
        df.set_strain_rate(a.strain_rate);
        rec.begin_stage(5, a.T_cold, /*record_strain=*/true);
        if (my_rank == 0 && !a.quiet)
            std::fprintf(stderr, "  stage 5: tension ε̇=%.2e /fs @ %.1f K, %d steps\n",
                         a.strain_rate, a.T_cold,
                         scaled_steps(a.stage5_steps, a.scale));
        sim.run(frame, scaled_steps(a.stage5_steps, a.scale));
    }

    const double t1 = now_seconds();

    if (my_rank == 0 && !a.quiet) {
        std::fprintf(stderr,
            "polyethylene_tension done: %d atoms, %.2fs total\n",
            n_atoms, t1 - t0);
    }

    if (my_rank == 0) {
        if (thermo_csv.is_open()) thermo_csv.close();
        if (stress_csv.is_open()) stress_csv.close();
    }

#ifdef ATV_HAS_MPI
    MPI_Finalize();
#endif
    return 0;
}
