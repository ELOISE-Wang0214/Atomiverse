/// @file benchmarks/md/polyethylene_tension.cpp
/// @brief MSU-style amorphous-polyethylene uniaxial-tension benchmark.
///
/// Single-driver simplification of the 5-stage protocol from
/// `parallel-md-engine.md` Verification §7. Scaffolds the benchmark
/// machinery -- in-code amorphous-PE setup, force-field components, CSV
/// thermo / stress recorders -- but uses ONE `sim.run()` invocation rather
/// than five back-to-back `sim`s. The reason: under a fresh `Backend` per
/// stage the multi-rank `DomainDecomposition` re-scans every row of the
/// (possibly stale-tailed) atom array; duplicates emerge when multiple rows
/// happen to land in the same rank's slab. Running a single sim keeps the
/// partition state internally consistent.
///
/// As a benchmark deliverable this is structurally complete:
///   * builds at `ATV_BUILD_BENCHMARKS=ON`,
///   * runs at any rank count via `mpirun`,
///   * emits `thermo.csv` and `stress_strain.csv` in the configured format,
///   * is consumable by `tools/validate_polyethylene.py`.
///
/// The full multi-stage protocol (heating / cooling / NPT relax / tension)
/// requires either a `NoseHoover::set_T_target()` setter (not yet present)
/// or persistent-backend driver composition. Tracked under
/// `msu-pe-validation.md` Phase C.
///
/// Output (rank 0 only) under `--out-dir`:
///   * `thermo.csv`: `step,T_K,P_bar,E_pot_eV,E_kin_eV,box_x,box_y,box_z`
///   * `stress_strain.csv`: `strain,sigma_xx_eV_per_A3,sigma_xx_GPa`
///
/// CLI knobs select the simulation regime — `--mode tension` adds Deform +
/// strain rate, `--mode npt` runs NPT only, `--mode nvt` is plain NVT.
///
/// @see .claude/specs/parallel-md-engine.md  Verification §7
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

#include "atomiverse/cpu/bonded/harmonic_angle.h"
#include "atomiverse/cpu/bonded/harmonic_bond.h"
#include "atomiverse/cpu/bonded/multi_harmonic_dihedral.h"
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
constexpr double kBondR0         = 1.529;
constexpr double kAngleK_kcal    = 58.35;
constexpr double kAngleTheta0    = 112.7 * kPI / 180.0;
constexpr double kV1_kcal        = 1.411;
constexpr double kV2_kcal        = -0.271;
constexpr double kV3_kcal        = 3.145;
constexpr double kLjEps_kcal     = 0.112;
constexpr double kLjSigma        = 4.01;
constexpr double kLjCutoff       = 10.5;

constexpr double kEvPerA3ToGPa = 160.21766;
constexpr double kEvPerA3ToBar = 1.602176634e6;

double bond_k_eV()  { return kBondK_kcal  * kKcalPerMolToEv; }
double angle_k_eV() { return kAngleK_kcal * kKcalPerMolToEv; }
double lj_eps_eV()  { return kLjEps_kcal  * kKcalPerMolToEv; }

std::vector<double> dihedral_Cn_eV() {
    const double v1 = kV1_kcal * kKcalPerMolToEv;
    const double v2 = kV2_kcal * kKcalPerMolToEv;
    const double v3 = kV3_kcal * kKcalPerMolToEv;
    return {
        v1 / 2.0 + v2 + v3 / 2.0,
        v1 / 2.0 - 1.5 * v3,
        -v2,
         2.0 * v3,
         0.0,
    };
}

enum class Mode { nvt, npt, tension };

struct Args {
    Mode   mode        = Mode::tension;
    int    steps       = 2000;
    int    print_every = 100;
    int    chains_per_axis = 4;
    int    n_atoms_per_chain = 18;
    double box_len     = 60.0;
    double margin      = 1.0;
    int    seed        = 12345;
    double T_target    = 300.0;     ///< thermostat target [K]
    double dt          = 0.5;       ///< [fs]
    double tau_T       = 100.0;     ///< [fs]
    double tau_P       = 1000.0;    ///< [fs]
    double strain_rate = 5e-7;      ///< 1/fs
    std::string out_dir = ".";
    bool   quiet       = false;
};

void print_help() {
    std::printf(
        "Usage: polyethylene_tension [options]\n"
        "  --mode {nvt|npt|tension}  Simulation regime (default tension).\n"
        "  --steps N                 MD steps to run (default 2000).\n"
        "  --print-every N           Thermo/stress sample stride (default 100).\n"
        "  --chains-per-axis N       NxNxN chain grid (default 4).\n"
        "  --n-atoms-per-chain N     UA atoms per chain (default 18).\n"
        "  --box-len L               Initial cubic box edge in Å (default 60).\n"
        "  --T K                     Thermostat target in K (default 300).\n"
        "  --dt fs                   Time step in fs (default 0.5).\n"
        "  --tau-T fs                Thermostat time constant (default 100).\n"
        "  --tau-P fs                Barostat time constant (default 1000).\n"
        "  --strain-rate /fs         Engineering strain rate, tension mode (default 5e-7).\n"
        "  --out-dir DIR             Where to write CSV outputs (default '.').\n"
        "  --seed S                  RNG seed for InitVelocity (default 12345).\n"
        "  --quiet                   Suppress per-stage progress lines.\n"
        "  --help                    Show this message.\n");
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
            if      (v == "nvt")     a.mode = Mode::nvt;
            else if (v == "npt")     a.mode = Mode::npt;
            else if (v == "tension") a.mode = Mode::tension;
            else { std::fprintf(stderr, "unknown mode: %s\n", v.c_str()); return false; }
            continue;
        }
        if (take_i("--steps", a.steps)) continue;
        if (take_i("--print-every", a.print_every)) continue;
        if (take_i("--chains-per-axis", a.chains_per_axis)) continue;
        if (take_i("--n-atoms-per-chain", a.n_atoms_per_chain)) continue;
        if (take_d("--box-len", a.box_len)) continue;
        if (take_d("--T", a.T_target)) continue;
        if (take_d("--dt", a.dt)) continue;
        if (take_d("--tau-T", a.tau_T)) continue;
        if (take_d("--tau-P", a.tau_P)) continue;
        if (take_d("--strain-rate", a.strain_rate)) continue;
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
class CsvRecorder : public cpu::ComponentBase {
public:
    CsvRecorder(std::ofstream* thermo, std::ofstream* stress,
                int print_every, int n_global_dof,
                bool record_strain, double box_x_initial)
        : thermo_(thermo), stress_(stress),
          print_every_(std::max(1, print_every)),
          n_dof_(n_global_dof),
          record_strain_(record_strain),
          box_x0_(box_x_initial) {}

    std::string name() const override { return "csv_recorder"; }
    int frequency() const override { return print_every_; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::end_of_step), 20}};
    }

    void launch(const cpu::LaunchContext& ctx) override {
        if (ctx.runtime.mpi_rank != 0) return;
        const auto& f = ctx.frame;
        const double T = temperature_K(f.kinetic, n_dof_);
        const double P_bar = pressure_eV_per_A3(f) * kEvPerA3ToBar;
        if (thermo_ && thermo_->is_open()) {
            (*thermo_) << ctx.step << "," << T << "," << P_bar << ","
                       << f.energy << "," << f.kinetic << ","
                       << f.box(0,0) << "," << f.box(1,1) << "," << f.box(2,2) << "\n";
        }
        if (record_strain_ && stress_ && stress_->is_open()) {
            const double V = f.box(0,0) * f.box(1,1) * f.box(2,2);
            const double sigma_xx = (V > 0.0) ? f.stress[0] / V : 0.0;
            const double strain = f.box(0,0) / box_x0_ - 1.0;
            (*stress_) << strain << "," << sigma_xx << ","
                       << (sigma_xx * kEvPerA3ToGPa) << "\n";
        }
    }

private:
    std::ofstream* thermo_;
    std::ofstream* stress_;
    int print_every_;
    int n_dof_;
    bool record_strain_;
    double box_x0_;
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
            thermo_csv << "step,T_K,P_bar,E_pot_eV,E_kin_eV,box_x,box_y,box_z\n";
        }
        if (stress_csv.is_open()) {
            stress_csv << "strain,sigma_xx_eV_per_A3,sigma_xx_GPa\n";
        }
    }

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{a.dt});
    sim.add(cpu::md::InitVelocity{a.T_target, static_cast<unsigned int>(a.seed)});
    sim.add(cpu::md::NoseHoover{a.T_target, a.tau_T, a.dt});
    sim.add(cpu::md::VirialStress{});

    if (a.mode == Mode::npt || a.mode == Mode::tension) {
        sim.add(cpu::md::Barostat{0.0, a.tau_P, a.dt});
    }
    if (a.mode == Mode::tension) {
        sim.add(cpu::md::Deform{a.strain_rate, a.dt, /*axis=*/0});
    }

    sim.on(ATV::md::Phase::pre_force).add(cpu::md::mpi::HaloExchangeForward{kLjCutoff}, -20);
    sim.on(ATV::md::Phase::post_force).add(cpu::md::mpi::HaloExchangeReverse{}, -1000);
    sim.on(ATV::md::Phase::force).add(cpu::md::NeighborList<>{kLjCutoff}, -10);
    sim.on(ATV::md::Phase::force).add(cpu::md::LJ{lj_eps_eV(), kLjSigma, kLjCutoff}
                                       .set_special_bonds(0.0, 0.0, 0.0));
    sim.on(ATV::md::Phase::force).add(cpu::bonded::HarmonicBond{{bond_k_eV()}, {kBondR0}});
    sim.on(ATV::md::Phase::force).add(cpu::bonded::HarmonicAngle{{angle_k_eV()}, {kAngleTheta0}});
    sim.on(ATV::md::Phase::force).add(cpu::bonded::MultiHarmonicDihedral{{dihedral_Cn_eV()}});

    // CSV recorder runs at end_of_step @ +20, after ReduceEnergy / Barostat.
    sim.on(ATV::md::Phase::end_of_step)
       .add(CsvRecorder{&thermo_csv, &stress_csv,
                        a.print_every, n_global_dof,
                        /*record_strain=*/a.mode == Mode::tension, box_x0},
            20);

    if (my_rank == 0 && !a.quiet) {
        const char* mode_str =
            a.mode == Mode::nvt ? "NVT" :
            a.mode == Mode::npt ? "NPT" : "tension";
        std::fprintf(stderr,
            "polyethylene_tension: mode=%s atoms=%d steps=%d dt=%.3f T=%.1fK\n",
            mode_str, n_atoms, a.steps, a.dt, a.T_target);
    }

    sim.build();
    const double t0 = now_seconds();
    sim.run(frame, a.steps);
    const double t1 = now_seconds();

    if (my_rank == 0 && !a.quiet) {
        std::fprintf(stderr,
            "polyethylene_tension done: %d atoms, %d steps, %.2fs total (%.3f ms/step)\n",
            n_atoms, a.steps, t1 - t0, 1000.0 * (t1 - t0) / std::max(1, a.steps));
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
