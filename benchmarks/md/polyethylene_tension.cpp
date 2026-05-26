/// @file benchmarks/md/polyethylene_tension.cpp
/// @brief MSU-style amorphous-polyethylene uniaxial-tension benchmark.
///
/// Drives the 5-stage MSU equilibration + uniaxial-tension protocol so
/// that each stage uses the *exact* integrator, ensemble and timestep of
/// the matching stage in `benchmarks/md/reference/polyethylene/lammps.in`.
/// The benchmark mirrors the LAMMPS deck stage-for-stage — it is the
/// precondition for the Phase D cross-engine comparison to be meaningful.
///
/// ── lammps.in → Atomiverse stage map (the binding contract) ──────────
///
///   | Stage | lammps.in                       | components                  | dt    |
///   |-------|---------------------------------|-----------------------------|-------|
///   | 1     | nve/limit 0.05 + langevin 500   | NVELimit+FixLangevin+InitVel| 1.0 fs|
///   | 2     | npt temp 500 500 iso 0          | isotropic NoseHooverNPT     | 0.5 fs|
///   | 3     | npt temp 500 100 iso 0          | ramped NoseHooverNPT 500→100| 0.5 fs|
///   | 4     | npt temp 100 100 iso 0          | isotropic NoseHooverNPT     | 0.5 fs|
///   | 5     | npt temp 100 100 y/z 0 + deform | aniso NoseHooverNPT + Deform| 1.0 fs|
///
/// ── Per-stage driver rebuild ─────────────────────────────────────────
///
/// `dt` is a constructor parameter of every integrator kernel and is
/// fundamental to the Velocity-Verlet / MTK splitting — there is
/// deliberately no `dt` setter. Stages 1 and 5 run at 1.0 fs, stages
/// 2-4 at 0.5 fs, so a single `sim.build()` spanning all stages cannot
/// honour the deck. The benchmark therefore constructs and `build()`s a
/// **fresh `cpu::md::MD` per stage**, runs that stage to completion,
/// then discards the driver. The `cpu::Frame` (positions, velocities,
/// box) is the carry-over state and persists across stages — correct,
/// since the Frame is the device-resident simulation state and the
/// driver owns only the loop. The shared `CsvRecorder` /
/// `TrajectoryRecorder` / optional `RDF` are lvalues re-registered into
/// every per-stage driver, so the CSV step column stays globally
/// monotonic across stage boundaries.
///
/// Stage 1 uses `NVELimit` rather than plain Velocity-Verlet because an
/// unrelaxed amorphous PE start-state carries severe atom overlap
/// (E_pair ~3.4e8 at step 0); the per-step displacement cap tames those
/// transient force spikes instead of letting the integrator blow up.
///
/// Stage 3 is a *single* `NoseHooverNPT` with `T_start=500, T_stop=100`
/// and a `set_run_span` spanning the stage — the continuous temperature
/// ramp interpolates the setpoint internally exactly as LAMMPS
/// `fix npt temp 500 100`. There is no cooling micro-segment loop.
///
/// `--mode {full|nvt|npt|tension|stage:<id>}` selects the run shape:
///
///   * `full`     -- run all 5 stages (default; the canonical benchmark);
///   * `nvt`      -- stage 1 only (kept for back-compat with
///                  `tests/mpi/test_pe_msu_pipeline_n8.cpp`);
///   * `npt`      -- stages 1+2;
///   * `tension`  -- stages 1+2+5 (skip cooling/equilibration);
///   * `stage:<id>` -- run a single stage by 1-indexed id.
///
/// Output (rank 0 only) under `--out-dir`:
///   * `thermo.csv`: 17-column extended schema (Phase G.2 + G.3,
///     msu-pe-cuda-components-01) — `step,stage,T_target_K,T_K,P_bar,
///     Pxx,Pyy,Pzz,E_pot_eV,E_pair_eV,E_bond_eV,E_angle_eV,E_dihedral_eV,
///     E_kin_eV,box_x,box_y,box_z`. Pxx/Pyy/Pzz are genuine per-axis
///     virials [bar]; the 4 per-term energies decompose E_pot_eV.
///   * `stress_strain.csv`: `step,stage,strain,sigma_xx_eV_per_A3,sigma_xx_GPa`
///   * `trajectory.xyz` (only with `--traj-every N`): plain-XYZ frames.
///   * `rdf.dat` (only with `--rdf RMAX,BINS[,STRIDE]`): radial g(r).
///
/// Both optional outputs default off, so a plain run is byte-identical
/// to the pre-flag CSV-only behaviour.
///
/// @see benchmarks/md/reference/polyethylene/lammps.in  (the reference deck)
/// @see .claude/specs/msu-pe-align-02-benchmark-protocol.md
/// @see .claude/specs/msu-pe-validation.md   Phase D

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/lammps_data.h"
#include "atomiverse/cpu/md/deform.h"
#include "atomiverse/cpu/md/fix_langevin.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/cpu/md/neighbor_list.h"
#include "atomiverse/cpu/md/nose_hoover_npt.h"
#include "atomiverse/cpu/md/nve_limit.h"
#include "atomiverse/cpu/observer/rdf.h"
#include "atomiverse/cpu/potentials/harmonic_angle.h"
#include "atomiverse/cpu/potentials/harmonic_bond.h"
#include "atomiverse/cpu/potentials/multi_harmonic_dihedral.h"
#include "atomiverse/cpu/recorder/trajectory.h"
#include "atomiverse/facade/stage.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#else
#include <chrono>
#endif

#ifdef ATV_BUILD_CUDA
#include "atomiverse/core/system_scalars.h"
#include "atomiverse/cuda/async/record_base.cuh"
#include "atomiverse/cuda/md/deform.cuh"
#include "atomiverse/cuda/md/fix_langevin.cuh"
#include "atomiverse/cuda/md/init.cuh"
#include "atomiverse/cuda/md/md.cuh"
#include "atomiverse/cuda/md/nose_hoover_npt.cuh"
#include "atomiverse/cuda/md/nve_limit.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cuda/observer/stress.cuh"
#include "atomiverse/cuda/potentials/angle_harmonic.cuh"
#include "atomiverse/cuda/potentials/bond_harmonic.cuh"
#include "atomiverse/cuda/potentials/dihedral_multi_harmonic.cuh"
#include "atomiverse/cuda/potentials/pair_lj_cut.cuh"
#include "atomiverse/cuda/types.cuh"
#include <cuda_runtime.h>
#endif

using namespace ATV;

namespace {

// ── Force-field constants (Dreiding-style UA polyethylene) ──────────
constexpr double kKcalPerMolToEv = 1.0 / 23.060547830619026;
constexpr double kBondK_kcal = 350.0;
constexpr double kBondR0 = 1.53;
constexpr double kAngleK_kcal = 60.0;
constexpr double kAngleTheta0 = 109.5 * kPI / 180.0;
constexpr double kLjEps_kcal = 0.112;
constexpr double kLjSigma = 4.01;
constexpr double kLjCutoff = 10.5;

constexpr double kEvPerA3ToGPa = 160.21766;
constexpr double kEvPerA3ToBar = 1.602176634e6;

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

enum class Mode { full, nvt, npt, tension, single_stage };

/// @brief Backend selector — CPU (default) or CUDA (single-GPU, requires
/// the binary to be compiled with `ATV_BUILD_CUDA`). msu-pe-cuda-components-06.
enum class Backend { cpu, cuda };

struct Args {
    Mode mode = Mode::full;
    Backend backend = Backend::cpu;
    int single_stage_id = 1; ///< 1..5 when mode == single_stage
    /// Per-stage step counts — byte-aligned with `lammps.in` /
    /// `lammps_msu_paper.in`, the immutable MSU reference protocol, so a
    /// plain run reproduces the reference exactly. `--scale F` multiplies
    /// every stage; pass F < 1 (e.g. 0.05) for a quick smoke run.
    int stage1_steps = 10000;  ///< NVELimit+Langevin @ 500 K  (run 10000)
    int stage2_steps = 50000;  ///< NPT @ 500 K                (run 50000)
    int stage3_steps = 50000;  ///< NPT cool 500 → 100 K       (run 50000)
    int stage4_steps = 50000;  ///< NPT @ 100 K                (run 50000)
    int stage5_steps = 171800; ///< tension along x            (run 171800)
    double scale = 1.0;        ///< multiplier on each stage's steps

    int print_every = 100;
    int chains_per_axis = 4;
    int n_atoms_per_chain = 18;
    double box_len = 60.0;
    double margin = 1.0;
    int seed = 12345;
    std::string data_path; ///< If non-empty, load this LAMMPS data file
                           ///< instead of generating the chain grid.
    double T_hot = 500.0;  ///< stages 1-3 start setpoint [K]
    double T_cold = 100.0; ///< stages 3-5 final setpoint [K]
    /// Per-stage timesteps [fs], byte-aligned with the `lammps.in`
    /// `timestep` budget: stages 1 and 5 at 1.0 fs, stages 2-4 at
    /// 0.5 fs. `dt` is fundamental to each integrator kernel (no
    /// setter), which is why every stage builds its own driver.
    double dt_stage1 = 1.0;       ///< [fs] -- lammps.in `timestep 1`
    double dt_stage2 = 0.5;       ///< [fs] -- lammps.in `timestep 0.5`
    double dt_stage3 = 0.5;       ///< [fs] -- lammps.in `timestep 0.5`
    double dt_stage4 = 0.5;       ///< [fs] -- lammps.in `timestep 0.5`
    double dt_stage5 = 1.0;       ///< [fs] -- lammps.in `timestep 1`
    double tau_T = 100.0;         ///< [fs]
    double tau_P = 1000.0;        ///< [fs]
    double langevin_damp = 10.0;  ///< [fs] stage-1 Langevin damp — matches
                                  ///  lammps.in `fix langevin 500 500 10`.
                                  ///  Was 100.0 (a weak-coupling default that
                                  ///  doesn't match the reference). At T=500
                                  ///  with NVELimit and 10000 atoms with
                                  ///  initial LJ overlap, damp=100 left the
                                  ///  Langevin too weak to absorb the
                                  ///  repulsive transient and stage 1 blew
                                  ///  up. damp=10 matches the LAMMPS GPU
                                  ///  reference both numerically (Pearson
                                  ///  block-average gate) and in stability.
    double strain_rate =
        1e-5; ///< 1/fs (stage 5) -- lammps.in `fix deform x erate 1e-5`
    std::string out_dir = ".";
    bool quiet = false;
    int traj_every = 0;   ///< trajectory XYZ stride (0 = off)
    std::string rdf_spec; ///< "RMAX,BINS[,STRIDE]" (empty = off)
    // ── Diagnostic flags: disable individual force kernels (CUDA path
    //    only) to isolate which kernel breaks energy conservation under
    //    pure NVE. Always-off in production.
    bool disable_lj = false;
    bool disable_bond = false;
    bool disable_angle = false;
    bool disable_dihedral = false;
    bool disable_langevin = false;
};

void print_help() {
    std::printf(
        "Usage: polyethylene_tension [options]\n"
        "  --mode {full|nvt|npt|tension|stage:N}   Run shape (default full).\n"
        "      full        : all 5 stages, one driver rebuilt per stage.\n"
        "      nvt         : stage 1 only.\n"
        "      npt         : stages 1+2.\n"
        "      tension     : stages 1+2+5 (no cooling/equilibration).\n"
        "      stage:N     : single stage by 1-indexed id (1..5).\n"
        "  --scale F                   Multiply every stage's step count by F "
        "(default 1).\n"
        "  --print-every N             Thermo/stress sample stride (default "
        "100).\n"
        "  --chains-per-axis N         NxNxN chain grid (default 4).\n"
        "  --n-atoms-per-chain N       UA atoms per chain (default 18).\n"
        "  --box-len L                 Initial cubic box edge in Å (default "
        "60).\n"
        "  --T-hot K                   Hot setpoint, stages 1-3 (default "
        "500).\n"
        "  --T-cold K                  Cold setpoint, stages 3-5 (default "
        "100).\n"
        "  --tau-T fs                  Thermostat time constant (default "
        "100).\n"
        "  --tau-P fs                  Barostat time constant (default 1000).\n"
        "  --langevin-damp fs          Stage-1 Langevin damp (default 10, "
        "matches LAMMPS).\n"
        "  --strain-rate /fs           Engineering strain rate stage 5 "
        "(default 1e-5).\n"
        "  --disable-{lj|bond|angle|dihedral|langevin}\n"
        "                              Diagnostic: drop the named force /\n"
        "                              thermostat term from the CUDA stack.\n"
        "  --out-dir DIR               Where to write CSV outputs (default "
        "'.').\n"
        "  --data PATH                 Load this LAMMPS data file as the\n"
        "                              initial system instead of generating a\n"
        "                              chain grid (e.g. the MSU PE deck).\n"
        "  --traj-every N              Append a plain-XYZ trajectory frame to\n"
        "                              <out-dir>/trajectory.xyz every N steps\n"
        "                              (and always the first frame). 0 = off.\n"
        "  --rdf RMAX,BINS[,STRIDE]    Accumulate a radial distribution\n"
        "                              function into <out-dir>/rdf.dat. "
        "STRIDE\n"
        "                              omitted = sample every step. Empty = "
        "off.\n"
        "  --seed S                    RNG seed for InitVelocity (default "
        "12345).\n"
        "  --quiet                     Suppress per-stage progress lines.\n"
        "  --help                      Show this message.\n");
}

bool parse_args(int argc, char** argv, Args& a) {
    for (int i = 1; i < argc; ++i) {
        const char* s = argv[i];
        auto take_d = [&](const char* name, double& dst) {
            if (std::strcmp(s, name) != 0)
                return false;
            if (i + 1 >= argc)
                std::exit(2);
            dst = std::atof(argv[++i]);
            return true;
        };
        auto take_i = [&](const char* name, int& dst) {
            if (std::strcmp(s, name) != 0)
                return false;
            if (i + 1 >= argc)
                std::exit(2);
            dst = std::atoi(argv[++i]);
            return true;
        };
        auto take_s = [&](const char* name, std::string& dst) {
            if (std::strcmp(s, name) != 0)
                return false;
            if (i + 1 >= argc)
                std::exit(2);
            dst = argv[++i];
            return true;
        };
        if (std::strcmp(s, "--backend") == 0) {
            if (i + 1 >= argc)
                std::exit(2);
            const std::string v = argv[++i];
            if (v == "cpu")
                a.backend = Backend::cpu;
            else if (v == "cuda")
                a.backend = Backend::cuda;
            else {
                std::fprintf(stderr, "unknown backend: %s (expected cpu|cuda)\n",
                             v.c_str());
                return false;
            }
            continue;
        }
        if (std::strcmp(s, "--mode") == 0) {
            if (i + 1 >= argc)
                std::exit(2);
            const std::string v = argv[++i];
            if (v == "full")
                a.mode = Mode::full;
            else if (v == "nvt")
                a.mode = Mode::nvt;
            else if (v == "npt")
                a.mode = Mode::npt;
            else if (v == "tension")
                a.mode = Mode::tension;
            else if (v.rfind("stage:", 0) == 0) {
                a.mode = Mode::single_stage;
                a.single_stage_id = std::atoi(v.c_str() + 6);
                if (a.single_stage_id < 1 || a.single_stage_id > 5) {
                    std::fprintf(stderr, "stage id must be 1..5: %s\n",
                                 v.c_str());
                    return false;
                }
            } else {
                std::fprintf(stderr, "unknown mode: %s\n", v.c_str());
                return false;
            }
            continue;
        }
        if (take_d("--scale", a.scale))
            continue;
        if (take_i("--print-every", a.print_every))
            continue;
        if (take_i("--chains-per-axis", a.chains_per_axis))
            continue;
        if (take_i("--n-atoms-per-chain", a.n_atoms_per_chain))
            continue;
        if (take_d("--box-len", a.box_len))
            continue;
        if (take_d("--T-hot", a.T_hot))
            continue;
        if (take_d("--T-cold", a.T_cold))
            continue;
        if (take_d("--tau-T", a.tau_T))
            continue;
        if (take_d("--tau-P", a.tau_P))
            continue;
        if (take_d("--langevin-damp", a.langevin_damp))
            continue;
        if (take_d("--strain-rate", a.strain_rate))
            continue;
        if (take_s("--out-dir", a.out_dir))
            continue;
        if (take_s("--data", a.data_path))
            continue;
        if (take_i("--traj-every", a.traj_every))
            continue;
        if (take_s("--rdf", a.rdf_spec))
            continue;
        if (take_i("--seed", a.seed))
            continue;
        if (std::strcmp(s, "--quiet") == 0) {
            a.quiet = true;
            continue;
        }
        if (std::strcmp(s, "--disable-lj") == 0) {
            a.disable_lj = true;
            continue;
        }
        if (std::strcmp(s, "--disable-bond") == 0) {
            a.disable_bond = true;
            continue;
        }
        if (std::strcmp(s, "--disable-angle") == 0) {
            a.disable_angle = true;
            continue;
        }
        if (std::strcmp(s, "--disable-dihedral") == 0) {
            a.disable_dihedral = true;
            continue;
        }
        if (std::strcmp(s, "--disable-langevin") == 0) {
            a.disable_langevin = true;
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

double temperature_K(double ke_eV, int n_dof) {
    if (n_dof <= 0)
        return 0.0;
    return 2.0 * ke_eV / (static_cast<double>(n_dof) * kBoltzmann);
}

double pressure_eV_per_A3(const cpu::Frame& f) {
    const double V = f.box(0, 0) * f.box(1, 1) * f.box(2, 2);
    if (V <= 0.0)
        return 0.0;
    // LAMMPS `compute pressure`: P = (2·KE + virial_trace) / (3V).
    // f.scalars.stress carries the pair+bonded virial; the kinetic part
    // (trace = 2·KE) lives in f.scalars.stress_kinetic. Barostat /
    // NoseHooverNPT already combine the two — this reporter must match.
    const double trace =
        f.scalars.stress[0] + f.scalars.stress[1] + f.scalars.stress[2]
        + f.scalars.stress_kinetic[0] + f.scalars.stress_kinetic[1]
        + f.scalars.stress_kinetic[2];
    return trace / (3.0 * V);
}

/// @brief Per-axis pressure tensor diagonal [bar], Voigt Pxx/Pyy/Pzz.
///
/// LAMMPS `compute pressure` (`compute_pressure.cpp:257`):
/// `P_αα = (Σ_i m_i v_iα² + virial_αα) / V`. Same sign convention as
/// `pressure_eV_per_A3` (positive ⇔ compressive). `f.scalars.stress`
/// carries the pair+bonded virial diagonal; `f.scalars.stress_kinetic`
/// carries Σ m·v_α². Reports in bar via `kEvPerA3ToBar`. Phase G.2
/// reporter (msu-pe-cuda-components-01).
void pressure_per_axis_bar(const cpu::Frame& f, double Pbar[3]) {
    const double V = f.box(0, 0) * f.box(1, 1) * f.box(2, 2);
    if (V <= 0.0) {
        Pbar[0] = Pbar[1] = Pbar[2] = 0.0;
        return;
    }
    for (int k = 0; k < 3; ++k) {
        const double t_axis =
            f.scalars.stress[k] + f.scalars.stress_kinetic[k];
        Pbar[k] = (t_axis / V) * kEvPerA3ToBar;
    }
}

/// @brief Inline observer: writes thermo + stress-strain CSV every step.
///
/// Bound to `end_of_step @ +20` -- after `ReduceEnergy(-50)` and the
/// integrator's `end_of_step` slot so the recorded values are
/// post-allreduce and reflect the post-step box. Only rank 0 writes
/// to disk.
///
/// This recorder is a shared lvalue re-registered into every per-stage
/// driver. `begin_stage()` is called between stages to tag the rows
/// with the active stage and keep the absolute step index globally
/// monotonic across the per-stage driver rebuilds.
class CsvRecorder : public cpu::ComponentBase {
  public:
    CsvRecorder(std::ofstream* thermo, std::ofstream* stress, int print_every,
                int n_global_dof, double box_x_initial)
        : thermo_(thermo), stress_(stress),
          print_every_(std::max(1, print_every)), n_dof_(n_global_dof),
          box_x0_(box_x_initial) {}

    std::string name() const override { return "csv_recorder"; }
    int frequency() const override { return print_every_; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::end_of_step), 20}};
    }

    /// @brief Enter a new stage. Called before each per-stage `sim.run`.
    /// Folds the previous stage's step count into the absolute offset so
    /// the CSV step column stays monotonic across the per-stage driver
    /// rebuilds (each fresh driver restarts its own step counter at 0).
    void begin_stage(int stage_id, double T_target, bool record_strain) {
        current_stage_ = stage_id;
        current_T_target_ = T_target;
        record_strain_ = record_strain;
        // `stage_steps_so_far_` counts launch() invocations, which fire
        // every `print_every_` integration steps. Multiply to recover
        // the actual step count of the completed stage so the absolute
        // step column stays globally monotonic across stages (assumes
        // each stage's step count is a multiple of print_every_, which
        // the MSU protocol satisfies).
        current_step_offset_ += stage_steps_so_far_ * print_every_;
        stage_steps_so_far_ = 0;
        // Re-arm the strain reference: it is snapshotted at the first
        // stress-recorded step so engineering strain is referenced to
        // the box edge at deformation onset, not the sim-initial box.
        if (record_strain)
            strain_ref_set_ = false;
    }

    void launch(const cpu::LaunchContext& ctx) override {
        ++stage_steps_so_far_;
        if (ctx.runtime.mpi_rank != 0)
            return;
        const auto& f = ctx.frame;
        const double T = temperature_K(f.scalars.kinetic, n_dof_);
        const double P_bar = pressure_eV_per_A3(f) * kEvPerA3ToBar;
        double Pbar[3];
        pressure_per_axis_bar(f, Pbar);
        const int abs_step = current_step_offset_ + ctx.step;
        if (thermo_ && thermo_->is_open()) {
            // Phase G.2 + G.3 extended schema (msu-pe-cuda-components-01):
            //   step,stage,T_target_K,T_K,P_bar,
            //   Pxx,Pyy,Pzz,
            //   E_pot_eV,E_pair_eV,E_bond_eV,E_angle_eV,E_dihedral_eV,E_kin_eV,
            //   box_x,box_y,box_z
            // Build the row in a local ostringstream and flush as ONE
            // ofstream insertion — keeps every CSV row a single atomic
            // write to the underlying file (the chained-<< form gets
            // shredded into ~16 tiny stream-buffer writes that can
            // interleave across step / phase / OMP boundaries on
            // shared-fs Lustre, producing ~0.8% malformed rows).
            std::ostringstream row;
            row << abs_step << "," << current_stage_ << ","
                << current_T_target_ << "," << T << "," << P_bar << ","
                << Pbar[0] << "," << Pbar[1] << "," << Pbar[2] << ","
                << f.scalars.energy << "," << f.scalars.energy_pair << ","
                << f.scalars.energy_bond << "," << f.scalars.energy_angle
                << "," << f.scalars.energy_dihedral << ","
                << f.scalars.kinetic << "," << f.box(0, 0) << ","
                << f.box(1, 1) << "," << f.box(2, 2) << "\n";
            (*thermo_) << row.str();
        }
        if (record_strain_ && stress_ && stress_->is_open()) {
            // Snapshot the engineering-strain reference at the first
            // stress-recorded step (deformation onset), mirroring LAMMPS
            // fix deform's L0-at-fix-issue rule. The sim-initial box is
            // wrong here: NPT stages 2-4 compress the box beforehand.
            if (!strain_ref_set_) {
                box_x0_ = f.box(0, 0);
                strain_ref_set_ = true;
            }
            const double V = f.box(0, 0) * f.box(1, 1) * f.box(2, 2);
            // Tensile stress σ_xx = -P_xx with LAMMPS pressure convention
            // P_αα = (kinetic_αα + virial_αα) / V. f.scalars.stress carries
            // the pair+bonded virial; f.scalars.stress_kinetic the kinetic
            // (Σ m·vα·vα). Both must be summed — omitting kinetic shifted
            // σ_xx by ~+0.05 GPa at 100 K (≈ NkT/V per axis), enough to
            // bias the small-strain elastic-modulus fit by several×.
            const double sigma_xx =
                (V > 0.0) ? -(f.scalars.stress[0]
                              + f.scalars.stress_kinetic[0]) / V
                          : 0.0;
            const double strain = f.box(0, 0) / box_x0_ - 1.0;
            std::ostringstream row;
            row << abs_step << "," << current_stage_ << "," << strain << ","
                << sigma_xx << "," << (sigma_xx * kEvPerA3ToGPa) << "\n";
            (*stress_) << row.str();
        }
    }

  private:
    std::ofstream* thermo_;
    std::ofstream* stress_;
    int print_every_;
    int n_dof_;
    double box_x0_;
    int current_stage_ = 0;
    double current_T_target_ = 0.0;
    bool record_strain_ = false;
    bool strain_ref_set_ = false;
    int current_step_offset_ = 0;
    int stage_steps_so_far_ = 0;
};

/// @brief Inline observer: appends a plain-XYZ trajectory frame.
///
/// Writes the very first frame and then one frame every `traj_every`
/// steps, so a run of *any* length yields at least one frame. Bound to
/// `end_of_step @ +25` -- after `CsvRecorder(+20)` so the recorded
/// geometry reflects the post-step box. Only rank 0 writes. Like
/// `CsvRecorder`, this is a shared lvalue re-registered into every
/// per-stage driver so the frame index is continuous.
class TrajectoryRecorder : public cpu::ComponentBase {
  public:
    TrajectoryRecorder(std::ofstream* xyz, int traj_every)
        : xyz_(xyz), traj_every_(std::max(1, traj_every)) {}

    std::string name() const override { return "trajectory_recorder"; }
    int frequency() const override { return 1; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::end_of_step), 25}};
    }

    void launch(const cpu::LaunchContext& ctx) override {
        ++abs_step_;
        if (ctx.runtime.mpi_rank != 0)
            return;
        if (!xyz_ || !xyz_->is_open())
            return;
        if (abs_step_ != 1 && abs_step_ % traj_every_ != 0)
            return;
        const auto& f = ctx.frame;
        cpu::trajectory_append(*xyz_, f.atoms.type.data(), f.atoms.x.data(),
                               f.atoms.y.data(), f.atoms.z.data(),
                               f.atoms.size(), abs_step_);
    }

  private:
    std::ofstream* xyz_;
    int traj_every_;
    int abs_step_ = 0;
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

// Determine which stages a given run mode actually executes.
bool stage_active(Mode m, int stage_id, int single_id) {
    switch (m) {
    case Mode::full:
        return stage_id >= 1 && stage_id <= 5;
    case Mode::nvt:
        return stage_id == 1;
    case Mode::npt:
        return stage_id == 1 || stage_id == 2;
    case Mode::tension:
        return stage_id == 1 || stage_id == 2 || stage_id == 5;
    case Mode::single_stage:
        return stage_id == single_id;
    }
    return false;
}

int scaled_steps(int base, double scale) {
    if (base <= 0)
        return 0;
    int s = static_cast<int>(base * scale);
    return std::max(s, 1);
}

/// @brief Parse a `--rdf` spec "RMAX,BINS[,STRIDE]" into RDF parameters.
///
/// STRIDE omitted (2 fields) yields freq 0 = sample every step. Returns
/// false on a malformed spec (wrong field count, non-positive RMAX/BINS,
/// negative STRIDE).
bool parse_rdf_spec(const std::string& spec, double& r_max, int& bins,
                    int& freq) {
    std::vector<std::string> parts;
    std::string cur;
    for (char c : spec) {
        if (c == ',') {
            parts.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    parts.push_back(cur);
    if (parts.size() != 2 && parts.size() != 3)
        return false;
    r_max = std::atof(parts[0].c_str());
    bins = std::atoi(parts[1].c_str());
    freq = (parts.size() == 3) ? std::atoi(parts[2].c_str()) : 0;
    return r_max > 0.0 && bins > 0 && freq >= 0;
}

/// @brief Shared evaluators + recorders threaded through every stage.
///
/// Held as lvalues for the whole `main()` so each per-stage `cpu::md::MD`
/// driver can store non-owning pointers to them. Re-registering the same
/// lvalue into a fresh driver is the documented carry-over pattern: the
/// `cpu::Frame` and these stateful recorders persist across the per-stage
/// driver rebuilds.
struct SharedComponents {
    cpu::md::NeighborList<> nl;
    cpu::md::LJ lj;
    HarmonicBond hb;
    HarmonicAngle ha;
    MultiHarmonicDihedral hd;
    CsvRecorder* rec;
    TrajectoryRecorder* trec;
    RDF* rdf; ///< nullptr when --rdf is off
    bool traj_on;
};

/// @brief Register the force-field evaluators, neighbor list and shared
/// recorders that every stage needs into a freshly-constructed driver.
void wire_shared(cpu::md::MD& sim, SharedComponents& sc) {
    sim.set_neighbor(sc.nl);
    sim.set_pair(sc.lj);
    sim.set_bond(sc.hb);
    sim.set_angle(sc.ha);
    sim.set_dihedral(sc.hd);
    sim.on(ATV::md::Phase::end_of_step).add(*sc.rec, 20);
    if (sc.traj_on) {
        sim.on(ATV::md::Phase::end_of_step).add(*sc.trec, 25);
    }
    if (sc.rdf) {
        sim.add(*sc.rdf);
    }
}

#ifdef ATV_BUILD_CUDA

/// @brief Snapshot POD that bundles every per-step scalar the CUDA CSV
///        recorder needs on the host. Trivially copyable so a single
///        async D2H copy can stage it device-to-host. Composed of the
///        canonical `SystemScalars` + the 3-element box diagonal that
///        the 17-column thermo schema needs (`box_x`, `box_y`, `box_z`).
struct SystemScalarsSnapshot {
    SystemScalars scalars;
    cuda::prec_t box_x = 0;
    cuda::prec_t box_y = 0;
    cuda::prec_t box_z = 0;
};

/// @brief CUDA-side recorder mirroring the 17-column extended-schema CPU
///        `CsvRecorder`. Inherits the zero-block sidecar contract from
///        `cuda::AsyncRecordBase<SystemScalarsSnapshot>`: the GPU
///        compute stream stages the snapshot into a device-resident
///        buffer (D2D); the base class enqueues a single
///        non-blocking D2H async copy on the defer stream and submits
///        the CSV-row write to a host worker thread. Hot-path
///        `cudaStreamSynchronize` count stays at zero
///        (msu-pe-cuda-components-07).
///
/// Set the environment variable `ATV_TEST_RECORD_SYNC=1` to force the
/// legacy synchronous D2H path for the chain link 07 correctness test
/// (`tests/cuda/async/test_async_record_correctness.cu`). Default is
/// the async path.
class CudaCsvRecorder
    : public cuda::AsyncRecordBase<SystemScalarsSnapshot> {
  public:
    CudaCsvRecorder(std::ofstream* thermo, std::ofstream* stress,
                    int print_every, int n_global_dof, double box_x_initial)
        : thermo_(thermo), stress_(stress),
          print_every_(std::max(1, print_every)), n_dof_(n_global_dof),
          box_x0_(box_x_initial) {
        if (const char* env = std::getenv("ATV_TEST_RECORD_SYNC"))
            sync_mode_ = (env[0] == '1');
    }

    std::string name() const override { return "cuda_csv_recorder"; }
    int frequency() const override { return print_every_; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::end_of_step), 20}};
    }

    /// @brief Enter a new stage. Mirrors `CsvRecorder::begin_stage`.
    void begin_stage(int stage_id, double T_target, bool record_strain) {
        current_stage_ = stage_id;
        current_T_target_ = T_target;
        record_strain_ = record_strain;
        current_step_offset_ += stage_steps_so_far_ * print_every_;
        stage_steps_so_far_ = 0;
        if (record_strain)
            strain_ref_set_ = false;
    }

    void launch(const cuda::LaunchContext& ctx) override {
        ++stage_steps_so_far_;
        const int abs_step = current_step_offset_ + ctx.step;

        if (sync_mode_) {
            // Legacy synchronous D2H path -- preserved behind a runtime
            // env-flag so chain link 07's correctness test can compare
            // async vs sync byte-by-byte. Production runs (env unset)
            // always take the async branch.
            SystemScalarsSnapshot snap{};
            cudaMemcpy(&snap.scalars, ctx.frame.scalars, sizeof(SystemScalars),
                       cudaMemcpyDeviceToHost);
            cuda::prec_t box_h[9];
            cudaMemcpy(box_h, ctx.frame.box, 9 * sizeof(cuda::prec_t),
                       cudaMemcpyDeviceToHost);
            snap.box_x = box_h[0];
            snap.box_y = box_h[4];
            snap.box_z = box_h[8];
            write_csv_row(snap, abs_step);
            return;
        }

        // Async path. Stage the snapshot on the compute stream so the
        // D2D copies are ordered after any kernel that wrote
        // ctx.frame.scalars / ctx.frame.box in this step (e.g.
        // KineticEnergy, VirialStress, NoseHooverNPT integrator).
        // `stage_device_field` wraps the async D2D copy so the recorder
        // never names the CUDA runtime API directly.
        auto compute_stream =
            static_cast<cudaStream_t>(ctx.runtime.stream);
        stage_device_field(&device_snapshot_->scalars, ctx.frame.scalars,
                           sizeof(SystemScalars), compute_stream);
        // Box diagonal -- 3 separate prec_t fields; the strided source
        // (box[0], box[4], box[8]) makes a single cudaMemcpy2D less
        // clear than three explicit element copies.
        stage_device_field(&device_snapshot_->box_x, &ctx.frame.box[0],
                           sizeof(cuda::prec_t), compute_stream);
        stage_device_field(&device_snapshot_->box_y, &ctx.frame.box[4],
                           sizeof(cuda::prec_t), compute_stream);
        stage_device_field(&device_snapshot_->box_z, &ctx.frame.box[8],
                           sizeof(cuda::prec_t), compute_stream);
        submit_async(static_cast<std::int64_t>(abs_step),
                     [this](const SystemScalarsSnapshot& snap,
                            std::int64_t step) {
                         write_csv_row(snap, static_cast<int>(step));
                     });
    }

  private:
    /// @brief Format one CSV row from a host-resident snapshot. Runs on
    ///        the AsyncRecordBase worker thread (or inline on the hot
    ///        thread when sync_mode_ is set).
    void write_csv_row(const SystemScalarsSnapshot& snap, int abs_step) {
        const double box_x = static_cast<double>(snap.box_x);
        const double box_y = static_cast<double>(snap.box_y);
        const double box_z = static_cast<double>(snap.box_z);
        const auto& s = snap.scalars;
        const double T = (n_dof_ > 0)
                             ? 2.0 * s.kinetic
                                   / (static_cast<double>(n_dof_) * kBoltzmann)
                             : 0.0;
        const double V = box_x * box_y * box_z;
        double P_bar = 0, Pxx = 0, Pyy = 0, Pzz = 0;
        if (V > 0.0) {
            const double trace = s.stress[0] + s.stress[1] + s.stress[2]
                                 + s.stress_kinetic[0] + s.stress_kinetic[1]
                                 + s.stress_kinetic[2];
            P_bar = trace / (3.0 * V) * kEvPerA3ToBar;
            Pxx = (s.stress[0] + s.stress_kinetic[0]) / V * kEvPerA3ToBar;
            Pyy = (s.stress[1] + s.stress_kinetic[1]) / V * kEvPerA3ToBar;
            Pzz = (s.stress[2] + s.stress_kinetic[2]) / V * kEvPerA3ToBar;
        }
        std::lock_guard<std::mutex> lk(io_mtx_);
        if (thermo_ && thermo_->is_open()) {
            (*thermo_) << abs_step << "," << current_stage_ << ","
                       << current_T_target_ << "," << T << "," << P_bar << ","
                       << Pxx << "," << Pyy << "," << Pzz << "," << s.energy
                       << "," << s.energy_pair << "," << s.energy_bond << ","
                       << s.energy_angle << "," << s.energy_dihedral << ","
                       << s.kinetic << "," << box_x << "," << box_y << ","
                       << box_z << "\n";
        }
        if (record_strain_ && stress_ && stress_->is_open()) {
            if (!strain_ref_set_) {
                box_x0_ = box_x;
                strain_ref_set_ = true;
            }
            const double sigma_xx =
                (V > 0.0) ? -(s.stress[0] + s.stress_kinetic[0]) / V : 0.0;
            const double strain = box_x / box_x0_ - 1.0;
            (*stress_) << abs_step << "," << current_stage_ << "," << strain
                       << "," << sigma_xx << "," << (sigma_xx * kEvPerA3ToGPa)
                       << "\n";
        }
    }

    std::ofstream* thermo_;
    std::ofstream* stress_;
    int print_every_;
    int n_dof_;
    double box_x0_;
    bool sync_mode_ = false; ///< ATV_TEST_RECORD_SYNC=1 → legacy sync D2H.
    int current_stage_ = 0;
    double current_T_target_ = 0.0;
    bool record_strain_ = false;
    bool strain_ref_set_ = false;
    int current_step_offset_ = 0;
    int stage_steps_so_far_ = 0;
    std::mutex io_mtx_; ///< Guards thermo_ / stress_ streams across worker + main.
};

/// @brief Run the 5-stage MSU PE protocol on CUDA. Mirrors the CPU
///        pipeline line for line via the categorical setter API on
///        `cuda::md::MD`. `cpu::Frame` ferries state across stages —
///        each per-stage driver materializes to device on `run()` and
///        syncs back when it returns, exactly the pattern used in
///        `tests/cuda/md/test_md.cu`.
void run_cuda_pipeline(const Args& a, cpu::Frame& frame,
                       std::ofstream& thermo_csv, std::ofstream& stress_csv,
                       int n_global_dof, double box_x0, bool quiet) {
    CudaCsvRecorder rec{&thermo_csv, &stress_csv, a.print_every, n_global_dof,
                        box_x0};

    // cuda::DihedralMultiHarmonic takes std::vector<std::array<double, 5>>
    // (one per type); the CPU helper returns std::vector<double>. Convert
    // once for both stage-1 and the NPT stages.
    auto Cn_vec = dihedral_Cn_eV();
    std::vector<std::array<double, 5>> cuda_Cn{
        {Cn_vec[0], Cn_vec[1], Cn_vec[2], Cn_vec[3], Cn_vec[4]}};

    // Per-stage timing log written to <out-dir>/stage_timing.csv so we
    // can compare wall time / ns-per-day with the matching LAMMPS log
    // breakdown. One line per executed stage; columns mirror LAMMPS's
    // "Performance" line.
    std::ofstream timing_csv(a.out_dir + "/stage_timing.csv");
    timing_csv << "stage,steps,dt_fs,wall_seconds,steps_per_second,ns_per_day\n";
    auto log_stage = [&](int stage_id, int steps, double dt_fs, double wall_s) {
        const double sps = (wall_s > 0.0) ? steps / wall_s : 0.0;
        const double ns_per_day =
            (wall_s > 0.0)
                ? (double(steps) * dt_fs * 1e-6) * (86400.0 / wall_s)
                : 0.0;
        timing_csv << stage_id << "," << steps << "," << dt_fs << ","
                   << wall_s << "," << sps << "," << ns_per_day << "\n";
        timing_csv.flush();
        if (!quiet)
            std::fprintf(stderr,
                         "  stage %d wall=%.2fs  %.0f step/s  %.2f ns/day\n",
                         stage_id, wall_s, sps, ns_per_day);
    };

    // Stage 1 — NVELimit + FixLangevin @ T_hot.
    if (stage_active(a.mode, 1, a.single_stage_id)) {
        const int steps = scaled_steps(a.stage1_steps, a.scale);
        if (!quiet)
            std::fprintf(stderr,
                         "  stage 1 (cuda): NVE/limit + Langevin @ %.1f K, "
                         "dt=%.2f fs, %d steps\n",
                         a.T_hot, a.dt_stage1, steps);
        rec.begin_stage(1, a.T_hot, /*record_strain=*/false);
        cuda::md::MD sim;
        sim.set_neighbor(cuda::NeighborList{kLjCutoff});
        if (!a.disable_lj)
            sim.set_pair(cuda::PairLJCut{lj_eps_eV(), kLjSigma, kLjCutoff}
                             .set_special_bonds(0.0, 0.0, 0.0));
        if (!a.disable_bond)
            sim.set_bond(cuda::BondHarmonic{{bond_k_eV()}, {kBondR0}});
        if (!a.disable_angle)
            sim.set_angle(cuda::AngleHarmonic{{angle_k_eV()}, {kAngleTheta0}});
        if (!a.disable_dihedral)
            sim.set_dihedral(cuda::DihedralMultiHarmonic{cuda_Cn});
        sim.set_integrator(cuda::md::NVELimit{a.dt_stage1, 0.05});
        if (!std::getenv("ATV_NO_INITVEL"))
            sim.add(cuda::md::InitVelocity{a.T_hot,
                                       static_cast<unsigned long long>(a.seed)});
        if (!a.disable_langevin)
            sim.set_thermostat(cuda::md::FixLangevin{
                a.T_hot, a.T_hot, a.langevin_damp, a.dt_stage1,
                static_cast<std::uint32_t>(a.seed + 1)});
        sim.add(rec);  // uses default_bindings: end_of_step @ +20
        sim.build();
        const double ts0 = now_seconds();
        sim.run(frame, steps);
        log_stage(1, steps, a.dt_stage1, now_seconds() - ts0);
        // DIAG: dump per-atom forces right after stage 1 step `steps`.
        // For step-0 comparison with LAMMPS, run with `--scale 0.0001`
        // (1 step). The frame is now host-synced (DynamicsDriver retrieves
        // on every run()), so the host-side fx/fy/fz hold device forces
        // from the LAST evaluated step.
        if (std::getenv("ATV_DUMP_FORCES")) {
            FILE* fp = std::fopen("/tmp/cuda_forces.dump", "w");
            std::fprintf(fp, "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n%d\n"
                             "ITEM: BOX BOUNDS pp pp pp\n0 1\n0 1\n0 1\n"
                             "ITEM: ATOMS id fx fy fz\n",
                         frame.atoms.size());
            for (int i = 0; i < frame.atoms.size(); ++i)
                std::fprintf(fp, "%lld %.10e %.10e %.10e\n",
                             (long long)(frame.atoms.tag[i] + 1),
                             frame.atoms.fx[i], frame.atoms.fy[i],
                             frame.atoms.fz[i]);
            std::fclose(fp);
        }
    }

    // Helper to build a per-stage NPT driver with the shared FF.
    auto wire_npt_stage = [&](cuda::md::MD& sim,
                              cuda::md::NoseHooverNPT& npt) {
        sim.set_neighbor(cuda::NeighborList{kLjCutoff});
        if (!a.disable_lj)
            sim.set_pair(cuda::PairLJCut{lj_eps_eV(), kLjSigma, kLjCutoff}
                             .set_special_bonds(0.0, 0.0, 0.0));
        if (!a.disable_bond)
            sim.set_bond(cuda::BondHarmonic{{bond_k_eV()}, {kBondR0}});
        if (!a.disable_angle)
            sim.set_angle(cuda::AngleHarmonic{{angle_k_eV()}, {kAngleTheta0}});
        if (!a.disable_dihedral)
            sim.set_dihedral(cuda::DihedralMultiHarmonic{cuda_Cn});
        sim.set_integrator(npt);
        sim.add(rec);  // uses default_bindings: end_of_step @ +20
    };

    // Stage 2 — iso NoseHooverNPT @ T_hot.
    if (stage_active(a.mode, 2, a.single_stage_id)) {
        const int steps = scaled_steps(a.stage2_steps, a.scale);
        if (!quiet)
            std::fprintf(stderr,
                         "  stage 2 (cuda): NPT (iso) @ %.1f K, dt=%.2f fs, "
                         "%d steps\n",
                         a.T_hot, a.dt_stage2, steps);
        rec.begin_stage(2, a.T_hot, /*record_strain=*/false);
        cuda::md::NoseHooverNPT::Config cfg{
            .T_target = a.T_hot, .tau_T = a.tau_T, .P_target = 0.0,
            .tau_P = a.tau_P, .dt = a.dt_stage2, .drag = 2.0};  // LAMMPS `drag 2`
        cuda::md::MD sim;
        cuda::md::NoseHooverNPT npt{cfg};
        wire_npt_stage(sim, npt);
        sim.add(cuda::md::ZeroMomentum{1});
        sim.build();
        const double ts0 = now_seconds();
        sim.run(frame, steps);
        log_stage(2, steps, a.dt_stage2, now_seconds() - ts0);
    }

    // Stage 3 — ramped NoseHooverNPT T_hot → T_cold.
    if (stage_active(a.mode, 3, a.single_stage_id)) {
        const int steps = scaled_steps(a.stage3_steps, a.scale);
        if (!quiet)
            std::fprintf(stderr,
                         "  stage 3 (cuda): NPT cool %.1f -> %.1f K (ramp), "
                         "dt=%.2f fs, %d steps\n",
                         a.T_hot, a.T_cold, a.dt_stage3, steps);
        rec.begin_stage(3, a.T_cold, /*record_strain=*/false);
        cuda::md::NoseHooverNPT::Config cfg{
            .T_target = a.T_hot, .tau_T = a.tau_T, .P_target = 0.0,
            .tau_P = a.tau_P, .dt = a.dt_stage3, .drag = 2.0,
            .T_stop = a.T_cold};
        cuda::md::MD sim;
        cuda::md::NoseHooverNPT npt{cfg};
        npt.set_run_span(0, std::max(0, steps - 1));
        wire_npt_stage(sim, npt);
        sim.add(cuda::md::ZeroMomentum{1});
        sim.build();
        const double ts0 = now_seconds();
        sim.run(frame, steps);
        log_stage(3, steps, a.dt_stage3, now_seconds() - ts0);
    }

    // Stage 4 — iso NoseHooverNPT @ T_cold.
    if (stage_active(a.mode, 4, a.single_stage_id)) {
        const int steps = scaled_steps(a.stage4_steps, a.scale);
        if (!quiet)
            std::fprintf(stderr,
                         "  stage 4 (cuda): NPT (iso) @ %.1f K, dt=%.2f fs, "
                         "%d steps\n",
                         a.T_cold, a.dt_stage4, steps);
        rec.begin_stage(4, a.T_cold, /*record_strain=*/false);
        cuda::md::NoseHooverNPT::Config cfg{
            .T_target = a.T_cold, .tau_T = a.tau_T, .P_target = 0.0,
            .tau_P = a.tau_P, .dt = a.dt_stage4, .drag = 2.0};
        cuda::md::MD sim;
        cuda::md::NoseHooverNPT npt{cfg};
        wire_npt_stage(sim, npt);
        sim.add(cuda::md::ZeroMomentum{1});
        sim.build();
        const double ts0 = now_seconds();
        sim.run(frame, steps);
        log_stage(4, steps, a.dt_stage4, now_seconds() - ts0);
    }

    // Stage 5 — aniso NoseHooverNPT (y/z) + Deform x.
    if (stage_active(a.mode, 5, a.single_stage_id)) {
        const int steps = scaled_steps(a.stage5_steps, a.scale);
        if (!quiet)
            std::fprintf(stderr,
                         "  stage 5 (cuda): tension epsdot=%.2e /fs @ %.1f K, "
                         "dt=%.2f fs, %d steps\n",
                         a.strain_rate, a.T_cold, a.dt_stage5, steps);
        rec.begin_stage(5, a.T_cold, /*record_strain=*/true);
        cuda::md::NoseHooverNPT::Config cfg{
            .T_target = a.T_cold, .tau_T = a.tau_T, .P_target = 0.0,
            .tau_P = a.tau_P, .dt = a.dt_stage5, .drag = 2.0,
            .couple_axis = {false, true, true}};
        cuda::md::MD sim;
        cuda::md::NoseHooverNPT npt{cfg};
        wire_npt_stage(sim, npt);
        sim.add(cuda::md::Deform{/*axis=*/0, a.strain_rate, a.dt_stage5});
        sim.build();
        const double ts0 = now_seconds();
        sim.run(frame, steps);
        log_stage(5, steps, a.dt_stage5, now_seconds() - ts0);
    }
}

#endif // ATV_BUILD_CUDA

} // namespace

int main(int argc, char** argv) {
    Args a;
    if (!parse_args(argc, argv, a))
        return 2;

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

    // The cpu::Frame is the carry-over state — it persists across every
    // per-stage driver rebuild, holding positions, velocities and the box.
    cpu::Frame frame = a.data_path.empty()
                           ? make_pe_box(a.chains_per_axis, a.n_atoms_per_chain,
                                         a.box_len, a.margin)
                           : cpu::load_lammps_data(a.data_path);
    const int n_atoms = frame.atoms.size();
    const int n_global_dof = 3 * n_atoms - 3;
    const double box_x0 = frame.box(0, 0);

    // `traj_xyz` is opened only when --traj-every is set, so a default
    // run leaves CSV-only output byte-identical to the pre-flag behaviour.
    std::ofstream thermo_csv, stress_csv, traj_xyz;
    if (my_rank == 0) {
        thermo_csv.open(a.out_dir + "/thermo.csv");
        stress_csv.open(a.out_dir + "/stress_strain.csv");
        if (thermo_csv.is_open()) {
            // Phase G.2 + G.3 17-column extended schema
            // (msu-pe-cuda-components-01). Pxx/Pyy/Pzz are real per-axis
            // virials; 4 per-term energies decompose E_pot.
            thermo_csv
                << "step,stage,T_target_K,T_K,P_bar,Pxx,Pyy,Pzz,E_pot_eV,"
                   "E_pair_eV,E_bond_eV,E_angle_eV,E_dihedral_eV,E_kin_eV,"
                   "box_x,box_y,box_z\n";
        }
        if (stress_csv.is_open()) {
            stress_csv << "step,stage,strain,sigma_xx_eV_per_A3,sigma_xx_GPa\n";
        }
        if (a.traj_every > 0) {
            traj_xyz.open(a.out_dir + "/trajectory.xyz");
        }
    }

    // ── Shared evaluators + recorders (lvalues; outlive every driver) ──
    cpu::md::NeighborList<> nl{kLjCutoff};
    cpu::md::LJ lj =
        cpu::md::LJ{lj_eps_eV(), kLjSigma, kLjCutoff}.set_special_bonds(
            0.0, 0.0, 0.0);
    HarmonicBond hb{{bond_k_eV()}, {kBondR0}};
    HarmonicAngle ha{{angle_k_eV()}, {kAngleTheta0}};
    MultiHarmonicDihedral hd{{dihedral_Cn_eV()}};
    CsvRecorder rec{&thermo_csv, &stress_csv, a.print_every, n_global_dof,
                    box_x0};
    TrajectoryRecorder trec{&traj_xyz, std::max(1, a.traj_every)};

    // RDF observer (optional). Heap-allocated so it is conditionally
    // created yet still outlives every per-stage driver, which stores a
    // non-owning pointer; finalize() writes <out-dir>/rdf.dat at the end
    // of each stage's driver teardown.
    std::unique_ptr<RDF> rdf;
    if (!a.rdf_spec.empty()) {
        double r_max = 0.0;
        int bins = 0, freq = 0;
        if (!parse_rdf_spec(a.rdf_spec, r_max, bins, freq)) {
            std::fprintf(
                stderr,
                "invalid --rdf spec '%s' (expected RMAX,BINS[,STRIDE])\n",
                a.rdf_spec.c_str());
            return 2;
        }
        rdf = std::make_unique<RDF>(r_max, bins, freq, a.out_dir + "/rdf.dat");
    }

    SharedComponents shared{
        nl, lj, hb, ha, hd, &rec, &trec, rdf.get(), a.traj_every > 0};

    if (my_rank == 0 && !a.quiet) {
        const char* mode_str = a.mode == Mode::full      ? "full"
                               : a.mode == Mode::nvt     ? "NVT"
                               : a.mode == Mode::npt     ? "NPT"
                               : a.mode == Mode::tension ? "tension"
                                                         : "single_stage";
        std::fprintf(
            stderr,
            "polyethylene_tension: mode=%s atoms=%d T_hot=%.1f T_cold=%.1f\n"
            "  per-stage dt[fs] = {%.2f, %.2f, %.2f, %.2f, %.2f}  (mirrors "
            "lammps.in)\n",
            mode_str, n_atoms, a.T_hot, a.T_cold, a.dt_stage1, a.dt_stage2,
            a.dt_stage3, a.dt_stage4, a.dt_stage5);
    }

    const double t0 = now_seconds();

    // ── Dispatch: CUDA pipeline (chain link 06) or CPU stages below ──
    if (a.backend == Backend::cuda) {
#ifdef ATV_BUILD_CUDA
        if (my_rank == 0) {
            try {
                run_cuda_pipeline(a, frame, thermo_csv, stress_csv,
                                  n_global_dof, box_x0, a.quiet);
            } catch (const std::exception& e) {
                std::fprintf(stderr,
                             "polyethylene_tension: CUDA pipeline threw: %s\n",
                             e.what());
                return 3;
            } catch (...) {
                std::fprintf(stderr,
                             "polyethylene_tension: CUDA pipeline threw "
                             "unknown exception\n");
                return 3;
            }
        }
#else
        std::fprintf(stderr,
                     "polyethylene_tension: --backend cuda not compiled in. "
                     "Rebuild with -DATV_BUILD_CUDA=ON.\n");
        return 2;
#endif
        // Close output streams + finalize MPI before returning.
        if (my_rank == 0) {
            if (thermo_csv.is_open())
                thermo_csv.close();
            if (stress_csv.is_open())
                stress_csv.close();
            if (traj_xyz.is_open())
                traj_xyz.close();
        }
#ifdef ATV_HAS_MPI
        MPI_Finalize();
#endif
        return 0;
    }

    // ── Per-stage DD handoff (multi-rank only). ─────────────────────────
    // Each per-stage `cpu::md::MD` instance owns its own Backend ⇒ its
    // own `DomainDecomposition`; a fresh DD has `n_local==0` so
    // `partition_atoms`'s `first_partition` branch fires and re-scans the
    // *whole* frame allocation, claiming stale tail rows from the
    // previous stage's compaction as owned atoms (LJ then finds duplicate
    // positions ⇒ step-0 blow-up). Threading the same DD through every
    // stage preserves `n_local` so subsequent `partition_atoms` iterates
    // only [0..n_local), which is the well-tested "Driver::run() twice"
    // path. cart_comm is recreated on each stage's `init_topology()`;
    // that is intentional (the comm shouldn't outlive its parent
    // Backend) and fast (microseconds).
    std::unique_ptr<cpu::md::mpi::DomainDecomposition> preserved_dd;

    // ── Stage 1 — NVELimit + FixLangevin @ T_hot, dt 1.0 fs ───────────
    // Mirrors lammps.in `fix nve/limit 0.05` + `fix langevin 500 500 10`
    // + `velocity create 500`. NVELimit's per-step displacement cap
    // tames the unrelaxed-config atom overlap (E_pair ~3.4e8 at step 0)
    // that plain Velocity-Verlet cannot survive.
    if (stage_active(a.mode, 1, a.single_stage_id)) {
        const int steps = scaled_steps(a.stage1_steps, a.scale);
        if (my_rank == 0 && !a.quiet)
            std::fprintf(stderr,
                         "  stage 1: NVE/limit + Langevin @ %.1f K, dt=%.2f "
                         "fs, %d steps\n",
                         a.T_hot, a.dt_stage1, steps);
        cpu::md::NVELimit nve{a.dt_stage1, /*xmax=*/0.05};
        cpu::md::InitVelocity init_v{a.T_hot,
                                     static_cast<unsigned int>(a.seed)};
        cpu::md::FixLangevin langevin{a.T_hot, /*T_stop=*/a.T_hot,
                                      a.langevin_damp, a.dt_stage1,
                                      static_cast<std::uint32_t>(a.seed + 1)};
        cpu::md::MD sim;
        sim.set_integrator(nve);
        sim.add(init_v);
        sim.set_thermostat(langevin);
        wire_shared(sim, shared);
        rec.begin_stage(1, a.T_hot, /*record_strain=*/false);
        if (preserved_dd)
            sim.adopt_dd(std::move(preserved_dd));
        sim.build();
        sim.run(frame, steps);
        preserved_dd = sim.release_dd();
    }

    // ── Stage 2 — isotropic NoseHooverNPT @ T_hot, P=0, dt 0.5 fs ─────
    // Mirrors lammps.in `fix npt temp 500 500 50 iso 0 0 1000 drag 2`.
    if (stage_active(a.mode, 2, a.single_stage_id)) {
        const int steps = scaled_steps(a.stage2_steps, a.scale);
        if (my_rank == 0 && !a.quiet)
            std::fprintf(
                stderr, "  stage 2: NPT (iso) @ %.1f K, dt=%.2f fs, %d steps\n",
                a.T_hot, a.dt_stage2, steps);
        cpu::md::NoseHooverNPT npt{{.T_target = a.T_hot,
                                    .tau_T = a.tau_T,
                                    .P_target = 0.0,
                                    .tau_P = a.tau_P,
                                    .dt = a.dt_stage2,
                                    .drag = 2.0}};
        // LAMMPS `fix momentum 1 linear 1 1 1`: zero net COM drift each
        // step under NPT + Langevin coupling (Phase G.1). Per-stage stack
        // lvalue — modulo counter starts fresh at each stage; freq=1 so
        // the counter doesn't matter for correctness.
        cpu::md::ZeroMomentum zm{1};
        cpu::md::MD sim;
        sim.set_integrator(npt);
        sim.add(zm);
        wire_shared(sim, shared);
        rec.begin_stage(2, a.T_hot, /*record_strain=*/false);
        if (preserved_dd)
            sim.adopt_dd(std::move(preserved_dd));
        sim.build();
        sim.run(frame, steps);
        preserved_dd = sim.release_dd();
    }

    // ── Stage 3 — ramped NoseHooverNPT 500 → 100 K, P=0, dt 0.5 fs ────
    // Mirrors lammps.in `fix npt temp 500 100 50 iso 0 0 1000 drag 2`:
    // a SINGLE NoseHooverNPT with T_start=T_hot, T_stop=T_cold and a
    // run span covering the whole stage. The integrator interpolates
    // the setpoint internally — no cooling micro-segment loop.
    if (stage_active(a.mode, 3, a.single_stage_id)) {
        const int steps = scaled_steps(a.stage3_steps, a.scale);
        if (my_rank == 0 && !a.quiet)
            std::fprintf(
                stderr,
                "  stage 3: NPT cool %.1f -> %.1f K (continuous ramp), "
                "dt=%.2f fs, %d steps\n",
                a.T_hot, a.T_cold, a.dt_stage3, steps);
        cpu::md::NoseHooverNPT npt{{.T_target = a.T_hot,
                                    .tau_T = a.tau_T,
                                    .P_target = 0.0,
                                    .tau_P = a.tau_P,
                                    .dt = a.dt_stage3,
                                    .drag = 2.0,
                                    .T_stop = a.T_cold}};
        // Run span [0, steps-1]: a fresh per-stage driver restarts its
        // step counter at 0, so the ramp delta sweeps 0 -> 1 over the
        // stage exactly as LAMMPS `update->beginstep / endstep`.
        npt.set_run_span(0, std::max(0, steps - 1));
        // LAMMPS `fix momentum 1 linear 1 1 1` (Phase G.1).
        cpu::md::ZeroMomentum zm{1};
        cpu::md::MD sim;
        sim.set_integrator(npt);
        sim.add(zm);
        wire_shared(sim, shared);
        rec.begin_stage(3, a.T_cold, /*record_strain=*/false);
        if (preserved_dd)
            sim.adopt_dd(std::move(preserved_dd));
        sim.build();
        sim.run(frame, steps);
        preserved_dd = sim.release_dd();
    }

    // ── Stage 4 — isotropic NoseHooverNPT @ T_cold, P=0, dt 0.5 fs ────
    // Mirrors lammps.in `fix npt temp 100 100 50 iso 0 0 1000 drag 2`.
    if (stage_active(a.mode, 4, a.single_stage_id)) {
        const int steps = scaled_steps(a.stage4_steps, a.scale);
        if (my_rank == 0 && !a.quiet)
            std::fprintf(
                stderr, "  stage 4: NPT (iso) @ %.1f K, dt=%.2f fs, %d steps\n",
                a.T_cold, a.dt_stage4, steps);
        cpu::md::NoseHooverNPT npt{{.T_target = a.T_cold,
                                    .tau_T = a.tau_T,
                                    .P_target = 0.0,
                                    .tau_P = a.tau_P,
                                    .dt = a.dt_stage4,
                                    .drag = 2.0}};
        // LAMMPS `fix momentum 1 linear 1 1 1` (Phase G.1).
        cpu::md::ZeroMomentum zm{1};
        cpu::md::MD sim;
        sim.set_integrator(npt);
        sim.add(zm);
        wire_shared(sim, shared);
        rec.begin_stage(4, a.T_cold, /*record_strain=*/false);
        if (preserved_dd)
            sim.adopt_dd(std::move(preserved_dd));
        sim.build();
        sim.run(frame, steps);
        preserved_dd = sim.release_dd();
    }

    // ── Stage 5 — anisotropic NoseHooverNPT + Deform x, dt 1.0 fs ─────
    // Mirrors lammps.in `fix npt temp 100 100 50 y 0 0 1000 z 0 0 1000`
    // + `fix deform 1 x erate 1e-5`: x held rigid by the NPT barostat
    // (couple_axis = {false, true, true}) while Deform drives the x box
    // edge at the engineering strain rate; y and z barostat to P=0.
    if (stage_active(a.mode, 5, a.single_stage_id)) {
        const int steps = scaled_steps(a.stage5_steps, a.scale);
        if (my_rank == 0 && !a.quiet)
            std::fprintf(
                stderr,
                "  stage 5: tension epsdot=%.2e /fs @ %.1f K, dt=%.2f fs, "
                "%d steps\n",
                a.strain_rate, a.T_cold, a.dt_stage5, steps);
        cpu::md::NoseHooverNPT npt{{.T_target = a.T_cold,
                                    .tau_T = a.tau_T,
                                    .P_target = 0.0,
                                    .tau_P = a.tau_P,
                                    .dt = a.dt_stage5,
                                    .drag = 2.0,
                                    .couple_axis = {false, true, true}}};
        cpu::md::Deform deform{cpu::md::Deform::Mode::erate, a.strain_rate,
                               a.dt_stage5, /*axis=*/0};
        cpu::md::MD sim;
        sim.set_integrator(npt);
        sim.add(deform);
        wire_shared(sim, shared);
        rec.begin_stage(5, a.T_cold, /*record_strain=*/true);
        if (preserved_dd)
            sim.adopt_dd(std::move(preserved_dd));
        sim.build();
        sim.run(frame, steps);
        preserved_dd = sim.release_dd();
    }

    const double t1 = now_seconds();

    if (my_rank == 0 && !a.quiet) {
        std::fprintf(stderr,
                     "polyethylene_tension done: %d atoms, %.2fs total\n",
                     n_atoms, t1 - t0);
        if (a.traj_every > 0)
            std::fprintf(stderr, "  trajectory : %s/trajectory.xyz\n",
                         a.out_dir.c_str());
        if (rdf)
            std::fprintf(stderr, "  rdf        : %s/rdf.dat\n",
                         a.out_dir.c_str());
    }

    if (my_rank == 0) {
        if (thermo_csv.is_open())
            thermo_csv.close();
        if (stress_csv.is_open())
            stress_csv.close();
        if (traj_xyz.is_open())
            traj_xyz.close();
    }

#ifdef ATV_HAS_MPI
    MPI_Finalize();
#endif
    return 0;
}
