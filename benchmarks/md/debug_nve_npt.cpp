/// @file benchmarks/md/debug_nve_npt.cpp
/// @brief Minimal single-rank NVE→NPT debug driver on the MSU 10k-atom PE
///        system, for fast step-by-step comparison against LAMMPS.
///
/// Rationale. The full 5-stage polyethylene_tension benchmark is a ~24h,
/// MPI, thermostatted run — too heavy to iterate on when a force/energy bug
/// must be pinned against LAMMPS. This driver loads the same `msu_pe.data`,
/// applies the same Dreiding-UA force field as `lammps.in`, and runs a SHORT
/// relaxation burst followed by a SHORT NPT burst, writing one CSV row per
/// step (energies per term, KE, T, box) so the run can be diffed against a
/// matching LAMMPS deck (`lammps_debug.in`) column-by-column.
///
/// Two comparison modes against LAMMPS:
///
///   1. **Single point (step 0).** With zero initial velocities and zeroed
///      forces, the first `step_once` performs no motion (initial_integrate
///      sees v=0, f=0) and simply evaluates the force field at the loaded
///      geometry. The step-0 CSV row is therefore the single-point energy of
///      `msu_pe.data` and must match LAMMPS `run 0` term-by-term. This is the
///      strongest bug detector and is validated to <0.02% on every term.
///
///   2. **Dynamics.** The raw MSU packing is far from equilibrium (genuine LJ
///      overlaps both codes agree on), so a *pure* NVE burst blows up — MSU
///      stage 1 relaxes with `fix nve/limit`. This driver mirrors that: the
///      relaxation stage uses `nve/limit xmax` (`--xmax`, default 0.1 Å; pass
///      `--xmax 0` for plain `fix nve`). The comparison against LAMMPS here is
///      energy conservation / ensemble behaviour, not a bit-for-bit
///      trajectory (AV takes t=0 as the state after the priming step, the
///      repo convention — cf. bench_lammps_npt_parity).
///
/// Usage:
///   debug_nve_npt [--data msu_pe.data] [--nve N=200] [--npt M=200]
///                 [--dt-nve 1.0] [--dt-npt 0.5] [--xmax 0.1]
///                 [--out debug_av.csv]
///
/// Units: AV is eV / Å / fs internally; energies are emitted BOTH in eV and
/// in kcal/mol (×23.060547830619026) so they line up with LAMMPS `real` units.

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/lammps_data.h"
#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/neighbor_list.h"
#include "atomiverse/cpu/md/nose_hoover_npt.h"
#include "atomiverse/cpu/md/nve_limit.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/md/virial_stress.h"
#include "atomiverse/cpu/potentials/bond_harmonic.h"
#include "atomiverse/cpu/potentials/angle_harmonic.h"
#include "atomiverse/cpu/potentials/dihedral_multi_harmonic.h"

using namespace ATV;

namespace {

// ── Dreiding-UA force field — verbatim from lammps.in / lammps_msu_paper.in ──
constexpr double kKcalToEv = 1.0 / 23.060547830619026;
constexpr double kBoltzEv = 8.617333262e-5; // eV/K
constexpr double kPI = 3.14159265358979323846;

constexpr double kBondK_kcal = 350.0;  // bond_coeff 1 350 1.53
constexpr double kBondR0 = 1.53;
constexpr double kAngleK_kcal = 60.0;   // angle_coeff 1 60 109.5
constexpr double kAngleTheta0 = 109.5 * kPI / 180.0;
constexpr double kLjEps_kcal = 0.112;   // pair_coeff 1 1 0.112 4.01 10.5
constexpr double kLjSigma = 4.01;
constexpr double kLjCut = 10.5;

std::vector<double> dihedral_Cn_eV() {
    return {1.73 * kKcalToEv, -4.49 * kKcalToEv, 0.776 * kKcalToEv,
            6.99 * kKcalToEv, 0.0};
}

struct Args {
    std::string data = "benchmarks/md/reference/polyethylene/msu_pe.data";
    std::string out = "debug_av.csv";
    int nve = 200;
    int npt = 200;
    double dt_nve = 1.0;
    double dt_npt = 0.5;
    double xmax = 0.1; // nve/limit cap [Å]; 0 ⇒ plain fix nve.
};

void wire_forcefield(cpu::md::MD& sim) {
    sim.on(ATV::md::Phase::force).add(cpu::md::NeighborList<>{kLjCut}, -10);
    sim.on(ATV::md::Phase::force)
        .add(cpu::md::LJ{kLjEps_kcal * kKcalToEv, kLjSigma, kLjCut}
                 .set_special_bonds(0.0, 0.0, 0.0))
        .add(cpu::BondHarmonic{{kBondK_kcal * kKcalToEv}, {kBondR0}})
        .add(cpu::AngleHarmonic{{kAngleK_kcal * kKcalToEv}, {kAngleTheta0}})
        .add(cpu::DihedralMultiHarmonic{{dihedral_Cn_eV()}});
    sim.add(cpu::md::VirialStress{});
}

// Per-step thermo row. Energies emitted in eV and kcal/mol; T from KE.
void write_header(std::ofstream& os) {
    os << "step,stage,pe_eV,e_pair_eV,e_bond_eV,e_angle_eV,e_dihedral_eV,"
          "ke_eV,pe_kcal,ke_kcal,temp_K,box_x,box_y,box_z\n";
}

void write_row(std::ofstream& os, int step, const char* stage,
               const cpu::Frame& f, int ndof) {
    const auto& s = f.scalars;
    const double T = (ndof > 0) ? 2.0 * s.kinetic / (ndof * kBoltzEv) : 0.0;
    char buf[512];
    std::snprintf(
        buf, sizeof(buf),
        "%d,%s,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.6f,%.6f,%.6f,"
        "%.6f\n",
        step, stage, s.energy, s.energy_pair, s.energy_bond, s.energy_angle,
        s.energy_dihedral, s.kinetic, s.energy / kKcalToEv,
        s.kinetic / kKcalToEv, T, f.box(0, 0), f.box(1, 1), f.box(2, 2));
    os << buf;
}

bool parse(int argc, char** argv, Args& a) {
    for (int i = 1; i < argc; ++i) {
        auto eq = [&](const char* f) { return std::strcmp(argv[i], f) == 0; };
        if (eq("--data") && i + 1 < argc)
            a.data = argv[++i];
        else if (eq("--out") && i + 1 < argc)
            a.out = argv[++i];
        else if (eq("--nve") && i + 1 < argc)
            a.nve = std::atoi(argv[++i]);
        else if (eq("--npt") && i + 1 < argc)
            a.npt = std::atoi(argv[++i]);
        else if (eq("--dt-nve") && i + 1 < argc)
            a.dt_nve = std::atof(argv[++i]);
        else if (eq("--dt-npt") && i + 1 < argc)
            a.dt_npt = std::atof(argv[++i]);
        else if (eq("--xmax") && i + 1 < argc)
            a.xmax = std::atof(argv[++i]);
        else {
            std::fprintf(stderr, "unknown arg: %s\n", argv[i]);
            return false;
        }
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    Args a;
    if (!parse(argc, argv, a))
        return 2;

    cpu::Frame frame = cpu::load_lammps_data(a.data);
    const int n = frame.atoms.size();
    const int ndof = 3 * n - 3; // LAMMPS default (subtract COM dof)
    // Deterministic start: zero velocities so the NVE trajectory is a pure
    // function of the forces and matches LAMMPS (no RNG-seeded velocities).
    for (int i = 0; i < n; ++i) {
        frame.atoms.vx[i] = 0.0;
        frame.atoms.vy[i] = 0.0;
        frame.atoms.vz[i] = 0.0;
    }
    std::fprintf(stderr, "[debug_nve_npt] loaded %d atoms from %s\n", n,
                 a.data.c_str());

    std::ofstream os(a.out);
    write_header(os);

    // ── Stage 1: relaxation (fix nve/limit, or plain fix nve if xmax==0) ──
    {
        cpu::md::MD sim;
        if (a.xmax > 0.0)
            sim.add(cpu::md::NVELimit{a.dt_nve, a.xmax});
        else
            sim.add(cpu::md::VelocityVerlet{a.dt_nve});
        wire_forcefield(sim);
        sim.build();
        // Priming step: with v=0 and f=0, initial_integrate moves nothing, so
        // this step just evaluates the force field at the loaded geometry and
        // leaves f(x0) in the frame (final_integrate's half-kick is discarded
        // below by re-zeroing v). The result — x=x0, v=0, f=f(x0) — is exactly
        // LAMMPS's post-`setup` state, so the step-0 row is the single point
        // (== LAMMPS `run 0`, KE=0) and step 1 onward matches `fix nve`
        // step-for-step.
        sim.run(frame, 1);
        for (int i = 0; i < n; ++i) {
            frame.atoms.vx[i] = 0.0;
            frame.atoms.vy[i] = 0.0;
            frame.atoms.vz[i] = 0.0;
        }
        frame.scalars.kinetic = 0.0;
        write_row(os, 0, "nve", frame, ndof);
        for (int step = 1; step <= a.nve; ++step) {
            sim.run(frame, 1);
            write_row(os, step, "nve", frame, ndof);
        }
        std::fprintf(stderr,
                     "[debug_nve_npt] relax %d steps done (xmax=%.3g)\n", a.nve,
                     a.xmax);
    }

    // ── Stage 2: NPT (fix npt temp T iso 0) ──────────────────────────────
    {
        cpu::md::NoseHooverNPT::Config cfg;
        cfg.T_target = 300.0;
        cfg.tau_T = 100.0 * a.dt_npt;
        cfg.P_target = 0.0;
        cfg.tau_P = 1000.0 * a.dt_npt;
        cfg.dt = a.dt_npt;
        cpu::md::MD sim;
        sim.add(cpu::md::VelocityVerlet{a.dt_npt});
        sim.add(cpu::md::NoseHooverNPT{cfg});
        wire_forcefield(sim);
        sim.build();
        for (int step = 0; step <= a.npt; ++step) {
            if (step > 0)
                sim.run(frame, 1);
            write_row(os, a.nve + step, "npt", frame, ndof);
        }
        std::fprintf(stderr, "[debug_nve_npt] NPT %d steps done -> %s\n", a.npt,
                     a.out.c_str());
    }
    return 0;
}
