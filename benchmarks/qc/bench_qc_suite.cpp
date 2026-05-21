/// @file benchmarks/qc/bench_qc_suite.cpp
/// @brief Phase 1.2 QC benchmark suite — Minimize + SCFEvaluator at PBE /
///        def2-SVP, Armijo backtracking on. Compared against
///        geomeTRIC-optimized PySCF references.
///
/// Build / run:
///     cmake --build build --target bench_qc_suite
///     ./build/benchmarks/bench_qc_suite           # all molecules with a ref
///     ./build/benchmarks/bench_qc_suite h2o nh3   # specific stems
///
/// Reference generation:
///     python benchmarks/qc/reference/generate_references.py \
///       --method rks --xc pbe --basis def2-svp
///
/// Data-on-disk contract:
///     benchmarks/qc/molecules/<stem>.xyz          - starting geometry
///     benchmarks/qc/reference/<stem>.json         - PySCF reference (geom + E_opt)
///
/// Paths resolve relative to the binary's cwd unless overridden by the
/// build system via `ATV_QC_DATA_DIR` (set by CMakeLists.txt).
#ifdef ATV_HAS_RHF

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;

namespace {

namespace fs = std::filesystem;

int z_from_symbol(const std::string& s) {
    static const std::pair<const char*, int> TBL[] = {
        {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5},
        {"C", 6}, {"N", 7},  {"O", 8},  {"F", 9},  {"Ne", 10},
    };
    for (const auto& [sym, z] : TBL)
        if (s == sym)
            return z;
    throw std::runtime_error("z_from_symbol: unknown element " + s);
}

// ── XYZ reader (simple 2-line header + atom lines) ──────────────────

struct XYZAtom {
    int Z;
    double x, y, z;
};

std::vector<XYZAtom> read_xyz(const fs::path& path) {
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("read_xyz: cannot open " + path.string());
    int natm;
    in >> natm;
    std::string rest;
    std::getline(in, rest);
    std::string comment;
    std::getline(in, comment);
    std::vector<XYZAtom> atoms;
    atoms.reserve(natm);
    for (int i = 0; i < natm; ++i) {
        std::string sym;
        double x, y, z;
        in >> sym >> x >> y >> z;
        atoms.push_back({z_from_symbol(sym), x, y, z});
    }
    return atoms;
}

// ── Reference-JSON slice (only the fields the bench needs) ──────────
//
// We avoid pulling a full JSON dependency; the PySCF generator writes a
// pretty-printed object, so a tiny field extractor is sufficient.

struct Reference {
    std::string method;        // "rhf" | "rks"
    std::string xc;            // "pbe", ... ("null" string if RHF)
    std::string basis;         // "def2-svp" ...
    int spin = 0;              // 2S (PySCF convention)
    double energy_hartree = 0; // optimized energy
};

std::string extract_string(const std::string& text, const std::string& key) {
    const auto pos = text.find("\"" + key + "\"");
    if (pos == std::string::npos)
        return {};
    const auto colon = text.find(':', pos);
    const auto q1 = text.find('"', colon);
    if (q1 == std::string::npos)
        return {};
    const auto q2 = text.find('"', q1 + 1);
    if (q2 == std::string::npos)
        return {};
    return text.substr(q1 + 1, q2 - q1 - 1);
}

double extract_number(const std::string& text, const std::string& key) {
    const auto pos = text.find("\"" + key + "\"");
    if (pos == std::string::npos)
        throw std::runtime_error("reference JSON missing key " + key);
    const auto colon = text.find(':', pos);
    size_t i = colon + 1;
    while (i < text.size() && (text[i] == ' ' || text[i] == '\t'))
        ++i;
    size_t j = i;
    while (j < text.size()
           && (std::isdigit(text[j]) || text[j] == '-' || text[j] == '+'
               || text[j] == '.' || text[j] == 'e' || text[j] == 'E'))
        ++j;
    return std::stod(text.substr(i, j - i));
}

Reference read_reference(const fs::path& path) {
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("read_reference: cannot open "
                                 + path.string());
    std::stringstream ss;
    ss << in.rdbuf();
    const std::string text = ss.str();
    Reference r;
    r.method = extract_string(text, "method");
    r.xc = extract_string(text, "xc");
    r.basis = extract_string(text, "basis");
    r.spin = static_cast<int>(extract_number(text, "spin"));
    r.energy_hartree = extract_number(text, "energy_hartree_optimized");
    return r;
}

cpu::scf::BasisSet basis_from_name(const std::string& b) {
    if (b == "sto-3g" || b == "STO-3G")
        return cpu::scf::BasisSet::STO3G;
    if (b == "6-31g" || b == "6-31G")
        return cpu::scf::BasisSet::Pople631G;
    if (b == "6-31g*" || b == "6-31G*")
        return cpu::scf::BasisSet::Pople631Gs;
    if (b == "def2-svp" || b == "def2-SVP")
        return cpu::scf::BasisSet::Def2SVP;
    if (b == "def2-tzvp" || b == "def2-TZVP")
        return cpu::scf::BasisSet::Def2TZVP;
    throw std::runtime_error("bench_qc_suite: unsupported basis '" + b + "'");
}

cpu::Frame atoms_to_frame(const std::vector<XYZAtom>& atoms) {
    cpu::Frame f(static_cast<int>(atoms.size()));
    for (size_t i = 0; i < atoms.size(); ++i) {
        f.atoms.type[i] = atoms[i].Z;
        f.atoms.x[i] = atoms[i].x;
        f.atoms.y[i] = atoms[i].y;
        f.atoms.z[i] = atoms[i].z;
    }
    f.box = Mat3::cubic(100.0);
    return f;
}

cpu::scf::SCFConfig make_config(const cpu::scf::IntegralPack& ints,
                                const Reference& ref) {
    cpu::scf::SCFConfig c;
    const bool open_shell = ref.spin > 0;
    c.spin = open_shell ? cpu::scf::SpinType::Unrestricted
                        : cpu::scf::SpinType::Restricted;
    if (ref.method == "rks" || ref.method == "uks") {
        c.method = cpu::scf::MeanFieldType::KS;
        c.xc_type =
            (ref.xc == "pbe") ? cpu::scf::XCType::PBE : cpu::scf::XCType::LSDA;
    } else {
        c.method = cpu::scf::MeanFieldType::HF;
    }
    // For open-shell, split electrons using n_occ ± spin/2.
    const int n_total = ints.n_electrons;
    c.n_alpha = (n_total + ref.spin) / 2;
    c.n_beta = n_total - c.n_alpha;
    c.tol = 1e-9;
    c.use_diis = true;
    c.diis_max_vectors = 8;
    c.diis_start_iter = open_shell ? 5 : 1;
    c.damping_factor = open_shell ? 0.5 : 0.0;
    c.init_guess =
        open_shell ? cpu::scf::InitGuess::Hcore : cpu::scf::InitGuess::MINAO;
    // High-accuracy grid: reproduces PySCF's default level=3 grid so the
    // ATV-vs-PySCF energy residue drops into the sub-microHartree regime.
    c.grid_level = cpu::scf::GridLevel::MatchPySCF;
    return c;
}

struct RunResult {
    bool converged;
    double wall_ms;
    double energy_ha;
    double delta_ha;
};

RunResult run_molecule(const std::string& stem, const fs::path& xyz_path,
                       const Reference& ref) {
    auto atoms = read_xyz(xyz_path);
    auto frame = atoms_to_frame(atoms);
    const auto basis_id = basis_from_name(ref.basis);

    auto ints = cpu::scf::build_integrals(frame, basis_id);
    auto cfg = make_config(ints, ref);

    cpu::minimizer::LBFGS step{
        /*memory=*/5,      /*alpha0=*/0.01,
        /*ftol=*/5e-4,     /*max_ls=*/20,
        /*rms_ftol=*/3e-4,
        /*max_dtol=*/5e-3,
        /*rms_dtol=*/3e-3,
        /*etol=*/1e-6,
    };
    step.line_search();

    cpu::minimizer::Minimize opt;
    opt.add(std::move(step));
    opt.on(ATV::minimizer::Phase::force)
        .add(cpu::scf::SCFEvaluator{cfg, basis_id});
    opt.build();

    const auto t0 = std::chrono::steady_clock::now();
    opt.run(frame, 200);
    const double wall_ms = std::chrono::duration<double, std::milli>(
                               std::chrono::steady_clock::now() - t0)
                               .count();

    auto ints_final = cpu::scf::build_integrals(frame, basis_id);
    auto cfg_final = make_config(ints_final, ref);
    cfg_final.tol = 1e-10;
    auto res = cpu::scf::scf_solve(ints_final, cfg_final);

    return {res.converged, wall_ms, res.energy,
            res.energy - ref.energy_hartree};
}

std::vector<std::string> discover_stems(const fs::path& ref_dir) {
    std::vector<std::string> out;
    for (const auto& e : fs::directory_iterator(ref_dir)) {
        if (!e.is_regular_file())
            continue;
        if (e.path().extension() != ".json")
            continue;
        out.push_back(e.path().stem().string());
    }
    std::sort(out.begin(), out.end());
    return out;
}

fs::path find_data_dir() {
#ifdef ATV_QC_DATA_DIR
    const fs::path candidate = ATV_QC_DATA_DIR;
    if (fs::exists(candidate))
        return candidate;
#endif
    // Fall back to cwd-relative.
    return fs::current_path() / "benchmarks" / "qc";
}

} // namespace

int main(int argc, char** argv) {
    const fs::path data_dir = find_data_dir();
    const fs::path mol_dir = data_dir / "molecules";
    const fs::path ref_dir = data_dir / "reference";

    std::vector<std::string> stems;
    if (argc > 1) {
        for (int i = 1; i < argc; ++i)
            stems.emplace_back(argv[i]);
    } else {
        stems = discover_stems(ref_dir);
    }

    std::printf("# Phase 1.2 QC benchmark suite\n"
                "# data_dir=%s\n"
                "# level of theory is per-reference (from the JSON)\n"
                "#\n"
                "# %-14s  %-10s  %-10s  %-16s  %-14s  %-10s\n",
                data_dir.string().c_str(), "mol", "conv", "wall_ms",
                "E_final [Ha]", "dE_vs_ref [Ha]", "status");

    // The GridLevel::MatchPySCF grid (Treutler-Ahlrichs M4 radial +
    // nwchem-style angular pruning bounded by the level-3 maximum Lebedev
    // order) reproduces PySCF's default level=3 grid. The ATV-vs-PySCF
    // def2-SVP energy residue on the H2O/NH3/CH4 triplet is sub-1e-5 Ha;
    // 5e-4 Ha is a conservative whole-suite acceptance bar that still
    // catches genuine integral / SCF regressions across the larger
    // benchmark molecules.
    const double threshold_ha = 5e-4;
    int n_fail = 0;
    int n_run = 0;

    for (const auto& stem : stems) {
        const auto xyz_path = mol_dir / (stem + ".xyz");
        const auto ref_path = ref_dir / (stem + ".json");
        if (!fs::exists(xyz_path)) {
            std::printf("  %-14s  skipped: missing %s\n", stem.c_str(),
                        xyz_path.string().c_str());
            continue;
        }
        if (!fs::exists(ref_path)) {
            std::printf("  %-14s  skipped: missing %s\n", stem.c_str(),
                        ref_path.string().c_str());
            continue;
        }

        Reference ref;
        try {
            ref = read_reference(ref_path);
        } catch (const std::exception& e) {
            std::printf("  %-14s  ERROR: %s\n", stem.c_str(), e.what());
            ++n_fail;
            continue;
        }

        RunResult r;
        try {
            r = run_molecule(stem, xyz_path, ref);
        } catch (const std::exception& e) {
            std::printf("  %-14s  ERROR: %s\n", stem.c_str(), e.what());
            ++n_fail;
            continue;
        }

        const bool within = std::abs(r.delta_ha) <= threshold_ha;
        const char* status =
            !r.converged ? "SCF FAIL" : (within ? "PASS" : "DRIFT");
        std::printf("  %-14s  %-10s  %10.1f  %16.10f  %+14.2e  %-10s\n",
                    stem.c_str(), r.converged ? "yes" : "no", r.wall_ms,
                    r.energy_ha, r.delta_ha, status);
        if (!r.converged || !within)
            ++n_fail;
        ++n_run;
    }

    std::printf("#\n# %d / %d molecules within %.0e Ha of PySCF reference.\n",
                n_run - n_fail, n_run, threshold_ha);
    return (n_fail == 0) ? 0 : 1;
}

#else
#include <cstdio>
int main() {
    std::printf("bench_qc_suite requires -DATV_BUILD_RHF=ON\n");
    return 1;
}
#endif
