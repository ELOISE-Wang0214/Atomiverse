/// @file tests/cpu/vib/test_vib_output_components.cpp
/// @brief Output-component round-trip tests for the vibrational-analysis
///        thermochemistry link (spec phase4-vibrations-03,
///        ac-004 / ac-005 / ac-006).
///
/// After a full `cpu::vib::Hessian` run on H2O, this exercises the two
/// Recorder components:
///   * `FrequencyReport` — its TSV must parse with one data row per
///     frequency and a trailer carrying the count and the ZPE.
///   * `NormalModeWriter` — its per-mode XYZ animations must load with
///     the correct atom count and the input atom symbols.
/// It also verifies the structural-role tags of all three new
/// components (ac-006).

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/vib/frequency_report.h"
#include "atomiverse/cpu/vib/hessian.h"
#include "atomiverse/cpu/vib/hessian_result.h"
#include "atomiverse/cpu/vib/normal_mode_writer.h"
#include "atomiverse/cpu/vib/thermochem.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

using namespace ATV;

namespace {

/// @brief H2O at its RHF/6-31G optimized geometry (Angstrom).
cpu::Frame make_water() {
    cpu::Frame f(3);
    f.atoms.x = {-0.000000000000, -0.000000000000, 0.000000000000};
    f.atoms.y = {-0.000000000000, 0.785166496369, -0.785166496369};
    f.atoms.z = {0.034902683452, 0.569048658274, 0.569048658274};
    f.atoms.type = {8, 1, 1}; // atomic numbers Z.
    f.box = Mat3::cubic(100.0);
    return f;
}

/// @brief Run cpu::vib::Hessian on H2O with custom output paths.
/// @param freq_path  TSV output path for the FrequencyReport.
/// @param mode_prefix XYZ prefix for the NormalModeWriter.
/// @return The HessianResult (frequencies + thermochemistry populated).
cpu::vib::HessianResult run_water(const std::string& freq_path,
                                  const std::string& mode_prefix) {
    cpu::Frame frame = make_water();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::Pople631G);
    auto cfg = cpu::scf::rhf_config(ints, /*tol=*/1e-12);
    cfg.conv_tol_grad = 1e-9;

    // Build the driver with H2O's symmetry number sigma = 2 and the
    // test-controlled output paths for the two Recorders.
    cpu::vib::Thermochem::Config tc;
    tc.symmetry_number = 2;
    cpu::vib::Hessian::OutputConfig out;
    out.frequency_tsv = freq_path;
    out.mode_xyz_prefix = mode_prefix;
    cpu::vib::Hessian driver{tc, out};
    driver.on(ATV::vib::Phase::force)
        .add(cpu::scf::SCFEvaluator{cfg, cpu::scf::BasisSet::Pople631G});
    driver.build();
    return driver.run(frame);
}

} // namespace

// ── ac-006: structural-role tags ────────────────────────────────────────

TEST(VibOutputComponents, ComponentRoles) {
    cpu::vib::Thermochem tc{cpu::vib::Thermochem::Config{}};
    cpu::vib::FrequencyReport fr{std::string("frequencies.tsv")};
    cpu::vib::NormalModeWriter nm{std::string("mode")};

    EXPECT_EQ(tc.role(), ComponentRole::Observer);
    EXPECT_EQ(fr.role(), ComponentRole::Recorder);
    EXPECT_EQ(nm.role(), ComponentRole::Recorder);
}

// ── ac-004: FrequencyReport writes a parseable TSV ──────────────────────

TEST(VibOutputComponents, FrequencyReportTsv) {
    namespace fs = std::filesystem;
    const std::string tsv =
        (fs::temp_directory_path() / "atv_vib_h2o_freqs.tsv").string();
    const std::string prefix =
        (fs::temp_directory_path() / "atv_vib_h2o_modeA").string();
    fs::remove(tsv);

    const auto result = run_water(tsv, prefix);

    ASSERT_TRUE(fs::exists(tsv)) << "FrequencyReport did not write " << tsv;

    std::ifstream is(tsv);
    ASSERT_TRUE(is.good());

    std::string line;
    // First line is the header.
    ASSERT_TRUE(std::getline(is, line));
    EXPECT_NE(line.find("frequency"), std::string::npos);

    int n_data_rows = 0;
    double trailer_zpe = 0.0;
    int trailer_n_modes = -1;
    bool saw_zpe = false, saw_count = false;
    while (std::getline(is, line)) {
        if (line.empty())
            continue;
        std::istringstream ss(line);
        if (line[0] == '#') {
            // Trailer comment: "# key<TAB>value".
            std::string hash, key;
            ss >> hash >> key;
            if (key == "zpe_hartree") {
                ss >> trailer_zpe;
                saw_zpe = true;
            }
            if (key == "n_modes") {
                ss >> trailer_n_modes;
                saw_count = true;
            }
        } else {
            // Data row: index<TAB>frequency<TAB>type.
            int idx;
            double freq;
            std::string type;
            ss >> idx >> freq >> type;
            EXPECT_FALSE(ss.fail()) << "unparseable TSV data row: " << line;
            ++n_data_rows;
        }
    }

    // One data row per frequency (3N-6 = 3 modes for H2O).
    EXPECT_EQ(n_data_rows, static_cast<int>(result.frequencies.size()));
    EXPECT_EQ(n_data_rows, 3);

    // Trailer carries the correct count and the ZPE.
    ASSERT_TRUE(saw_count) << "TSV trailer missing n_modes";
    ASSERT_TRUE(saw_zpe) << "TSV trailer missing zpe_hartree";
    EXPECT_EQ(trailer_n_modes, static_cast<int>(result.frequencies.size()));
    EXPECT_NEAR(trailer_zpe, result.zpe, 1e-12);

    fs::remove(tsv);
}

// ── ac-005: NormalModeWriter writes a valid loadable XYZ animation ──────

TEST(VibOutputComponents, NormalModeWriterXyz) {
    namespace fs = std::filesystem;
    const std::string tsv =
        (fs::temp_directory_path() / "atv_vib_h2o_freqs_b.tsv").string();
    const std::string prefix =
        (fs::temp_directory_path() / "atv_vib_h2o_modeB").string();

    // Clean any stale per-mode files.
    for (int k = 0; k < 12; ++k)
        fs::remove(prefix + "_mode" + std::to_string(k) + ".xyz");

    const auto result = run_water(tsv, prefix);
    const int n_modes = static_cast<int>(result.frequencies.size());
    ASSERT_EQ(n_modes, 3) << "H2O must yield 3N-6 = 3 real modes";

    // Each real mode produced a multi-frame XYZ file.
    for (int k = 0; k < n_modes; ++k) {
        const std::string xyz = prefix + "_mode" + std::to_string(k) + ".xyz";
        ASSERT_TRUE(fs::exists(xyz))
            << "NormalModeWriter did not write " << xyz;

        std::ifstream is(xyz);
        ASSERT_TRUE(is.good());

        int n_frames = 0;
        std::string line;
        while (std::getline(is, line)) {
            if (line.empty())
                continue;
            // Each XYZ frame: count line, comment line, N atom lines.
            const int natoms = std::stoi(line);
            EXPECT_EQ(natoms, 3) << "each frame must have the H2O atom count";

            ASSERT_TRUE(std::getline(is, line)); // comment line.

            for (int i = 0; i < natoms; ++i) {
                ASSERT_TRUE(std::getline(is, line));
                std::istringstream ss(line);
                std::string sym;
                double x, y, z;
                ss >> sym >> x >> y >> z;
                EXPECT_FALSE(ss.fail())
                    << "unparseable XYZ atom line: " << line;
                // Atom symbols must match the H2O input (O, H, H).
                if (i == 0)
                    EXPECT_EQ(sym, "O");
                else
                    EXPECT_EQ(sym, "H");
            }
            ++n_frames;
        }
        EXPECT_GT(n_frames, 1)
            << "mode " << k << " must be a multi-frame animation";
    }

    for (int k = 0; k < n_modes; ++k)
        fs::remove(prefix + "_mode" + std::to_string(k) + ".xyz");
    fs::remove(tsv);
}
