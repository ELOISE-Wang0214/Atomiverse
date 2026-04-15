/// @file examples/tutorial/03_uhf.cpp
/// @brief Tutorial example for unrestricted Hartree-Fock.
///
/// Demonstrates the unrestricted Hartree-Fock capability on three cases:
///
///   1. H atom (n_alpha=1, n_beta=0): only solvable with UHF.
///   2. H2 at equilibrium (0.74 A): RHF ~ UHF — no symmetry breaking.
///   3. H2 bond-length scan: UHF breaks spin symmetry at stretched geometries,
///      giving a lower energy than RHF and nonzero spin contamination.
///
/// Molecular visualization output (always produced):
///   h_atom.xyz          : H atom extXYZ with Mulliken population
///   h2_eq_rhf.xyz       : H2 @ 0.74 A, RHF
///   h2_eq_uhf.xyz       : H2 @ 0.74 A, UHF
///   h2_scan_rhf.xyz     : multi-frame extXYZ: RHF PES scan (open in OVITO/Avogadro)
///   h2_scan_uhf.xyz     : multi-frame extXYZ: UHF PES scan
///
/// With ATV_USE_MATPLOT=ON, saves two figures:
///   h2_pes.png          : RHF vs UHF energy curves with exact dissociation limit
///   h2_spin_contam.png  : spin contamination <S^2> - S*(S+1) vs bond length
///
/// Build without plots:
///   cmake -B build -DATV_BUILD_CUDA=OFF -DATV_BUILD_RHF=ON -DATV_BUILD_EXAMPLES=ON
///   cmake --build build --target 03_uhf
///
/// Build with plots (requires gnuplot):
///   cmake -B build -DATV_BUILD_CUDA=OFF -DATV_BUILD_RHF=ON \
///         -DATV_BUILD_EXAMPLES=ON -DATV_USE_MATPLOT=ON
///   cmake --build build --target 03_uhf

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/molrec.h"
#include "atomiverse/cpu/observer/mulliken.h"
#include "atomiverse/cpu/recorder/xyz.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cpu/scf/uhf_solve.h"

#ifdef ATV_HAS_MATPLOT
#include <matplot/matplot.h>
#endif

using namespace ATV;

// ── Geometry helpers ─────────────────────────────────────────────────────────

static cpu::Frame single_atom(int z) {
    cpu::Frame f(1);
    f.atoms.x[0] = double(0); f.atoms.y[0] = double(0); f.atoms.z[0] = double(0);
    f.atoms.type[0] = z;
    f.box = Mat3::cubic(100.0);
    return f;
}

static cpu::Frame diatomic_h2(double r_angstrom) {
    cpu::Frame f(2);
    f.atoms.x[0] = double(0);           f.atoms.y[0] = double(0); f.atoms.z[0] = double(0);
    f.atoms.x[1] = double(r_angstrom);  f.atoms.y[1] = double(0); f.atoms.z[1] = double(0);
    f.atoms.type[0] = 1; f.atoms.type[1] = 1;
    f.box = Mat3::cubic(100.0);
    return f;
}

// ── Plotting ─────────────────────────────────────────────────────────────────

#ifdef ATV_HAS_MATPLOT

static void plot_pes(const std::vector<double>& r_vals,
                     const std::vector<double>& e_rhf,
                     const std::vector<double>& e_uhf)
{
    auto fig = matplot::figure(true);
    fig->size(780, 600);

    auto p_rhf = matplot::plot(r_vals, e_rhf);
    p_rhf->line_width(2.0);
    p_rhf->color({0.5f, 0.5f, 0.5f, 1.0f});
    p_rhf->display_name("RHF");

    matplot::hold(matplot::on);

    auto p_uhf = matplot::plot(r_vals, e_uhf);
    p_uhf->line_width(2.0);
    p_uhf->color({0.2f, 0.4f, 0.8f, 1.0f});
    p_uhf->display_name("UHF");

    constexpr double kExactLimit = 2.0 * (-0.4666);
    std::vector<double> ref_y(r_vals.size(), kExactLimit);
    auto p_ref = matplot::plot(r_vals, ref_y);
    p_ref->line_width(1.5);
    p_ref->line_style("--");
    p_ref->color({0.8f, 0.2f, 0.2f, 1.0f});
    p_ref->display_name("2 \\times E_H (exact limit)");

    auto ax = matplot::gca();
    ax->font_size(11);
    matplot::xlabel("R(H-H) / \\AA");
    matplot::ylabel("E / Hartree");
    matplot::title("H_{2} / STO-3G: RHF vs UHF Potential Energy Surface");
    matplot::legend()->location(matplot::legend::general_alignment::topright);

    matplot::hold(matplot::off);
    matplot::save("h2_pes.png");
    std::cout << "  Saved h2_pes.png\n";
}

static void plot_spin_contamination(const std::vector<double>& r_vals,
                                    const std::vector<double>& spin_contam)
{
    auto fig = matplot::figure(true);
    fig->size(780, 600);

    auto p = matplot::plot(r_vals, spin_contam);
    p->line_width(2.0);
    p->color({0.2f, 0.4f, 0.8f, 1.0f});
    p->marker(matplot::line_spec::marker_style::circle);
    p->marker_size(4);
    p->marker_face_color({0.2f, 0.4f, 0.8f, 1.0f});

    matplot::hold(matplot::on);
    std::vector<double> zero_line(r_vals.size(), 0.0);
    std::vector<double> one_line(r_vals.size(), 1.0);
    auto p0 = matplot::plot(r_vals, zero_line);
    p0->line_width(1.0); p0->line_style("--"); p0->color({0.5f, 0.5f, 0.5f, 0.6f});
    auto p1 = matplot::plot(r_vals, one_line);
    p1->line_width(1.0); p1->line_style("--"); p1->color({0.8f, 0.2f, 0.2f, 0.6f});

    auto ax = matplot::gca();
    ax->font_size(11);
    matplot::xlabel("R(H-H) / \\AA");
    matplot::ylabel("\\langle S^2 \\rangle - S(S+1)");
    matplot::title("H_{2} / STO-3G: UHF Spin Contamination");

    matplot::hold(matplot::off);
    matplot::save("h2_spin_contam.png");
    std::cout << "  Saved h2_spin_contam.png\n";
}

static void plot_mulliken_scan(const std::vector<double>& r_vals,
                               const std::vector<double>& q_rhf,
                               const std::vector<double>& q_uhf_a,
                               const std::vector<double>& spin_a)
{
    // Fig 3: Mulliken charge on atom 0 vs R — shows charge delocalization
    auto fig = matplot::figure(true);
    fig->size(780, 600);

    auto p1 = matplot::plot(r_vals, q_rhf);
    p1->line_width(2.0); p1->color({0.5f, 0.5f, 0.5f, 1.0f});
    p1->display_name("RHF q(H_A)");

    matplot::hold(matplot::on);

    auto p2 = matplot::plot(r_vals, q_uhf_a);
    p2->line_width(2.0); p2->color({0.2f, 0.4f, 0.8f, 1.0f});
    p2->display_name("UHF q(H_A)");

    auto p3 = matplot::plot(r_vals, spin_a);
    p3->line_width(2.0); p3->line_style("-."); p3->color({0.8f, 0.4f, 0.0f, 1.0f});
    p3->display_name("UHF spin(H_A)");

    // Neutral reference
    std::vector<double> zero(r_vals.size(), 0.0);
    auto p0 = matplot::plot(r_vals, zero);
    p0->line_width(1.0); p0->line_style("--"); p0->color({0.7f, 0.7f, 0.7f, 0.6f});

    auto ax = matplot::gca();
    ax->font_size(11);
    matplot::xlabel("R(H-H) / \\AA");
    matplot::ylabel("Mulliken charge / spin [e]");
    matplot::title("H_{2} / STO-3G: Mulliken Population vs Bond Length");
    matplot::legend()->location(matplot::legend::general_alignment::topright);

    matplot::hold(matplot::off);
    matplot::save("h2_mulliken.png");
    std::cout << "  Saved h2_mulliken.png\n";
}

#endif // ATV_HAS_MATPLOT

// ── Main ─────────────────────────────────────────────────────────────────────

int main() {
    const double tol = (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-8);
    const int max_iter_atoms = 50;
    const int max_iter_h2    = 200;
    const int max_iter_scan  = 400;

    // ── 1. H atom ────────────────────────────────────────────────────────────
    std::cout << "=== Case 1: H atom / STO-3G / UHF ===\n";
    {
        auto frame  = single_atom(1);
        auto ints   = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
        auto result = cpu::scf::uhf_solve(ints, 1, 0, tol, max_iter_atoms);

        std::cout << "  Converged:       " << (result.converged ? "yes" : "no")
                  << " (" << result.iterations << " iterations)\n"
                  << std::fixed << std::setprecision(8)
                  << "  E_total:         " << static_cast<double>(result.energy)
                  << " Hartree  (reference: -0.4666)\n"
                  << std::setprecision(6)
                  << "  Spin contam.:    " << static_cast<double>(result.spin_contamination) << "\n";

        cpu::MolRec rec;
        rec.set_geometry(frame);
        rec.add_string("method", "UHF");
        rec.add_scalar("energy",            result.energy);
        rec.add_scalar("iterations",        double(result.iterations));
        rec.add_scalar("n_alpha",           double(result.n_alpha));
        rec.add_scalar("n_beta",            double(result.n_beta));
        rec.add_scalar("spin_contamination", result.spin_contamination);
        xt::xtensor<double,2> P_total = result.P_alpha + result.P_beta;
        rec.add_field("charge",       cpu::mulliken_charge(ints, P_total));
        rec.add_field("spin_density", cpu::mulliken_spin(ints, result.P_alpha, result.P_beta));
        rec.commit_frame();
        cpu::print_summary(rec);
        cpu::XYZ("h_atom.xyz").write_ext(rec);
        std::cout << "  Wrote h_atom.xyz\n\n";
    }

    // ── 2. H2 at equilibrium: RHF vs UHF ─────────────────────────────────────
    std::cout << "=== Case 2: H2 (0.74 Ang) / STO-3G / RHF vs UHF ===\n";
    {
        auto frame = diatomic_h2(0.74);
        auto ints  = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
        auto rhf   = cpu::scf::rhf_solve(ints, tol, 50);
        auto uhf   = cpu::scf::uhf_solve(ints, 1, 1, tol, max_iter_h2);

        std::cout << std::fixed << std::setprecision(8)
                  << "  RHF:  E = " << static_cast<double>(rhf.energy)
                  << " Ha  iter=" << rhf.iterations << "\n"
                  << "  UHF:  E = " << static_cast<double>(uhf.energy)
                  << " Ha  iter=" << uhf.iterations
                  << std::setprecision(4)
                  << "  spin_contam=" << static_cast<double>(uhf.spin_contamination) << "\n"
                  << std::scientific << std::setprecision(2)
                  << "  |dE|: " << std::fabs(static_cast<double>(uhf.energy - rhf.energy))
                  << " Ha  (should be ~0 at equilibrium)\n\n";

        {
            cpu::MolRec rhf_rec;
            rhf_rec.set_geometry(frame);
            rhf_rec.add_string("method", "RHF");
            rhf_rec.add_scalar("energy",     rhf.energy);
            rhf_rec.add_scalar("iterations", double(rhf.iterations));
            rhf_rec.add_scalar("n_alpha",    double(ints.n_occ()));
            rhf_rec.add_scalar("n_beta",     double(ints.n_occ()));
            rhf_rec.add_scalar("spin_contamination", double(0));
            rhf_rec.add_field("charge",       cpu::mulliken_charge(ints, rhf.P));
            rhf_rec.add_field("spin_density", std::vector<double>(frame.atoms.size(), double(0)));
            rhf_rec.commit_frame();
            cpu::print_summary(rhf_rec);
            cpu::XYZ("h2_eq_rhf.xyz").write_ext(rhf_rec);
        }
        {
            cpu::MolRec uhf_rec;
            uhf_rec.set_geometry(frame);
            uhf_rec.add_string("method", "UHF");
            uhf_rec.add_scalar("energy",     uhf.energy);
            uhf_rec.add_scalar("iterations", double(uhf.iterations));
            uhf_rec.add_scalar("n_alpha",    double(uhf.n_alpha));
            uhf_rec.add_scalar("n_beta",     double(uhf.n_beta));
            uhf_rec.add_scalar("spin_contamination", uhf.spin_contamination);
            xt::xtensor<double,2> P_total = uhf.P_alpha + uhf.P_beta;
            uhf_rec.add_field("charge",       cpu::mulliken_charge(ints, P_total));
            uhf_rec.add_field("spin_density", cpu::mulliken_spin(ints, uhf.P_alpha, uhf.P_beta));
            uhf_rec.commit_frame();
            cpu::print_summary(uhf_rec);
            cpu::XYZ("h2_eq_uhf.xyz").write_ext(uhf_rec);
        }
        std::cout << "  Wrote h2_eq_rhf.xyz, h2_eq_uhf.xyz\n\n";
    }

    // ── 3. H2 PES scan: RHF vs UHF ───────────────────────────────────────────
    std::cout << "=== Case 3: H2 PES scan / STO-3G ===\n"
              << std::left
              << std::setw(8) << "R (Ang)"
              << std::setw(15) << "E_RHF"
              << std::setw(15) << "E_UHF"
              << std::setw(13) << "dE"
              << std::setw(13) << "<S^2>_c"
              << std::setw(11) << "q_H(RHF)"
              << std::setw(11) << "q_H(UHF)" << "\n"
              << std::setw(8) << "--------"
              << std::setw(15) << "--------"
              << std::setw(15) << "--------"
              << std::setw(13) << "--------"
              << std::setw(13) << "--------"
              << std::setw(11) << "----------"
              << std::setw(11) << "----------" << "\n"
              << std::right;

    std::vector<double> r_vals, e_rhf_vals, e_uhf_vals, spin_contam_vals;
    std::vector<double> q_rhf_a, q_uhf_a, spin_uhf_a;  // Mulliken on atom 0

    cpu::MolRec rhf_scan, uhf_scan;
    cpu::XYZ    rhf_xyz("h2_scan_rhf.xyz");
    cpu::XYZ    uhf_xyz("h2_scan_uhf.xyz");

    for (double r = 0.5; r <= 4.01; r += 0.25) {
        auto frame = diatomic_h2(r);
        auto ints  = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
        auto rhf   = cpu::scf::rhf_solve(ints, tol, 100);
        auto uhf   = cpu::scf::uhf_solve(ints, 1, 1, tol, max_iter_scan);

        // Build RHF frame and stream to XYZ
        rhf_scan.set_geometry(frame);
        rhf_scan.add_string("method", "RHF");
        rhf_scan.add_scalar("energy",     rhf.energy);
        rhf_scan.add_scalar("iterations", double(rhf.iterations));
        rhf_scan.add_scalar("n_alpha",    double(ints.n_occ()));
        rhf_scan.add_scalar("n_beta",     double(ints.n_occ()));
        rhf_scan.add_scalar("spin_contamination", double(0));
        auto rhf_charge = cpu::mulliken_charge(ints, rhf.P);
        rhf_scan.add_field("charge",       rhf_charge);
        rhf_scan.add_field("spin_density", std::vector<double>(frame.atoms.size(), double(0)));
        rhf_scan.commit_frame();
        rhf_xyz.append_ext(rhf_scan);

        // Build UHF frame and stream to XYZ
        xt::xtensor<double,2> P_total = uhf.P_alpha + uhf.P_beta;
        auto uhf_charge = cpu::mulliken_charge(ints, P_total);
        auto uhf_spin   = cpu::mulliken_spin(ints, uhf.P_alpha, uhf.P_beta);
        uhf_scan.set_geometry(frame);
        uhf_scan.add_string("method", "UHF");
        uhf_scan.add_scalar("energy",     uhf.energy);
        uhf_scan.add_scalar("iterations", double(uhf.iterations));
        uhf_scan.add_scalar("n_alpha",    double(uhf.n_alpha));
        uhf_scan.add_scalar("n_beta",     double(uhf.n_beta));
        uhf_scan.add_scalar("spin_contamination", uhf.spin_contamination);
        uhf_scan.add_field("charge",       uhf_charge);
        uhf_scan.add_field("spin_density", uhf_spin);
        uhf_scan.commit_frame();
        uhf_xyz.append_ext(uhf_scan);

        double dE = static_cast<double>(rhf.energy - uhf.energy);
        double qr = static_cast<double>(rhf_charge[0]);
        double qu = static_cast<double>(uhf_charge[0]);
        std::cout << std::fixed
                  << std::setw(8) << std::setprecision(3) << r
                  << std::setw(15) << std::setprecision(8) << static_cast<double>(rhf.energy)
                  << std::setw(15) << static_cast<double>(uhf.energy)
                  << std::setw(13) << std::setprecision(6) << dE
                  << std::setw(13) << static_cast<double>(uhf.spin_contamination)
                  << std::setw(11) << qr
                  << std::setw(11) << qu << "\n";

        r_vals.push_back(r);
        e_rhf_vals.push_back(static_cast<double>(rhf.energy));
        e_uhf_vals.push_back(static_cast<double>(uhf.energy));
        spin_contam_vals.push_back(static_cast<double>(uhf.spin_contamination));
        q_rhf_a.push_back(qr);
        q_uhf_a.push_back(qu);
        spin_uhf_a.push_back(static_cast<double>(uhf_spin[0]));
    }
    std::cout << "\n  Wrote h2_scan_rhf.xyz (" << r_vals.size() << " frames)\n"
              << "  Wrote h2_scan_uhf.xyz (" << r_vals.size() << " frames)\n"
              << "  Open in OVITO or Avogadro to animate the PES scan.\n"
              << std::fixed << std::setprecision(8)
              << "  Exact dissociation limit (2 x E_H): " << (2.0 * (-0.4666)) << " Hartree\n\n";

    // ── Plots ─────────────────────────────────────────────────────────────────
#ifdef ATV_HAS_MATPLOT
    std::cout << "Generating figures...\n";
    plot_pes(r_vals, e_rhf_vals, e_uhf_vals);
    plot_spin_contamination(r_vals, spin_contam_vals);
    plot_mulliken_scan(r_vals, q_rhf_a, q_uhf_a, spin_uhf_a);
    std::cout << "\nFigure guide:\n"
              << "  h2_pes.png       : RHF/UHF energy curves + exact dissociation limit\n"
              << "  h2_spin_contam.png: UHF spin contamination (0 at eq, ->1 at dissociation)\n"
              << "  h2_mulliken.png  : Mulliken charges + spin density on H_A vs R\n"
              << "                     Charges are ~0 at eq (symmetric); UHF spin density\n"
              << "                     rises to ~+0.5 as alpha localizes on one atom.\n";
#else
    std::cout << "(Build with -DATV_USE_MATPLOT=ON for h2_pes.png, h2_spin_contam.png, h2_mulliken.png)\n";
#endif

    std::cout << "\n=== Done ===\n";
    return 0;
}
