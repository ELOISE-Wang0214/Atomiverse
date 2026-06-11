/// @file cpu/io/lammps_data.cpp
/// @brief Self-contained C++ reader for LAMMPS atom_style=molecular data
///        files, sufficient for the MSU PE benchmark.
///
/// We parse the file directly here rather than going through molrs because
/// molrs's atom-style auto-detection misclassifies the MSU file: it sees a
/// 6-column atoms line whose third token parses as a float, and routes it
/// through the `charge` branch, which then mislabels mol_id as `type` in
/// the resulting frame. The MSU file format is simple enough (single atom
/// type, orthogonal box, standard sections) that an in-house parser is
/// cleaner than fighting that heuristic.

#include "atomiverse/cpu/io/lammps_data.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame/frame.h"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace ATV {
namespace cpu {

namespace {

bool starts_with(const std::string& s, const char* prefix) {
    const std::size_t n = std::char_traits<char>::length(prefix);
    return s.size() >= n && std::strncmp(s.data(), prefix, n) == 0;
}

std::string ltrim(const std::string& s) {
    std::size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i])))
        ++i;
    return s.substr(i);
}

// Throw with file context.
[[noreturn]] void die(const std::string& path, const std::string& msg) {
    throw std::runtime_error("load_lammps_data(" + path + "): " + msg);
}

// Read the next non-blank, non-comment line. Returns false at EOF.
bool next_data_line(std::ifstream& in, std::string& out) {
    while (std::getline(in, out)) {
        const std::string t = ltrim(out);
        if (t.empty() || t[0] == '#')
            continue;
        out = t;
        return true;
    }
    return false;
}

} // namespace

Frame load_lammps_data(const std::string& path) {
    std::ifstream in(path);
    if (!in)
        die(path, "cannot open file");

    // ── Pass 1: header ─────────────────────────────────────────────────
    int n_atoms = 0, n_bonds = 0, n_angles = 0, n_dihedrals = 0;
    double xlo = 0.0, xhi = 0.0, ylo = 0.0, yhi = 0.0, zlo = 0.0, zhi = 0.0;

    std::string line;
    // Skip title line.
    std::getline(in, line);

    while (next_data_line(in, line)) {
        // Section headers end the header block.
        if (starts_with(line, "Masses") || starts_with(line, "Atoms")
            || starts_with(line, "Bonds") || starts_with(line, "Angles")
            || starts_with(line, "Dihedrals") || starts_with(line, "Impropers")
            || starts_with(line, "Velocities")) {
            // Push the section line back by seeking; simpler: process it as
            // the section dispatcher below.
            break;
        }
        std::istringstream iss(line);
        if (line.find("atoms") != std::string::npos
            && line.find("atom types") == std::string::npos) {
            iss >> n_atoms;
        } else if (line.find("bonds") != std::string::npos
                   && line.find("bond types") == std::string::npos) {
            iss >> n_bonds;
        } else if (line.find("angles") != std::string::npos
                   && line.find("angle types") == std::string::npos) {
            iss >> n_angles;
        } else if (line.find("dihedrals") != std::string::npos
                   && line.find("dihedral types") == std::string::npos) {
            iss >> n_dihedrals;
        } else if (line.find("xlo xhi") != std::string::npos) {
            iss >> xlo >> xhi;
        } else if (line.find("ylo yhi") != std::string::npos) {
            iss >> ylo >> yhi;
        } else if (line.find("zlo zhi") != std::string::npos) {
            iss >> zlo >> zhi;
        }
        // Ignore "N atom types", "N bond types", etc. — single type assumed.
    }

    if (n_atoms <= 0)
        die(path, "header reports zero atoms");

    Frame frame(n_atoms);

    std::unordered_map<int, double> type_to_mass;

    // ── Pass 2: sections — `line` currently holds the first section header
    auto consume_n = [&](int n, auto handler) {
        int read = 0;
        while (read < n) {
            if (!next_data_line(in, line))
                die(path, "unexpected EOF in section");
            handler(line);
            ++read;
        }
    };

    auto is_section = [&](const std::string& s) {
        return starts_with(line, s.c_str());
    };

    // LAMMPS atom-id (1-indexed) → 0-indexed position in `frame.atoms.*`.
    std::unordered_map<int, int> id_to_idx;
    id_to_idx.reserve(n_atoms);
    int next_idx = 0;

    do {
        if (is_section("Masses")) {
            // Masses: <type> <mass>
            // Type count is unbounded in the header; consume until next header.
            while (next_data_line(in, line)) {
                if (starts_with(line, "Atoms") || starts_with(line, "Bonds")
                    || starts_with(line, "Angles")
                    || starts_with(line, "Dihedrals")
                    || starts_with(line, "Velocities")
                    || starts_with(line, "Impropers")) {
                    break;
                }
                std::istringstream iss(line);
                int t;
                double m;
                if (iss >> t >> m)
                    type_to_mass[t] = m;
            }
            continue; // `line` now holds the next section header
        }
        if (is_section("Atoms")) {
            // atom_style=molecular: id mol-id type x y z
            consume_n(n_atoms, [&](const std::string& l) {
                std::istringstream iss(l);
                int id, mol_id, atype;
                double x, y, z;
                if (!(iss >> id >> mol_id >> atype >> x >> y >> z)) {
                    die(path, "malformed Atoms row: " + l);
                }
                const int idx = next_idx++;
                id_to_idx[id] = idx;
                frame.atoms.x[idx] = x;
                frame.atoms.y[idx] = y;
                frame.atoms.z[idx] = z;
                frame.atoms.type[idx] = atype - 1; // 0-indexed
                auto it = type_to_mass.find(atype);
                frame.atoms.mass[idx] =
                    (it != type_to_mass.end()) ? it->second : 0.0;
            });
            // Fetch next header.
            if (!next_data_line(in, line))
                line.clear();
            continue;
        }
        // Topology sections — scatter into per-atom layout on `frame.atoms`.
        // `id_to_idx` maps the LAMMPS atom-id (1-indexed) to the row in
        // frame.atoms (0-indexed); the row IS the global tag (Frame ctor
        // sets tag[i] == i). The `add_*` helpers route the entry to the
        // LAMMPS owner atom's slot via `atoms.map`.
        if (is_section("Bonds")) {
            frame.atoms.map_clear();
            frame.atoms.map_set(0, frame.atoms.nlocal);
            consume_n(n_bonds, [&](const std::string& l) {
                std::istringstream iss(l);
                int bid, btype, ai, aj;
                if (!(iss >> bid >> btype >> ai >> aj)) {
                    die(path, "malformed Bonds row: " + l);
                }
                frame.atoms.add_bond(id_to_idx[ai], id_to_idx[aj], btype - 1);
            });
            if (!next_data_line(in, line))
                line.clear();
            continue;
        }
        if (is_section("Angles")) {
            consume_n(n_angles, [&](const std::string& l) {
                std::istringstream iss(l);
                int aid, atype, ai, aj, ak;
                if (!(iss >> aid >> atype >> ai >> aj >> ak)) {
                    die(path, "malformed Angles row: " + l);
                }
                frame.atoms.add_angle(id_to_idx[ai], id_to_idx[aj],
                                      id_to_idx[ak], atype - 1);
            });
            if (!next_data_line(in, line))
                line.clear();
            continue;
        }
        if (is_section("Dihedrals")) {
            consume_n(n_dihedrals, [&](const std::string& l) {
                std::istringstream iss(l);
                int did, dtype, ai, aj, ak, al;
                if (!(iss >> did >> dtype >> ai >> aj >> ak >> al)) {
                    die(path, "malformed Dihedrals row: " + l);
                }
                frame.atoms.add_dihedral(id_to_idx[ai], id_to_idx[aj],
                                         id_to_idx[ak], id_to_idx[al],
                                         dtype - 1);
            });
            if (!next_data_line(in, line))
                line.clear();
            continue;
        }
        // Skip Velocities / Impropers — not used in this benchmark.
        if (is_section("Velocities") || is_section("Impropers")) {
            // Best-effort skip: advance until the next section header or EOF.
            while (next_data_line(in, line)) {
                if (starts_with(line, "Atoms") || starts_with(line, "Bonds")
                    || starts_with(line, "Angles")
                    || starts_with(line, "Dihedrals")
                    || starts_with(line, "Masses")
                    || starts_with(line, "Velocities")
                    || starts_with(line, "Impropers"))
                    break;
            }
            continue;
        }
        // Unknown section header — bail out.
        break;
    } while (!line.empty());

    // ── Box: orthogonal ────────────────────────────────────────────────
    Mat3 cell;
    cell(0, 0) = xhi - xlo;
    cell(1, 1) = yhi - ylo;
    cell(2, 2) = zhi - zlo;
    frame.box = cell;

    // Build 1-2 / 1-3 / 1-4 exclusion list once topology is scattered.
    frame.atoms.build_special_list();

    return frame;
}

} // namespace cpu
} // namespace ATV
