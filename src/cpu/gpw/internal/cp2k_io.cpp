/// @file src/cpu/gpw/internal/cp2k_io.cpp
/// @brief Implementation of CP2K GTH potential file parsing (see header).

#include "atomiverse/cpu/gpw/internal/cp2k_io.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ATV {
namespace cpu {
namespace gpw {

namespace {

/// @brief Read the next content line, skipping blank and '#'-comment lines.
/// @return false at end of file.
bool next_content_line(std::istream& in, std::string& line) {
    while (std::getline(in, line)) {
        std::size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue; // blank
        if (line[first] == '#') continue;         // comment
        return true;
    }
    return false;
}

std::vector<std::string> tokens(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream iss(line);
    std::string t;
    while (iss >> t) out.push_back(t);
    return out;
}

[[noreturn]] void fail(const std::string& msg) {
    throw std::runtime_error("parse_gth_potential: " + msg);
}

double to_d(const std::string& s) { return std::stod(s); }
int to_i(const std::string& s) { return std::stoi(s); }

} // namespace

GthPotential parse_gth_potential(const std::string& path, const std::string& element,
                                 const std::string& name) {
    std::ifstream in(path);
    if (!in) fail("cannot open file '" + path + "'");

    std::string line;
    // Locate the header line: first token == element AND `name` is one of the
    // following (alias) tokens. Data lines start with a number, so they never
    // match an element symbol.
    bool found = false;
    while (next_content_line(in, line)) {
        auto tk = tokens(line);
        if (tk.empty() || tk[0] != element) continue;
        for (std::size_t i = 1; i < tk.size(); ++i) {
            if (tk[i] == name) {
                found = true;
                break;
            }
        }
        if (found) break;
    }
    if (!found) fail("no entry for element '" + element + "' name '" + name + "' in " + path);

    GthPotential pot;
    pot.element = element;
    pot.name = name;

    // Occupation line: valence electrons per l-shell; sum is Zion.
    if (!next_content_line(in, line)) fail("unexpected EOF reading occupation");
    for (const auto& t : tokens(line)) pot.zion += to_i(t);

    // Local part: rloc, nexp, then nexp coefficients (may wrap, but CP2K keeps
    // them on one line for the parameterizations we target).
    if (!next_content_line(in, line)) fail("unexpected EOF reading local part");
    {
        auto tk = tokens(line);
        if (tk.size() < 2) fail("malformed local line: '" + line + "'");
        pot.rloc = to_d(tk[0]);
        int nexp = to_i(tk[1]);
        if (static_cast<int>(tk.size()) < 2 + nexp)
            fail("local line lists " + std::to_string(nexp) +
                 " coeffs but has fewer fields: '" + line + "'");
        for (int i = 0; i < nexp; ++i) pot.cloc.push_back(to_d(tk[2 + i]));
    }

    // Number of nonlocal projector channels.
    if (!next_content_line(in, line)) fail("unexpected EOF reading channel count");
    int nchannels = to_i(tokens(line).at(0));

    for (int l = 0; l < nchannels; ++l) {
        if (!next_content_line(in, line)) fail("unexpected EOF reading channel header");
        auto tk = tokens(line);
        if (tk.size() < 2) fail("malformed channel line: '" + line + "'");
        GthNonlocalChannel ch;
        ch.l = l;
        ch.r = to_d(tk[0]);
        ch.nproj = to_i(tk[1]);
        ch.h.assign(static_cast<std::size_t>(ch.nproj) * ch.nproj, 0.0);

        // Upper triangle, row by row. Row 0's entries follow (r, nproj) on the
        // same line; rows 1..nproj-1 each occupy their own continuation line.
        for (int i = 0; i < ch.nproj; ++i) {
            const std::vector<std::string>* row = &tk;
            std::vector<std::string> cont;
            int base = 0;
            if (i == 0) {
                base = 2; // skip r_l and nproj on the first line
            } else {
                if (!next_content_line(in, line)) fail("unexpected EOF reading h-matrix row");
                cont = tokens(line);
                row = &cont;
                base = 0;
            }
            int expected = ch.nproj - i; // entries h(i,i)..h(i,nproj-1)
            if (static_cast<int>(row->size()) - base < expected)
                fail("h-matrix row " + std::to_string(i) + " short: '" + line + "'");
            for (int j = 0; j < expected; ++j) {
                double v = to_d((*row)[base + j]);
                int col = i + j;
                ch.h[static_cast<std::size_t>(i) * ch.nproj + col] = v;
                ch.h[static_cast<std::size_t>(col) * ch.nproj + i] = v; // symmetric
            }
        }
        pot.channels.push_back(std::move(ch));
    }

    return pot;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
