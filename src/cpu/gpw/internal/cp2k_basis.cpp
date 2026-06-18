/// @file src/cpu/gpw/internal/cp2k_basis.cpp
/// @brief Implementation of CP2K Gaussian basis-set file parsing (see header).

#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ATV {
namespace cpu {
namespace gpw {

namespace {

/// @brief Read the next content line, skipping blank and '#'-comment lines.
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
    throw std::runtime_error("parse_gth_basis: " + msg);
}

double to_d(const std::string& s) { return std::stod(s); }
int to_i(const std::string& s) { return std::stoi(s); }

} // namespace

GaussianBasis parse_gth_basis(const std::string& path, const std::string& element,
                              const std::string& name) {
    std::ifstream in(path);
    if (!in) fail("cannot open file '" + path + "'");

    std::string line;
    // Locate the header: first token == element AND `name` among the alias
    // tokens. Data lines start with a number, so they never match an element.
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

    GaussianBasis basis;
    basis.element = element;
    basis.name = name;

    // Number of contraction sets.
    if (!next_content_line(in, line)) fail("unexpected EOF reading set count");
    int nset = to_i(tokens(line).at(0));

    for (int s = 0; s < nset; ++s) {
        if (!next_content_line(in, line)) fail("unexpected EOF reading set header");
        auto tk = tokens(line);
        // Header: n lmin lmax nexp nshell_lmin ... nshell_lmax
        if (tk.size() < 4) fail("malformed set header: '" + line + "'");
        GaussianBasisSet set;
        set.n = to_i(tk[0]);
        set.lmin = to_i(tk[1]);
        set.lmax = to_i(tk[2]);
        int nexp = to_i(tk[3]);
        int nl = set.lmax - set.lmin + 1;
        if (nl <= 0) fail("invalid lmin/lmax in set header: '" + line + "'");
        if (static_cast<int>(tk.size()) < 4 + nl)
            fail("set header missing per-l shell counts: '" + line + "'");
        set.ncoef = 0;
        for (int i = 0; i < nl; ++i) {
            int ns = to_i(tk[4 + i]);
            set.nshell.push_back(ns);
            set.ncoef += ns;
        }

        // nexp rows: <exp> <ncoef coefficients>.
        set.exps.reserve(nexp);
        set.coef.reserve(static_cast<std::size_t>(nexp) * set.ncoef);
        for (int e = 0; e < nexp; ++e) {
            if (!next_content_line(in, line)) fail("unexpected EOF reading exponent row");
            auto row = tokens(line);
            if (static_cast<int>(row.size()) < 1 + set.ncoef)
                fail("exponent row has too few coefficients: '" + line + "'");
            set.exps.push_back(to_d(row[0]));
            for (int c = 0; c < set.ncoef; ++c) set.coef.push_back(to_d(row[1 + c]));
        }
        basis.sets.push_back(std::move(set));
    }

    return basis;
}

} // namespace gpw
} // namespace cpu
} // namespace ATV
