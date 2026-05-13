/// @file cpu/pw/internal/upf_reader.cpp
/// @brief UPF v1 XML-like format parser.

#include "atomiverse/cpu/pw/internal/upf_reader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstring>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

namespace {

std::string extract_attr(const std::string& line, const std::string& attr) {
    auto p = line.find(attr + "=\"");
    if (p == std::string::npos) return "";
    p += attr.size() + 2; // skip past attr="
    auto q = line.find("\"", p);
    if (q == std::string::npos) return "";
    return line.substr(p, q - p);
}

bool tag_open(const std::string& line, const std::string& tag) {
    return line.find("<" + tag) != std::string::npos && line.find("/>") == std::string::npos;
}
bool tag_close(const std::string& line, const std::string& tag) {
    return line.find("</" + tag) != std::string::npos;
}

std::vector<Real> parse_values(std::istream& in, const std::string& end_tag) {
    std::vector<Real> vals;
    std::string line;
    while (std::getline(in, line)) {
        if (line.find(end_tag) != std::string::npos) break;
        // trim
        auto s = line;
        while (!s.empty() && (s[0]==' '||s[0]=='\t')) s.erase(0,1);
        if (s.empty()) continue;
        std::istringstream ss(s);
        double d;
        while (ss >> d) vals.push_back(static_cast<Real>(d));
    }
    return vals;
}

} // anonymous

UPFData read_upf_v1(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("UPF: cannot open " + path);

    UPFData upf;
    std::string line;

    while (std::getline(f, line)) {
        // PP_HEADER
        if (line.find("<PP_HEADER") != std::string::npos) {
            auto zv = extract_attr(line, "z_valence");
            auto lm = extract_attr(line, "l_max");
            if (!zv.empty()) upf.z_val = std::stoi(zv);
            if (!lm.empty()) upf.l_max = std::stoi(lm);
        }
        // PP_MESH
        if (tag_open(line, "PP_MESH")) {
            // read PP_R
            while (std::getline(f, line)) {
                if (tag_open(line, "PP_R")) {
                    upf.r = parse_values(f, "PP_R");
                    break;
                }
                if (tag_close(line, "PP_MESH")) break;
            }
        }
        // PP_LOCAL
        if (tag_open(line, "PP_LOCAL")) {
            upf.v_loc = parse_values(f, "PP_LOCAL");
        }
        // PP_NLCC
        if (tag_open(line, "PP_NLCC")) {
            upf.rho_atom = parse_values(f, "PP_NLCC");
        }
        // PP_NONLOCAL
        if (tag_open(line, "PP_NONLOCAL")) {
            while (std::getline(f, line)) {
                if (tag_close(line, "PP_NONLOCAL")) break;
                if (tag_open(line, "PP_BETA")) {
                    int l = 0;
                    auto lang = extract_attr(line, "angular_momentum");
                    if (!lang.empty()) l = std::stoi(lang);
                    auto beta_vals = parse_values(f, "PP_BETA");
                    // Read PP_DIJ for this l
                    Real d = 0;
                    auto pos = f.tellg();
                    std::string next_line;
                    if (std::getline(f, next_line)) {
                        if (tag_open(next_line, "PP_DIJ")) {
                            auto dij_vals = parse_values(f, "PP_DIJ");
                            if (!dij_vals.empty()) d = dij_vals[0];
                        } else {
                            f.seekg(pos);
                        }
                    }
                    UPFData::KBProjector kb;
                    kb.l = l;
                    kb.beta = std::move(beta_vals);
                    kb.d = d;
                    upf.kb.push_back(std::move(kb));
                }
            }
        }
    }

    if (upf.r.empty())
        throw std::runtime_error("UPF: no mesh found in " + path);
    return upf;
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
