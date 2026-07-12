/// @file cpu/replica/replica_expansion.cpp
/// @brief LinearInterpolation / IdentityReplication — the two shipped expansions.

#include "atomiverse/cpu/replica/replica_expansion.h"

#include <stdexcept>
#include <string>

namespace ATV {
namespace cpu {

LinearInterpolation::LinearInterpolation(cpu::Frame final_state, int n_images)
    : final_(std::move(final_state)), n_images_(n_images) {
    if (n_images < 3) {
        throw std::invalid_argument(
            "LinearInterpolation: n_images must be >= 3 (got "
            + std::to_string(n_images)
            + ") — with two images there is nothing between the endpoints to relax, "
              "which is not a path calculation");
    }
}

std::vector<Replica> LinearInterpolation::expand(const cpu::Frame& initial) const {
    const int n = initial.atoms.size();
    if (final_.atoms.size() != n) {
        throw std::invalid_argument(
            "LinearInterpolation: the final state has "
            + std::to_string(final_.atoms.size()) + " atoms but the initial state has "
            + std::to_string(n) + " — a band between different systems is not a band");
    }

    std::vector<Replica> replicas;
    replicas.reserve(n_images_);
    for (int im = 0; im < n_images_; ++im) {
        // t = 0 at the initial state, 1 at the final state. This is also the replica's
        // parameter: for a band, the thing that distinguishes an image from its siblings
        // IS its position along the path.
        const double t = double(im) / double(n_images_ - 1);

        Replica r{initial, t}; // copies topology, masses, types, box — all shared
        for (int i = 0; i < n; ++i) {
            r.frame.atoms.x[i] = (1.0 - t) * initial.atoms.x[i] + t * final_.atoms.x[i];
            r.frame.atoms.y[i] = (1.0 - t) * initial.atoms.y[i] + t * final_.atoms.y[i];
            r.frame.atoms.z[i] = (1.0 - t) * initial.atoms.z[i] + t * final_.atoms.z[i];
        }
        replicas.push_back(std::move(r));
    }
    return replicas;
}

void LinearInterpolation::collapse(const std::vector<Replica>& replicas,
                                   cpu::Frame& out) const {
    if (replicas.empty())
        throw std::logic_error("LinearInterpolation::collapse: no replicas");

    int best = 0;
    for (int i = 1; i < static_cast<int>(replicas.size()); ++i) {
        if (replicas[i].frame.scalars.energy > replicas[best].frame.scalars.energy)
            best = i;
    }
    const cpu::Frame& saddle = replicas[best].frame;
    for (int i = 0; i < out.atoms.size(); ++i) {
        out.atoms.x[i] = saddle.atoms.x[i];
        out.atoms.y[i] = saddle.atoms.y[i];
        out.atoms.z[i] = saddle.atoms.z[i];
    }
    out.scalars = saddle.scalars;
}

IdentityReplication::IdentityReplication(std::vector<double> parameters)
    : parameters_(std::move(parameters)) {
    if (parameters_.size() < 2) {
        throw std::invalid_argument(
            "IdentityReplication: needs >= 2 replicas (got "
            + std::to_string(parameters_.size())
            + ") — one replica is not a replica set, it is a Frame");
    }
}

std::vector<Replica> IdentityReplication::expand(const cpu::Frame& initial) const {
    std::vector<Replica> replicas;
    replicas.reserve(parameters_.size());
    for (const double p : parameters_)
        replicas.push_back(Replica{initial, p});
    return replicas;
}

void IdentityReplication::collapse(const std::vector<Replica>& replicas,
                                   cpu::Frame& out) const {
    if (replicas.empty())
        throw std::logic_error("IdentityReplication::collapse: no replicas");
    const cpu::Frame& r0 = replicas[0].frame;
    for (int i = 0; i < out.atoms.size(); ++i) {
        out.atoms.x[i] = r0.atoms.x[i];
        out.atoms.y[i] = r0.atoms.y[i];
        out.atoms.z[i] = r0.atoms.z[i];
    }
    out.scalars = r0.scalars;
}

} // namespace cpu
} // namespace ATV
