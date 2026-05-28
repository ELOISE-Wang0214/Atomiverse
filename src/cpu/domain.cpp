/// @file src/cpu/domain.cpp
/// @brief cpu::Domain -- impl. Cartesian-topology setup + slab math.

#include "atomiverse/cpu/domain.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>

#include "atomiverse/cpu/box.h"

namespace ATV {
namespace cpu {

void Domain::compute_sublo_subhi(const Box& box) {
    for (int d = 0; d < 3; ++d) {
        const double L = box.length(d);
        if (static_cast<int>(cuts[d].size()) == dims[d] + 1) {
            // Non-uniform cut planes (set by RCBBalancer). cuts[d] is in
            // box-length fractions, so the partition survives NPT/Deform
            // box rescaling without recomputation.
            sublo[d] = cuts[d][static_cast<std::size_t>(coords[d])] * L;
            subhi[d] = cuts[d][static_cast<std::size_t>(coords[d] + 1)] * L;
        } else {
            const double slab = L / static_cast<double>(dims[d]);
            sublo[d] = slab * static_cast<double>(coords[d]);
            subhi[d] = sublo[d] + slab;
        }
    }
}

#ifdef ATV_HAS_MPI

void Domain::release_comm(bool /*mpi_is_finalizing*/) {
    // We do not own a Cartesian sub-comm in this implementation (cart_comm
    // aliases the world communicator, see init_topology). Just clear the
    // handle.
    cart_comm = MPI_COMM_NULL;
}

void Domain::init_topology(MPI_Comm world) {
    int world_size = 0;
    int my_rank = 0;
    MPI_Comm_size(world, &world_size);
    MPI_Comm_rank(world, &my_rank);

    if (world_size <= 1) {
        // Single rank: keep defaults (dims={1,1,1}, no cart_comm).
        return;
    }

    // Do NOT use MPI_Cart_create -- on OpenMPI 5.0.x (macOS) the resulting
    // sub-communicator triggers a segfault inside MPI_Finalize when any
    // user MPI_Comm_free has been issued, or when the same process
    // repeatedly creates Cartesian comms (e.g. across GoogleTest cases).
    // Compute the brick decomposition directly on top of the world
    // communicator using simple rank arithmetic; this is functionally
    // equivalent and avoids the OpenMPI bug entirely.
    cart_comm = world;

    // `MPI_Dims_create` is safe -- it's a pure algorithm with no comm state.
    // Optional env-var override `ATV_FORCE_DIMS=Nx,Ny,Nz` skips the
    // automatic 3D decomposition and uses the given factorization.
    dims[0] = dims[1] = dims[2] = 0;
    if (const char* env = std::getenv("ATV_FORCE_DIMS")) {
        int dx, dy, dz;
        if (std::sscanf(env, "%d,%d,%d", &dx, &dy, &dz) == 3
            && dx * dy * dz == world_size) {
            dims[0] = dx;
            dims[1] = dy;
            dims[2] = dz;
        }
    }
    if (dims[0] == 0) {
        MPI_Dims_create(world_size, 3, dims);
    }

    // Convention: rank index = ix * (dy * dz) + iy * dz + iz.
    coords[0] = my_rank / (dims[1] * dims[2]);
    coords[1] = (my_rank / dims[2]) % dims[1];
    coords[2] = my_rank % dims[2];

    auto rank_at = [&](int cx, int cy, int cz) {
        cx = ((cx % dims[0]) + dims[0]) % dims[0];
        cy = ((cy % dims[1]) + dims[1]) % dims[1];
        cz = ((cz % dims[2]) + dims[2]) % dims[2];
        return cx * (dims[1] * dims[2]) + cy * dims[2] + cz;
    };

    // Face order: -x, +x, -y, +y, -z, +z.
    face_neighbors[0] = rank_at(coords[0] - 1, coords[1], coords[2]);
    face_neighbors[1] = rank_at(coords[0] + 1, coords[1], coords[2]);
    face_neighbors[2] = rank_at(coords[0], coords[1] - 1, coords[2]);
    face_neighbors[3] = rank_at(coords[0], coords[1] + 1, coords[2]);
    face_neighbors[4] = rank_at(coords[0], coords[1], coords[2] - 1);
    face_neighbors[5] = rank_at(coords[0], coords[1], coords[2] + 1);
}

#endif // ATV_HAS_MPI

} // namespace cpu
} // namespace ATV
