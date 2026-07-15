/// @file cpu/comm/bcast_frame.cpp
/// @brief Whole-Frame broadcast over the registry-driven exchange codec.

#include "atomiverse/cpu/comm/bcast_frame.h"

#include "atomiverse/cpu/frame/frame.h"

#ifdef ATV_HAS_MPI

#include <mpi.h>

#include <climits>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "atomiverse/core/precision.h"
#include "atomiverse/cpu/comm/atom_pack.h"

namespace ATV {
namespace cpu {
namespace comm {

namespace {

/// Fixed-size wire header, broadcast as raw bytes ahead of the row payload.
/// `has_bonded` travels explicitly because the receiving Atom's bonded flag
/// is only ever set by `ensure_bonded_storage()` -- the defensive stride
/// growth inside `unpack_exchange` does not flip it.
struct BcastHeader {
    BigInt n_atoms = 0;
    BigInt payload_doubles = 0;
    int has_bonded = 0;
    Box box{};
    SystemScalars scalars{};
};

// The header travels as MPI_BYTE. `cpu::Box` is not formally trivially
// copyable (its Mat3 members hold xtensor_fixed storage), but its byte
// transit over MPI_BYTE is the established precedent — atom_migration.cpp
// broadcasts `frame.box` the same way (inline storage, no indirection).
// Assert the parts with no such precedent.
static_assert(std::is_trivially_copyable<SystemScalars>::value,
              "SystemScalars is broadcast as MPI_BYTE");

} // namespace

void bcast_frame(Frame& frame, int root) {
    int inited = 0;
    MPI_Initialized(&inited);
    if (!inited)
        return;
    int rank = 0, world = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    if (world <= 1)
        return;

    BcastHeader hdr{};
    std::vector<double> payload;
    if (rank == root) {
        const int n = frame.atoms.nlocal;
        payload.resize(static_cast<std::size_t>(n) *
                       static_cast<std::size_t>(
                           max_exchange_doubles(frame.atoms)));
        long long m = 0;
        for (int i = 0; i < n; ++i)
            m += pack_exchange(frame.atoms, i, payload.data() + m);
        hdr.n_atoms = n;
        hdr.payload_doubles = m;
        hdr.has_bonded = frame.atoms.has_bonded() ? 1 : 0;
        hdr.box = frame.box;
        hdr.scalars = frame.scalars;
    }
    MPI_Bcast(&hdr, sizeof(BcastHeader), MPI_BYTE, root, MPI_COMM_WORLD);

    // Checked after the header collective so every rank throws together
    // (a pre-broadcast root throw would deadlock the others in MPI_Bcast).
    if (hdr.payload_doubles > INT_MAX)
        throw std::runtime_error(
            "bcast_frame: payload of " + std::to_string(hdr.payload_doubles) +
            " doubles exceeds the MPI_Bcast count limit (INT_MAX); no "
            "chunked fallback by contract");

    if (rank == root) {
        if (hdr.payload_doubles > 0)
            MPI_Bcast(payload.data(), static_cast<int>(hdr.payload_doubles),
                      MPI_DOUBLE, root, MPI_COMM_WORLD);
        return; // the root frame is never modified
    }

    // Non-root: receive with zero-filled slack so a divergent (wider) local
    // registration over-consumes into zeroed memory instead of out of
    // bounds; the count check below then rejects it.
    const int n = static_cast<int>(hdr.n_atoms);
    if (hdr.has_bonded)
        frame.atoms.ensure_bonded_storage();
    if (n > 0)
        frame.atoms.set_capacity(n);
    std::vector<double> buf(
        static_cast<std::size_t>(hdr.payload_doubles) +
            static_cast<std::size_t>(max_exchange_doubles(frame.atoms)),
        0.0);
    if (hdr.payload_doubles > 0)
        MPI_Bcast(buf.data(), static_cast<int>(hdr.payload_doubles),
                  MPI_DOUBLE, root, MPI_COMM_WORLD);

    long long m = 0;
    for (int slot = 0; slot < n; ++slot) {
        m += unpack_exchange(frame.atoms, slot, buf.data() + m);
        if (m > hdr.payload_doubles)
            break; // divergence detected; the check below names it
    }
    if (m != hdr.payload_doubles)
        throw std::runtime_error(
            "bcast_frame: decode consumed " + std::to_string(m) +
            " doubles but the root announced " +
            std::to_string(hdr.payload_doubles) +
            " -- exchanged extension-column registration differs between "
            "this rank and the root");

    frame.atoms.nlocal = n;
    frame.atoms.nghost = 0;
    frame.atoms.map_clear();
    frame.atoms.map_set(0, n);
    frame.box = hdr.box;
    frame.scalars = hdr.scalars;
}

} // namespace comm
} // namespace cpu
} // namespace ATV

#else // !ATV_HAS_MPI

namespace ATV {
namespace cpu {
namespace comm {

void bcast_frame(Frame& /*frame*/, int /*root*/) {}

} // namespace comm
} // namespace cpu
} // namespace ATV

#endif // ATV_HAS_MPI
