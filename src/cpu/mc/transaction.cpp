/// @file src/cpu/mc/transaction.cpp
/// @brief cpu::mc::Transaction implementation — Frame-patch record/rollback
///        over the migration row serializers.

#include "atomiverse/cpu/mc/transaction.h"

#include <stdexcept>

#include "atomiverse/cpu/comm/atom_pack.h"

namespace ATV {
namespace cpu {
namespace mc {

void Transaction::begin(Verb v, const cpu::Frame& f) {
    clear();
    verb = v;
    nlocal_before_ = f.atoms.nlocal;
}

void Transaction::record_touch(const cpu::Frame& f, int row) {
    const int k = static_cast<int>(touched_.size());
    if (patch_.atoms.nmax < k + 1)
        patch_.atoms.grow(k + 1);
    buf_.resize(static_cast<std::size_t>(comm::max_exchange_doubles(f.atoms)));
    comm::pack_exchange(f.atoms, row, buf_.data());
    comm::unpack_exchange(patch_.atoms, k, buf_.data());
    touched_.push_back(static_cast<LocalInt>(row));
    patch_.atoms.nlocal = k + 1;
}

void Transaction::record_box(const cpu::Frame& f) {
    patch_.box = f.box;
    box_recorded_ = true;
}

void Transaction::rollback(cpu::Frame& f) const {
    if (nlocal_before_ < 0)
        throw std::logic_error(
            "mc::Transaction::rollback: no proposal recorded (begin() never "
            "ran)");
    // Population first (an Insert's appended rows fall off the end), then the
    // touched rows (a Delete's victim gets its original content back — the
    // compaction survivor copied into its slot is overwritten, and the
    // survivor's own end-row content was never destroyed by copy_row).
    f.atoms.nlocal = nlocal_before_;
    for (std::size_t k = 0; k < touched_.size(); ++k) {
        buf_.resize(
            static_cast<std::size_t>(comm::max_exchange_doubles(patch_.atoms)));
        comm::pack_exchange(patch_.atoms, static_cast<int>(k), buf_.data());
        comm::unpack_exchange(f.atoms, touched_[k], buf_.data());
    }
    if (box_recorded_)
        f.box = patch_.box;
}

void Transaction::clear() noexcept {
    verb = Verb::None;
    touched_.clear();
    patch_.atoms.nlocal = 0;
    nlocal_before_ = -1;
    box_recorded_ = false;
}

} // namespace mc
} // namespace cpu
} // namespace ATV
