/// @file cpu/md/cell_list.cpp
/// @brief CellList implementation -- LAMMPS-style linked-cell neighbor list.

#include "atomiverse/cpu/md/cell_list.h"
#include "atomiverse/core/log.h"
#include "atomiverse/core/mat3.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

// ── CellList ─────────────────────────────────────────────────────

CellList::CellList(double cutoff, double skin)
    : cutoff_(cutoff),
      skin_(skin < 0 ? double(0.3) * cutoff : skin),
      cutoff_skin_(cutoff + (skin < 0 ? double(0.3) * cutoff : skin)) {}

void CellList::init(const cpu::Frame& frame, double /*cutoff_skin*/) {
    cutoff_skin_ = cutoff_ + skin_;
    // Single-rank: no axis extension; whole-box grid; ghosts wrap into [0, L).
    for (int d = 0; d < 3; ++d) {
        axis_extended_[d] = false;
        origin_[d] = 0.0;
    }
    capacity_ = frame.atoms.size();
    compute_grid(frame);
    build_stencil();

    atom_cell_ = xt::zeros<int>({capacity_});
    binned_atoms_ = xt::zeros<int>({capacity_});
    x_save_ = xt::zeros<double>({capacity_});
    y_save_ = xt::zeros<double>({capacity_});
    z_save_ = xt::zeros<double>({capacity_});
}

void CellList::init_extended(const cpu::Frame& frame, double /*cutoff_skin*/,
                              const int dims[3], int capacity) {
    cutoff_skin_ = cutoff_ + skin_;
    // Multi-rank: extend each axis with dims>1 to cover ghost halo.
    for (int d = 0; d < 3; ++d) {
        axis_extended_[d] = (dims[d] > 1);
        origin_[d] = axis_extended_[d] ? -cutoff_skin_ : 0.0;
    }
    capacity_ = capacity;
    compute_grid(frame);
    build_stencil();

    atom_cell_ = xt::zeros<int>({capacity_});
    binned_atoms_ = xt::zeros<int>({capacity_});
    x_save_ = xt::zeros<double>({capacity_});
    y_save_ = xt::zeros<double>({capacity_});
    z_save_ = xt::zeros<double>({capacity_});
}

void CellList::compute_grid(const cpu::Frame& frame) {
    // Per-axis effective grid length: box edge plus halo skin on each
    // extended side. Single-rank axes use the box edge as-is.
    double lx_eff = frame.box.length(0) + (axis_extended_[0] ? double(2) * cutoff_skin_ : double(0));
    double ly_eff = frame.box.length(1) + (axis_extended_[1] ? double(2) * cutoff_skin_ : double(0));
    double lz_eff = frame.box.length(2) + (axis_extended_[2] ? double(2) * cutoff_skin_ : double(0));

    ncx_ = std::max(int(1), static_cast<int>(std::floor(lx_eff / cutoff_skin_)));
    ncy_ = std::max(int(1), static_cast<int>(std::floor(ly_eff / cutoff_skin_)));
    ncz_ = std::max(int(1), static_cast<int>(std::floor(lz_eff / cutoff_skin_)));

    csx_ = lx_eff / ncx_;
    csy_ = ly_eff / ncy_;
    csz_ = lz_eff / ncz_;

    total_cells_ = ncx_ * ncy_ * ncz_;
    cell_count_ = xt::zeros<int>({total_cells_});
    cell_offset_ = xt::zeros<int>({total_cells_});
}

void CellList::build_stencil() {
    half_stencil_.clear();
    full_stencil_.clear();

    for (int di = -1; di <= 1; ++di) {
        for (int dj = -1; dj <= 1; ++dj) {
            for (int dk = -1; dk <= 1; ++dk) {
                full_stencil_.push_back({di, dj, dk});
                // Half stencil: include self (0,0,0), and any neighbor
                // whose flat offset > 0 in the original integer ordering
                // (lexicographic on (di, dj, dk)). The full-list filter
                // catches every pair twice; the half-list filter is the
                // same pair pruned to one direction.
                const int flat = di * ncy_ * ncz_ + dj * ncz_ + dk;
                if (flat > 0 || (di == 0 && dj == 0 && dk == 0))
                    half_stencil_.push_back({di, dj, dk});
            }
        }
    }
}

void CellList::bin_atoms(const cpu::Frame& frame, int n_active) {
    int n = (n_active < 0) ? frame.atoms.size() : n_active;
    double lx = frame.box.length(0);
    double ly = frame.box.length(1);
    double lz = frame.box.length(2);

    // Clear cell counts
    std::fill(cell_count_.begin(), cell_count_.end(), int(0));

    // Assign atoms to cells. Single-rank axes wrap into [0, L); extended
    // axes leave the position alone so PBC-shifted ghost coordinates land
    // in the dedicated halo cells around the slab.
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        double xi = frame.atoms.x[i];
        double yi = frame.atoms.y[i];
        double zi = frame.atoms.z[i];
        if (!axis_extended_[0]) {
            xi -= lx * std::floor(xi / lx);
        }
        if (!axis_extended_[1]) {
            yi -= ly * std::floor(yi / ly);
        }
        if (!axis_extended_[2]) {
            zi -= lz * std::floor(zi / lz);
        }

        int cx = static_cast<int>(std::floor((xi - origin_[0]) / csx_));
        int cy = static_cast<int>(std::floor((yi - origin_[1]) / csy_));
        int cz = static_cast<int>(std::floor((zi - origin_[2]) / csz_));
        cx = std::max(0, std::min(cx, ncx_ - 1));
        cy = std::max(0, std::min(cy, ncy_ - 1));
        cz = std::max(0, std::min(cz, ncz_ - 1));

        atom_cell_[i] = cx * ncy_ * ncz_ + cy * ncz_ + cz;
    }

    // Count atoms per cell (serial -- fast, O(n_cells) is small)
    for (int i = 0; i < n; ++i)
        cell_count_[atom_cell_[i]]++;

    // Prefix sum
    cell_offset_[0] = 0;
    for (int c = 1; c < total_cells_; ++c)
        cell_offset_[c] = cell_offset_[c - 1] + cell_count_[c - 1];

    // Scatter atoms into bins (use temporary write pointers)
    auto write_ptr = cell_offset_;  // copy
    for (int i = 0; i < n; ++i) {
        int c = atom_cell_[i];
        binned_atoms_[write_ptr[c]] = i;
        write_ptr[c]++;
    }
}

void CellList::save_positions(const cpu::Frame& frame, int n_active) {
    int n = (n_active < 0) ? frame.atoms.size() : n_active;
    for (int i = 0; i < n; ++i) {
        x_save_[i] = frame.atoms.x[i];
        y_save_[i] = frame.atoms.y[i];
        z_save_[i] = frame.atoms.z[i];
    }
}

void CellList::build(const cpu::Frame& frame, int n_active) {
    bin_atoms(frame, n_active);
    save_positions(frame, n_active);
    n_at_last_build_ = (n_active < 0) ? frame.atoms.size() : n_active;
    box_at_last_build_ = frame.box.length(0);
}

bool CellList::needs_rebuild(const cpu::Frame& frame, int n_active) const {
    int n = (n_active < 0) ? frame.atoms.size() : n_active;

    // Multi-rank: AtomMigration / HaloExchange may change n_local + n_ghost
    // between launches.  When that happens, `binned_atoms_` from the
    // previous build still references the OLD atom indices, including ones
    // beyond the current `n_active`. Reading those produces an out-of-
    // bounds access in `build_half`.  Force a rebuild whenever the active
    // count moves.
    if (n != n_at_last_build_) return true;

    // NPT: barostat may have scaled the box, invalidating the cell grid
    // (`csx_/csy_/csz_` were sized for the old box edge). Force a rebuild.
    if (frame.box.length(0) != box_at_last_build_) return true;

    double max_disp_sq = 0;

    for (int i = 0; i < n; ++i) {
        double dx = frame.atoms.x[i] - x_save_[i];
        double dy = frame.atoms.y[i] - y_save_[i];
        double dz = frame.atoms.z[i] - z_save_[i];
        // Use minimum-image convention so atom wrapping across PBC
        // boundaries doesn't falsely signal a ~box-length displacement.
        frame.box.min_image(dx, dy, dz);
        double dsq = dx * dx + dy * dy + dz * dz;
        if (dsq > max_disp_sq) max_disp_sq = dsq;
    }

    // Rebuild when 2 * max_disp > skin
    return max_disp_sq * double(4) > skin_ * skin_;
}

void CellList::update(const cpu::Frame& frame, int n_active) {
    if (needs_rebuild(frame, n_active))
        build(frame, n_active);
}

namespace {

/// @brief Grow-never-shrink resize of a std::vector<int>. When
/// `n <= current size` this is a no-op; when `n > current size` it
/// hands off to `std::vector::resize`, which amortises allocation and
/// preserves existing contents. Used by `refresh_*` so the per-atom
/// CSR arrays keep a monotonically-non-decreasing capacity across
/// steady-state launches.
inline void ensure_at_least(std::vector<int>& arr, std::size_t n) {
    if (arr.size() < n) arr.resize(n);
}

} // namespace

void CellList::refresh_half(const cpu::Frame& frame, int n_active,
                            PairList& result) const {
    const int n = (n_active < 0) ? static_cast<int>(frame.atoms.size())
                                 : n_active;
    const double rc_sq = cutoff_skin_ * cutoff_skin_;

    ensure_at_least(result.ilist, static_cast<std::size_t>(n));
    ensure_at_least(result.numneigh, static_cast<std::size_t>(n));
    ensure_at_least(result.firstneigh, static_cast<std::size_t>(n));

    for (int i = 0; i < n; ++i) result.ilist[i] = i;
    result.inum = n;

    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
#endif

    // ── Serial fast path ────────────────────────────────────────────
    // Single-thread builds write the CSR directly in one pass: for each
    // atom i we reserve its slot in `neighbors[]`, iterate the stencil,
    // append each surviving j to a running write cursor, then commit
    // `numneigh[i]`. No per-thread pair buffers, no pair-count pass, no
    // prefix-sum pass, no scatter pass — the old path made 3 trips over
    // ~250k pairs plus an `std::vector<std::vector<pair<int,int>>>`
    // allocation every call. `neighbors` grows amortized O(1) the first
    // few times and then stays at capacity (~250k ints = 1 MB for 10k
    // atoms). `firstneigh[i]` is the live cursor at function entry, so
    // consumers that later index `neighbors[firstneigh[i] + k]` see the
    // exact slot we wrote. Matches LAMMPS `NPairHalfBin::build`.
    if (n_threads == 1) {
        result.neighbors.clear();

        const double* __restrict x_ptr = frame.atoms.x.data();
        const double* __restrict y_ptr = frame.atoms.y.data();
        const double* __restrict z_ptr = frame.atoms.z.data();
        const double lx = frame.box.length(0);
        const double ly = frame.box.length(1);
        const double lz = frame.box.length(2);
        const double half_lx = 0.5 * lx;
        const double half_ly = 0.5 * ly;
        const double half_lz = 0.5 * lz;
        const int* __restrict atom_cell_ptr = atom_cell_.data();
        const int* __restrict cell_count_ptr = cell_count_.data();
        const int* __restrict cell_offset_ptr = cell_offset_.data();
        const int* __restrict binned_ptr = binned_atoms_.data();

        int visited[27];
        int cursor = 0;

        for (int i = 0; i < n; ++i) {
            result.firstneigh[i] = cursor;
            const int ci = atom_cell_ptr[i];
            const int cz0 = ci % ncz_;
            const int cy0 = (ci / ncz_) % ncy_;
            const int cx0 = ci / (ncy_ * ncz_);
            const double xi = x_ptr[i];
            const double yi = y_ptr[i];
            const double zi = z_ptr[i];

            int visited_n = 0;

            for (auto s : full_stencil_) {
                int cx = cx0 + s.di;
                int cy = cy0 + s.dj;
                int cz = cz0 + s.dk;

                if (axis_extended_[0]) {
                    if (cx < 0 || cx >= ncx_) continue;
                } else {
                    if (cx <  0)    cx += ncx_;
                    if (cx >= ncx_) cx -= ncx_;
                }
                if (axis_extended_[1]) {
                    if (cy < 0 || cy >= ncy_) continue;
                } else {
                    if (cy <  0)    cy += ncy_;
                    if (cy >= ncy_) cy -= ncy_;
                }
                if (axis_extended_[2]) {
                    if (cz < 0 || cz >= ncz_) continue;
                } else {
                    if (cz <  0)    cz += ncz_;
                    if (cz >= ncz_) cz -= ncz_;
                }

                const int cj = cx * ncy_ * ncz_ + cy * ncz_ + cz;

                // Small-grid PBC-wrap dedup. Stored as a fixed-size
                // int[27] so the compiler keeps it in a register file.
                bool seen = false;
                for (int v = 0; v < visited_n; ++v) {
                    if (visited[v] == cj) { seen = true; break; }
                }
                if (seen) continue;
                visited[visited_n++] = cj;

                const int start = cell_offset_ptr[cj];
                const int end   = start + cell_count_ptr[cj];

                for (int idx = start; idx < end; ++idx) {
                    int j = binned_ptr[idx];
                    if (j <= i) continue;

                    double dx = xi - x_ptr[j];
                    double dy = yi - y_ptr[j];
                    double dz = zi - z_ptr[j];
                    if (dx >  half_lx) dx -= lx;
                    if (dx < -half_lx) dx += lx;
                    if (dy >  half_ly) dy -= ly;
                    if (dy < -half_ly) dy += ly;
                    if (dz >  half_lz) dz -= lz;
                    if (dz < -half_lz) dz += lz;

                    double r2 = dx * dx + dy * dy + dz * dz;
                    if (r2 < rc_sq) {
                        // push_back on an empty-then-growing vector
                        // amortizes to O(1); after a few rebuilds
                        // capacity saturates and push_back is a
                        // pointer-bump store.
                        result.neighbors.push_back(j);
                        ++cursor;
                    }
                }
            }
            result.numneigh[i] = cursor - result.firstneigh[i];
        }
        result.total_pairs = cursor;
        return;
    }

    // ── Multi-thread path ───────────────────────────────────────────
    // Thread-private pair buffers feed a pair-count + prefix-sum +
    // scatter finalizer. Slightly more work per pair than the serial
    // fast path but keeps scaling on OMP ≥ 2. The old single-path
    // implementation paid this overhead even with OMP=1 — that cost is
    // what the serial branch above reclaims.
    std::vector<std::vector<std::pair<int, int>>> thread_pairs(n_threads);
    for (int i = 0; i < n; ++i) result.numneigh[i] = 0;

    #pragma omp parallel
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        auto& my_pairs = thread_pairs[tid];

        std::vector<int> visited;
        visited.reserve(27);

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            const int ci = atom_cell_[i];
            const int cz0 = ci % ncz_;
            const int cy0 = (ci / ncz_) % ncy_;
            const int cx0 = ci / (ncy_ * ncz_);

            visited.clear();

            for (auto s : full_stencil_) {
                int cx = cx0 + s.di;
                int cy = cy0 + s.dj;
                int cz = cz0 + s.dk;

                if (axis_extended_[0]) {
                    if (cx < 0 || cx >= ncx_) continue;
                } else {
                    if (cx <  0)    cx += ncx_;
                    if (cx >= ncx_) cx -= ncx_;
                }
                if (axis_extended_[1]) {
                    if (cy < 0 || cy >= ncy_) continue;
                } else {
                    if (cy <  0)    cy += ncy_;
                    if (cy >= ncy_) cy -= ncy_;
                }
                if (axis_extended_[2]) {
                    if (cz < 0 || cz >= ncz_) continue;
                } else {
                    if (cz <  0)    cz += ncz_;
                    if (cz >= ncz_) cz -= ncz_;
                }

                const int cj = cx * ncy_ * ncz_ + cy * ncz_ + cz;

                bool seen = false;
                for (int v : visited) if (v == cj) { seen = true; break; }
                if (seen) continue;
                visited.push_back(cj);

                const int start = cell_offset_[cj];
                const int end   = start + cell_count_[cj];

                for (int idx = start; idx < end; ++idx) {
                    int j = binned_atoms_[idx];
                    if (j <= i) continue;

                    double dx = frame.atoms.x[i] - frame.atoms.x[j];
                    double dy = frame.atoms.y[i] - frame.atoms.y[j];
                    double dz = frame.atoms.z[i] - frame.atoms.z[j];
                    frame.box.min_image(dx, dy, dz);

                    double r2 = dx * dx + dy * dy + dz * dz;
                    if (r2 < rc_sq) {
                        my_pairs.emplace_back(i, j);
                    }
                }
            }
        }
    }

    for (auto& tp : thread_pairs)
        for (auto& [i, j] : tp)
            result.numneigh[i]++;

    result.firstneigh[0] = 0;
    for (int i = 1; i < n; ++i)
        result.firstneigh[i] = result.firstneigh[i - 1] + result.numneigh[i - 1];

    result.total_pairs = (n > 0) ? result.firstneigh[n - 1] + result.numneigh[n - 1] : 0;
    result.neighbors.resize(static_cast<std::size_t>(result.total_pairs));

    auto write_ptr = result.firstneigh;
    for (auto& tp : thread_pairs) {
        for (auto& [i, j] : tp) {
            result.neighbors[write_ptr[i]] = j;
            write_ptr[i]++;
        }
    }

}

PairList CellList::build_half(const cpu::Frame& frame, int n_active) const {
    PairList out;
    refresh_half(frame, n_active, out);
    return out;
}

void CellList::refresh_full(const cpu::Frame& frame, int n_active,
                            PairList& result) const {
    const int n = (n_active < 0) ? static_cast<int>(frame.atoms.size())
                                 : n_active;
    const double rc_sq = cutoff_skin_ * cutoff_skin_;

    ensure_at_least(result.ilist, static_cast<std::size_t>(n));
    ensure_at_least(result.numneigh, static_cast<std::size_t>(n));
    ensure_at_least(result.firstneigh, static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        result.ilist[i] = i;
        result.numneigh[i] = 0;
    }
    result.inum = n;

    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
#endif

    std::vector<std::vector<std::pair<int, int>>> thread_pairs(n_threads);

    #pragma omp parallel
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        auto& my_pairs = thread_pairs[tid];

        std::vector<int> visited;
        visited.reserve(27);

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            const int ci = atom_cell_[i];
            const int cz0 = ci % ncz_;
            const int cy0 = (ci / ncz_) % ncy_;
            const int cx0 = ci / (ncy_ * ncz_);

            visited.clear();

            for (auto s : full_stencil_) {
                int cx = cx0 + s.di;
                int cy = cy0 + s.dj;
                int cz = cz0 + s.dk;

                if (axis_extended_[0]) {
                    if (cx < 0 || cx >= ncx_) continue;
                } else {
                    if (cx <  0)    cx += ncx_;
                    if (cx >= ncx_) cx -= ncx_;
                }
                if (axis_extended_[1]) {
                    if (cy < 0 || cy >= ncy_) continue;
                } else {
                    if (cy <  0)    cy += ncy_;
                    if (cy >= ncy_) cy -= ncy_;
                }
                if (axis_extended_[2]) {
                    if (cz < 0 || cz >= ncz_) continue;
                } else {
                    if (cz <  0)    cz += ncz_;
                    if (cz >= ncz_) cz -= ncz_;
                }

                const int cj = cx * ncy_ * ncz_ + cy * ncz_ + cz;

                bool seen = false;
                for (int v : visited) if (v == cj) { seen = true; break; }
                if (seen) continue;
                visited.push_back(cj);

                const int start = cell_offset_[cj];
                const int end   = start + cell_count_[cj];

                for (int idx = start; idx < end; ++idx) {
                    int j = binned_atoms_[idx];
                    if (j == i) continue;

                    double dx = frame.atoms.x[i] - frame.atoms.x[j];
                    double dy = frame.atoms.y[i] - frame.atoms.y[j];
                    double dz = frame.atoms.z[i] - frame.atoms.z[j];
                    frame.box.min_image(dx, dy, dz);

                    double r2 = dx * dx + dy * dy + dz * dz;
                    if (r2 < rc_sq) {
                        my_pairs.emplace_back(i, j);
                    }
                }
            }
        }
    }

    for (auto& tp : thread_pairs)
        for (auto& [i, j] : tp)
            result.numneigh[i]++;

    result.firstneigh[0] = 0;
    for (int i = 1; i < n; ++i)
        result.firstneigh[i] = result.firstneigh[i - 1] + result.numneigh[i - 1];

    result.total_pairs = (n > 0) ? result.firstneigh[n - 1] + result.numneigh[n - 1] : 0;
    result.neighbors.resize(static_cast<std::size_t>(result.total_pairs));

    auto write_ptr = result.firstneigh;
    for (auto& tp : thread_pairs) {
        for (auto& [i, j] : tp) {
            result.neighbors[write_ptr[i]] = j;
            write_ptr[i]++;
        }
    }
}

PairList CellList::build_full(const cpu::Frame& frame, int n_active) const {
    PairList out;
    refresh_full(frame, n_active, out);
    return out;
}

} // namespace md
} // namespace cpu
} // namespace ATV
