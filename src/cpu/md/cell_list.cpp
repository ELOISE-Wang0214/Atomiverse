/// @file cpu/md/cell_list.cpp
/// @brief CellList implementation -- LAMMPS-style linked-cell neighbor list.

#include "atomiverse/cpu/md/cell_list.h"
#include "atomiverse/core/mat3.h"

#include <algorithm>
#include <cmath>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

// ── Minimum image helper ─────────────────────────────────────────

static inline void minimum_image(double& dx, double& dy, double& dz, const Mat3& box) {
    // Orthorhombic only: box(0,0), box(1,1), box(2,2) are edge lengths
    double lx = box(0, 0), ly = box(1, 1), lz = box(2, 2);
    if (dx >  lx * double(0.5)) dx -= lx;
    if (dx < -lx * double(0.5)) dx += lx;
    if (dy >  ly * double(0.5)) dy -= ly;
    if (dy < -ly * double(0.5)) dy += ly;
    if (dz >  lz * double(0.5)) dz -= lz;
    if (dz < -lz * double(0.5)) dz += lz;
}

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
    double lx_eff = frame.box(0, 0) + (axis_extended_[0] ? double(2) * cutoff_skin_ : double(0));
    double ly_eff = frame.box(1, 1) + (axis_extended_[1] ? double(2) * cutoff_skin_ : double(0));
    double lz_eff = frame.box(2, 2) + (axis_extended_[2] ? double(2) * cutoff_skin_ : double(0));

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
    double lx = frame.box(0, 0);
    double ly = frame.box(1, 1);
    double lz = frame.box(2, 2);

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
    box_at_last_build_ = frame.box(0, 0);
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

    double max_disp_sq = 0;

    for (int i = 0; i < n; ++i) {
        double dx = frame.atoms.x[i] - x_save_[i];
        double dy = frame.atoms.y[i] - y_save_[i];
        double dz = frame.atoms.z[i] - z_save_[i];
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

PairList CellList::build_half(const cpu::Frame& frame, int n_active) const {
    int n = (n_active < 0) ? frame.atoms.size() : n_active;
    double rc_sq = cutoff_skin_ * cutoff_skin_;

    PairList result;
    result.ilist = xt::arange<int>(n);
    result.inum = n;
    result.numneigh = xt::zeros<int>({n});
    result.firstneigh = xt::zeros<int>({n});

    // Per-thread pair buffers
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

        // Per-thread deduplication buffer for visited neighbor cells.
        // When a small grid (ncx <= 2 etc.) has the FULL stencil wrap
        // multiple distinct offsets onto the same cj, we'd visit that
        // cell more than once; dedup ensures each (i, cj) lookup is
        // unique so pair counts stay correct.
        std::vector<int> visited;
        visited.reserve(27);

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            const int ci = atom_cell_[i];
            const int cz0 = ci % ncz_;
            const int cy0 = (ci / ncz_) % ncy_;
            const int cx0 = ci / (ncy_ * ncz_);

            visited.clear();

            // Iterate the FULL stencil (27 cells incl. self) and filter
            // each candidate by `j > i`. Using the half-stencil instead
            // would skip cells reachable only via PBC wrap (e.g. cell
            // ncx-1 from cell 0 -- offset di=-1 wraps to the same cell
            // di=+(ncx-1) but the half-stencil only stores small flat
            // offsets), causing missed pairs at every box edge.
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

                // Skip if we already processed this neighbor cell (only
                // happens when ncN <= 2 and PBC wraps stencil offsets
                // onto the same cell; for typical large grids this loop
                // does at most one O(visited.size()) check per stencil
                // entry, which dwarfs the O(atoms-per-cell) work below).
                bool seen = false;
                for (int v : visited) if (v == cj) { seen = true; break; }
                if (seen) continue;
                visited.push_back(cj);

                const int start = cell_offset_[cj];
                const int end   = start + cell_count_[cj];

                for (int idx = start; idx < end; ++idx) {
                    int j = binned_atoms_[idx];
                    if (j <= i) continue;  // Half list: i < j

                    double dx = frame.atoms.x[i] - frame.atoms.x[j];
                    double dy = frame.atoms.y[i] - frame.atoms.y[j];
                    double dz = frame.atoms.z[i] - frame.atoms.z[j];
                    minimum_image(dx, dy, dz, frame.box);

                    double r2 = dx * dx + dy * dy + dz * dz;
                    if (r2 < rc_sq) {
                        my_pairs.emplace_back(i, j);
                    }
                }
            }
        }
    }

    // Merge per-thread pairs into CSR
    // First count per atom
    for (auto& tp : thread_pairs)
        for (auto& [i, j] : tp)
            result.numneigh[i]++;

    // Prefix sum for firstneigh
    result.firstneigh[0] = 0;
    for (int i = 1; i < n; ++i)
        result.firstneigh[i] = result.firstneigh[i - 1] + result.numneigh[i - 1];

    result.total_pairs = (n > 0) ? result.firstneigh[n - 1] + result.numneigh[n - 1] : 0;
    result.neighbors = xt::zeros<int>({result.total_pairs});

    // Write pointers (copy of firstneigh)
    auto write_ptr = result.firstneigh;
    for (auto& tp : thread_pairs) {
        for (auto& [i, j] : tp) {
            result.neighbors[write_ptr[i]] = j;
            write_ptr[i]++;
        }
    }

    return result;
}

PairList CellList::build_full(const cpu::Frame& frame, int n_active) const {
    int n = (n_active < 0) ? frame.atoms.size() : n_active;
    double rc_sq = cutoff_skin_ * cutoff_skin_;

    PairList result;
    result.ilist = xt::arange<int>(n);
    result.inum = n;
    result.numneigh = xt::zeros<int>({n});
    result.firstneigh = xt::zeros<int>({n});

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
                    minimum_image(dx, dy, dz, frame.box);

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
    result.neighbors = xt::zeros<int>({result.total_pairs});

    auto write_ptr = result.firstneigh;
    for (auto& tp : thread_pairs) {
        for (auto& [i, j] : tp) {
            result.neighbors[write_ptr[i]] = j;
            write_ptr[i]++;
        }
    }

    return result;
}

} // namespace md
} // namespace cpu
} // namespace ATV
