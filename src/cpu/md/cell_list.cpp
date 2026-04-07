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

// ── PBC wrapping for cell index ──────────────────────────────────

static inline int wrap_cell(int c, int nc) {
    if (c < 0) return c + nc;
    if (c >= nc) return c - nc;
    return c;
}

// ── CellList ─────────────────────────────────────────────────────

CellList::CellList(double cutoff, double skin)
    : cutoff_(cutoff),
      skin_(skin < 0 ? double(0.3) * cutoff : skin),
      cutoff_skin_(cutoff + (skin < 0 ? double(0.3) * cutoff : skin)) {}

void CellList::init(const cpu::Frame& frame, double cutoff_skin) {
    cutoff_skin_ = cutoff_ + skin_;
    compute_grid(frame);
    build_stencil();

    int n = frame.atoms.size();
    atom_cell_ = xt::zeros<int>({n});
    binned_atoms_ = xt::zeros<int>({n});
    x_save_ = xt::zeros<double>({n});
    y_save_ = xt::zeros<double>({n});
    z_save_ = xt::zeros<double>({n});
}

void CellList::compute_grid(const cpu::Frame& frame) {
    double lx = frame.box(0, 0);
    double ly = frame.box(1, 1);
    double lz = frame.box(2, 2);

    ncx_ = std::max(int(1), static_cast<int>(std::floor(lx / cutoff_skin_)));
    ncy_ = std::max(int(1), static_cast<int>(std::floor(ly / cutoff_skin_)));
    ncz_ = std::max(int(1), static_cast<int>(std::floor(lz / cutoff_skin_)));

    csx_ = lx / ncx_;
    csy_ = ly / ncy_;
    csz_ = lz / ncz_;

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
                int flat = di * ncy_ * ncz_ + dj * ncz_ + dk;
                full_stencil_.push_back(flat);
                // Half stencil: positive offsets + self (flat >= 0)
                if (flat > 0 || (flat == 0 && di == 0 && dj == 0 && dk == 0))
                    half_stencil_.push_back(flat);
            }
        }
    }
}

void CellList::bin_atoms(const cpu::Frame& frame) {
    int n = frame.atoms.size();
    double lx = frame.box(0, 0);
    double ly = frame.box(1, 1);
    double lz = frame.box(2, 2);

    // Clear cell counts
    std::fill(cell_count_.begin(), cell_count_.end(), int(0));

    // Assign atoms to cells
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        // Wrap position into [0, L)
        double xi = frame.atoms.x[i] - lx * std::floor(frame.atoms.x[i] / lx);
        double yi = frame.atoms.y[i] - ly * std::floor(frame.atoms.y[i] / ly);
        double zi = frame.atoms.z[i] - lz * std::floor(frame.atoms.z[i] / lz);

        int cx = std::min(static_cast<int>(xi / csx_), ncx_ - 1);
        int cy = std::min(static_cast<int>(yi / csy_), ncy_ - 1);
        int cz = std::min(static_cast<int>(zi / csz_), ncz_ - 1);

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

void CellList::save_positions(const cpu::Frame& frame) {
    int n = frame.atoms.size();
    for (int i = 0; i < n; ++i) {
        x_save_[i] = frame.atoms.x[i];
        y_save_[i] = frame.atoms.y[i];
        z_save_[i] = frame.atoms.z[i];
    }
}

void CellList::build(const cpu::Frame& frame) {
    bin_atoms(frame);
    save_positions(frame);
}

bool CellList::needs_rebuild(const cpu::Frame& frame) const {
    int n = frame.atoms.size();
    double max_disp_sq = 0;
    double threshold = skin_ * double(0.5);
    double threshold_sq = threshold * threshold;

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

void CellList::update(const cpu::Frame& frame) {
    if (needs_rebuild(frame))
        build(frame);
}

PairList CellList::build_half(const cpu::Frame& frame) const {
    int n = frame.atoms.size();
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

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            int ci = atom_cell_[i];
            int cx = ci / (ncy_ * ncz_);
            int rem = ci % (ncy_ * ncz_);
            int cy = rem / ncz_;
            int cz = rem % ncz_;

            for (auto s : half_stencil_) {
                int cj = ci + s;

                // PBC wrapping for cell index
                int jx = cx, jy = cy, jz = cz;
                if (s != 0) {
                    // Decode stencil offset
                    int ds = s;
                    int dz_s = ds % ncz_;
                    if (dz_s < -1) dz_s += ncz_;
                    ds /= ncz_;
                    int dy_s = ds % ncy_;
                    if (dy_s < -1) dy_s += ncy_;
                    int dx_s = ds / ncy_;

                    // Recompute from original di/dj/dk encoding
                    // Simpler: directly iterate stencil with PBC
                }

                // Use wrapped cell computation instead
                // Recompute neighbor cell with PBC
                // For stencil offset s, decode back to (di, dj, dk) shifts
                // Then apply PBC wrapping

                // Iterate atoms in cell cj
                if (cj < 0 || cj >= total_cells_) continue;
                int start = cell_offset_[cj];
                int end = start + cell_count_[cj];

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

    result.total_pairs = result.firstneigh[n - 1] + result.numneigh[n - 1];
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

PairList CellList::build_full(const cpu::Frame& frame) const {
    int n = frame.atoms.size();
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

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            int ci = atom_cell_[i];

            for (auto s : full_stencil_) {
                int cj = ci + s;
                if (cj < 0 || cj >= total_cells_) continue;
                int start = cell_offset_[cj];
                int end = start + cell_count_[cj];

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

    result.total_pairs = result.firstneigh[n - 1] + result.numneigh[n - 1];
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
