/// @file cpu/scf/shell_pair.cpp
/// @brief Build canonical shell pair list with Cauchy-Schwarz screening.

#include "atomiverse/cpu/scf/shell_pair.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/core/log.h"

#include <cmath>
#include <algorithm>

namespace ATV {
namespace cpu {
namespace scf {

ShellPairList build_shell_pair_list(const IntegralPack& pack, double thr) {
    const int n = pack.n_basis;

    ShellPairList spl;
    spl.n_basis   = n;
    spl.threshold = thr;

    // Copy backing storage and fix up shell pointers to own data.
    // pack.shells[i].alpha/coeff point into pack.shell_data; after copy they
    // must be redirected into spl.shell_data.
    spl.shell_data = pack.shell_data;
    spl.atoms      = pack.atoms;
    spl.shells.resize(n);
    for (int i = 0; i < n; ++i) {
        spl.shells[i]       = pack.shells[i];
        const int np        = spl.shells[i].n_prim;
        spl.shells[i].alpha = spl.shell_data[i].data();
        spl.shells[i].coeff = spl.shell_data[i].data() + np;
    }

    // Compute Schwarz bounds Q_ij = sqrt((ij|ij)) and collect screened pairs.
    spl.pairs.reserve(static_cast<std::size_t>(n * (n + 1) / 2));
    int n_screened = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            // (ij|ij) is non-negative by positive semi-definiteness of ERIs.
            const double diag = native::eri(spl.shells[i], spl.shells[j],
                                            spl.shells[i], spl.shells[j]);
            const double schwarz = std::sqrt(std::max(0.0, diag));
            if (schwarz > thr) {
                spl.pairs.push_back({i, j, schwarz});
            } else {
                ++n_screened;
            }
        }
    }

    ATV_DEBUG(
        "[shell_pair] built n_basis={} pairs={} screened={} thr={:.1e}",
        n, spl.pairs.size(), n_screened, thr);

    return spl;
}

ShellBlockPairList build_shell_block_pair_list(const IntegralPack& pack,
                                               double thr) {
    const int n = pack.n_basis;
    const int n_raw_shells = static_cast<int>(pack.shells.size());

    ShellBlockPairList sbl;
    sbl.n_basis   = n;
    sbl.threshold = thr;

    // ── Step 1: Collect raw block descriptors (two-pass to avoid realloc bugs)
    // Rules (matching integrals.cpp shell layout):
    //   l=0 shell → one s-block (n_cart=1).
    //   l=1 shell (px) followed by py, pz → one p-block (n_cart=3).
    //     All three share the same alpha/coeff; use the px shell's storage.
    struct RawBlock { int shell_data_idx, l, n_cart, ao_start; };
    std::vector<RawBlock> raw_blocks;
    raw_blocks.reserve(n_raw_shells); // upper bound

    {
        int shell_idx = 0, ao_idx = 0;
        while (shell_idx < n_raw_shells && ao_idx < n) {
            const native::Shell& sh = pack.shells[shell_idx];
            if (sh.l == 0) {
                raw_blocks.push_back({shell_idx, 0, 1, ao_idx});
                ++shell_idx;
                ao_idx += 1;
            } else {
                // l=1: consume px, py, pz (same exponents/coefficients)
                raw_blocks.push_back({shell_idx, 1, 3, ao_idx});
                shell_idx += 3;
                ao_idx    += 3;
            }
        }
    }

    const int n_blocks = static_cast<int>(raw_blocks.size());

    // ── Step 2: Allocate block_data all at once, then set up ShellBlock pointers.
    // block_data must be stable (no further push_back after this point) before
    // the ShellBlock alpha/coeff pointers into it are stored.
    sbl.block_data.resize(n_blocks);
    for (int i = 0; i < n_blocks; ++i)
        sbl.block_data[i] = pack.shell_data[raw_blocks[i].shell_data_idx];

    sbl.blocks.resize(n_blocks);
    for (int i = 0; i < n_blocks; ++i) {
        const RawBlock& rb = raw_blocks[i];
        const native::Shell& sh = pack.shells[rb.shell_data_idx];
        native::ShellBlock& blk = sbl.blocks[i];
        blk.n_prim    = sh.n_prim;
        blk.center[0] = sh.center[0];
        blk.center[1] = sh.center[1];
        blk.center[2] = sh.center[2];
        blk.l         = rb.l;
        blk.n_cart    = rb.n_cart;
        blk.ao_start  = rb.ao_start;
        blk.alpha     = sbl.block_data[i].data();
        blk.coeff     = sbl.block_data[i].data() + sh.n_prim;
    }

    ATV_DEBUG("[shell_pair] block grouping done n_basis={} n_blocks={}", n, n_blocks);

    // ── Step 2: Canonical block pairs + Schwarz bounds ────────────────
    // For each canonical pair (si >= sj), compute the block-level Schwarz
    // bound Q_{AB} = sqrt(max_{ia,ib} |(ia ib|ia ib)|) using eri_block().
    sbl.pairs.reserve(static_cast<std::size_t>(n_blocks * (n_blocks + 1) / 2));
    int n_screened = 0;

    // Temporary buffer for the self-ERI block: max (3*3*3*3)=81 doubles
    double blk_buf[81];

    for (int si = 0; si < n_blocks; ++si) {
        for (int sj = 0; sj <= si; ++sj) {
            const native::ShellBlock& Bsi = sbl.blocks[si];
            const native::ShellBlock& Bsj = sbl.blocks[sj];
            const int na = Bsi.n_cart, nb = Bsj.n_cart;

            // Compute eri_block(Bsi, Bsj, Bsi, Bsj)
            const int buf_size = na * nb * na * nb;
            for (int k = 0; k < buf_size; ++k) blk_buf[k] = 0.0;
            native::eri_block(Bsi, Bsj, Bsi, Bsj, blk_buf);

            // Block-level Schwarz bound: sqrt(max_{ia,ib} blk[ia,ib,ia,ib])
            double max_sq = 0.0;
            for (int ia = 0; ia < na; ++ia)
                for (int ib = 0; ib < nb; ++ib) {
                    // Index: [ia,ib,ia,ib] in layout (na, nb, na, nb)
                    const int idx = ((ia * nb + ib) * na + ia) * nb + ib;
                    if (blk_buf[idx] > max_sq) max_sq = blk_buf[idx];
                }

            const double schwarz = std::sqrt(std::max(0.0, max_sq));
            if (schwarz > thr) {
                sbl.pairs.push_back({si, sj, schwarz});
            } else {
                ++n_screened;
            }
        }
    }

    ATV_DEBUG(
        "[shell_pair] block pair list built n_blocks={} pairs={} screened={} thr={:.1e}",
        n_blocks, sbl.pairs.size(), n_screened, thr);

    return sbl;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
