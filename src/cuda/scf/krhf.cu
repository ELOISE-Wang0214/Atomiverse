/// @file cuda/scf/krhf.cu
/// @brief Device SCF step for cuda::scf::KRHF.
///
/// Setup (mesh, AO, Hcore, canonical X, Coulomb kernels) is
/// `cpu::scf::internal::krhf_assemble` — geometry and basis, not the
/// density. The iteration builds J (`fft_jk.get_j_kpts`), K
/// (`get_k_kpts`), the Fock, `zheevd`, and the density on device.

#include "atomiverse/cuda/scf/krhf.cuh"

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/scf/internal/kpoint_hf.h"

#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

namespace ATV {
namespace cuda {
namespace scf {

namespace {

struct Z {
    double r;
    double i;
};

void check_cuda(cudaError_t e, const char* what) {
    if (e != cudaSuccess)
        throw std::runtime_error(std::string("cuda::scf::KRHF: ") + what + ": "
                                 + cudaGetErrorString(e));
}

void check_sol(cusolverStatus_t s, const char* what) {
    if (s != CUSOLVER_STATUS_SUCCESS)
        throw std::runtime_error(std::string("cuda::scf::KRHF: ") + what
                                 + " cusolver status="
                                 + std::to_string(int(s)));
}

__global__ void k_dft_line(const Z* in, Z* out, int n, int nlines, int kind,
                           int nx, int ny, int nz, int inverse) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nlines * n)
        return;
    const int line = idx / n;
    const int k = idx - line * n;
    int base = 0;
    int stride = 1;
    if (kind == 0) {
        base = line * nz;
        stride = 1;
    } else if (kind == 1) {
        const int ix = line / nz;
        const int iz = line - ix * nz;
        base = iz + nz * ny * ix;
        stride = nz;
    } else {
        const int iy = line / nz;
        const int iz = line - iy * nz;
        base = iz + nz * iy;
        stride = ny * nz;
    }
    const double sign = inverse ? 1.0 : -1.0;
    const double inv = inverse ? 1.0 / static_cast<double>(n) : 1.0;
    double sr = 0.0, si = 0.0;
    for (int j = 0; j < n; ++j) {
        const Z v = in[base + j * stride];
        const double ang = sign * 2.0 * 3.14159265358979323846
                           * static_cast<double>(k * j)
                           / static_cast<double>(n);
        const double c = cos(ang);
        const double s = sin(ang);
        sr += v.r * c - v.i * s;
        si += v.r * s + v.i * c;
    }
    out[base + k * stride] = Z{sr * inv, si * inv};
}

void fft3_dev(Z* data, Z* tmp, int nx, int ny, int nz, bool inverse,
              cudaStream_t stream) {
    auto launch = [&](int kind, int n, int nlines) {
        const int ntot = nlines * n;
        const int blk = 128;
        const int grid = (ntot + blk - 1) / blk;
        k_dft_line<<<grid, blk, 0, stream>>>(data, tmp, n, nlines, kind, nx, ny,
                                             nz, inverse ? 1 : 0);
        check_cuda(cudaMemcpyAsync(data, tmp, sizeof(Z) * nx * ny * nz,
                                   cudaMemcpyDeviceToDevice, stream),
                   "fft copy");
    };
    launch(0, nz, nx * ny);
    launch(1, ny, nx * nz);
    launch(2, nx, ny * nz);
}

__global__ void k_rho(const Z* ao, const Z* dm, Z* rho, int nk, int nao,
                      int ng) {
    const int g = blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= ng)
        return;
    double rr = 0.0;
    for (int k = 0; k < nk; ++k) {
        const Z* a = ao + (static_cast<long long>(k) * ng + g) * nao;
        const Z* d = dm + static_cast<long long>(k) * nao * nao;
        for (int nu = 0; nu < nao; ++nu) {
            double ar = 0.0, ai = 0.0;
            for (int mu = 0; mu < nao; ++mu) {
                const Z am = a[mu];
                const Z dj = d[mu * nao + nu];
                ar += am.r * dj.r - am.i * dj.i;
                ai += am.r * dj.i + am.i * dj.r;
            }
            const Z an = a[nu];
            rr += ar * an.r + ai * an.i;
        }
    }
    rho[g] = Z{rr / static_cast<double>(nk), 0.0};
}

__global__ void k_scale(Z* v, const double* coul, int ng) {
    const int g = blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= ng)
        return;
    v[g].r *= coul[g];
    v[g].i *= coul[g];
}

__global__ void k_j_contract(const Z* ao, const Z* vr, Z* vj, int nk, int nao,
                             int ng, double weight) {
    const int job = blockIdx.x * blockDim.x + threadIdx.x;
    const int njob = nk * nao * nao;
    if (job >= njob)
        return;
    const int nu = job % nao;
    const int mu = (job / nao) % nao;
    const int k = job / (nao * nao);
    double sr = 0.0, si = 0.0;
    for (int g = 0; g < ng; ++g) {
        const Z a_mu = ao[(static_cast<long long>(k) * ng + g) * nao + mu];
        const Z a_nu = ao[(static_cast<long long>(k) * ng + g) * nao + nu];
        const double w = vr[g].r * weight;
        // conj(a_mu) * w * a_nu
        sr += (a_mu.r * a_nu.r + a_mu.i * a_nu.i) * w;
        si += (a_mu.r * a_nu.i - a_mu.i * a_nu.r) * w;
    }
    vj[static_cast<long long>(k) * nao * nao + mu * nao + nu] = Z{sr, si};
}

__global__ void k_hermitize(Z* a, int nk, int nao) {
    const int k = blockIdx.x;
    if (k >= nk)
        return;
    Z* m = a + static_cast<long long>(k) * nao * nao;
    for (int i = threadIdx.x; i < nao; i += blockDim.x) {
        m[i * nao + i].i = 0.0;
        for (int j = i + 1; j < nao; ++j) {
            const Z ij = m[i * nao + j];
            const Z ji = m[j * nao + i];
            const double hr = 0.5 * (ij.r + ji.r);
            const double hi = 0.5 * (ij.i - ji.i);
            m[i * nao + j] = Z{hr, hi};
            m[j * nao + i] = Z{hr, -hi};
        }
    }
}

__global__ void k_vhf(const Z* vj, const Z* vk, Z* v, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;
    v[i] = Z{vj[i].r - 0.5 * vk[i].r, vj[i].i - 0.5 * vk[i].i};
}

__global__ void k_fock(const Z* H, const Z* vj, const Z* vk, Z* F, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
        return;
    F[i] =
        Z{H[i].r + vj[i].r - 0.5 * vk[i].r, H[i].i + vj[i].i - 0.5 * vk[i].i};
}

__global__ void k_trace_energy(const Z* dm, const Z* H, const Z* vhf, int nk,
                               int nao, double* e_out) {
    // One thread: nk and nao are tiny on this path.
    double e1 = 0.0, e2 = 0.0;
    for (int k = 0; k < nk; ++k) {
        const Z* d = dm + static_cast<long long>(k) * nao * nao;
        const Z* h = H + static_cast<long long>(k) * nao * nao;
        const Z* v = vhf + static_cast<long long>(k) * nao * nao;
        for (int i = 0; i < nao; ++i)
            for (int j = 0; j < nao; ++j) {
                const Z dij = d[i * nao + j];
                const Z hij = h[j * nao + i];
                const Z vij = v[j * nao + i];
                e1 += dij.r * hij.r - dij.i * hij.i;
                e2 += dij.r * vij.r - dij.i * vij.i;
            }
    }
    e_out[0] = (e1 + 0.5 * e2) / static_cast<double>(nk);
}

__global__ void k_aodm(const Z* ao, const Z* dm, Z* aodm, int nao, int ng,
                       int k) {
    const int g = blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= ng)
        return;
    const Z* a = ao + (static_cast<long long>(k) * ng + g) * nao;
    const Z* d = dm + static_cast<long long>(k) * nao * nao;
    for (int mu = 0; mu < nao; ++mu) {
        double sr = 0.0, si = 0.0;
        for (int nu = 0; nu < nao; ++nu) {
            const Z dj = d[mu * nao + nu];
            const Z an = a[nu];
            // dm * conj(ao)
            sr += dj.r * an.r + dj.i * an.i;
            si += dj.i * an.r - dj.r * an.i;
        }
        aodm[static_cast<long long>(mu) * ng + g] = Z{sr, si};
    }
}

__global__ void k_rho1(const Z* ao1, const Z* ao2, const double* coords,
                       double kdx, double kdy, double kdz, int zero, Z* rho,
                       int nao, int ng, int k1, int k2) {
    const int g = blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= ng)
        return;
    const Z* a1 = ao1 + (static_cast<long long>(k1) * ng + g) * nao;
    const Z* a2 = ao2 + (static_cast<long long>(k2) * ng + g) * nao;
    double pr = 1.0, pi = 0.0;
    if (!zero) {
        const double dot = coords[g * 3] * kdx + coords[g * 3 + 1] * kdy
                           + coords[g * 3 + 2] * kdz;
        pr = cos(-dot);
        pi = sin(-dot);
    }
    for (int i = 0; i < nao; ++i)
        for (int j = 0; j < nao; ++j) {
            // conj(ao1[i]) * phase * ao2[j]
            const Z u = a1[i];
            const Z v = a2[j];
            const double cr = u.r * pr + u.i * pi;
            const double ci = u.r * pi - u.i * pr;
            const double sr = cr * v.r - ci * v.i;
            const double si = cr * v.i + ci * v.r;
            rho[(static_cast<long long>(i * nao + j) * ng) + g] = Z{sr, si};
        }
}

__global__ void k_vrdm(const Z* rho, const Z* aodm, const double* coords,
                       double kdx, double kdy, double kdz, int zero, Z* vrdm,
                       int nao, int ng) {
    const int g = blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= ng)
        return;
    double pr = 1.0, pi = 0.0;
    if (!zero) {
        const double dot = coords[g * 3] * kdx + coords[g * 3 + 1] * kdy
                           + coords[g * 3 + 2] * kdz;
        pr = cos(dot);
        pi = sin(dot);
    }
    for (int i = 0; i < nao; ++i) {
        double sr = 0.0, si = 0.0;
        for (int j = 0; j < nao; ++j) {
            const Z r = rho[(static_cast<long long>(i * nao + j) * ng) + g];
            const Z a = aodm[static_cast<long long>(j) * ng + g];
            sr += r.r * a.r - r.i * a.i;
            si += r.r * a.i + r.i * a.r;
        }
        const double tr = sr * pr - si * pi;
        const double ti = sr * pi + si * pr;
        vrdm[static_cast<long long>(i) * ng + g] = Z{tr, ti};
    }
}

__global__ void k_k_accum(const Z* vrdm, const Z* ao, Z* vk, int nao, int ng,
                          int k1, double weight) {
    const int job = blockIdx.x * blockDim.x + threadIdx.x;
    if (job >= nao * nao)
        return;
    const int p = job % nao;
    const int i = job / nao;
    double sr = 0.0, si = 0.0;
    for (int g = 0; g < ng; ++g) {
        const Z v = vrdm[static_cast<long long>(i) * ng + g];
        const Z a = ao[(static_cast<long long>(k1) * ng + g) * nao + p];
        sr += v.r * a.r - v.i * a.i;
        si += v.r * a.i + v.i * a.r;
    }
    Z* dst = vk + static_cast<long long>(k1) * nao * nao + i * nao + p;
    dst->r += weight * sr;
    dst->i += weight * si;
}

__global__ void k_zero(Z* a, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
        a[i] = Z{0.0, 0.0};
}

__global__ void k_xform(const Z* X, const Z* F, Z* Hp, int nao, int nmo,
                        int k) {
    // Hp = X^H F X, one block.
    const Z* x = X + static_cast<long long>(k) * nao * nao;
    const Z* f = F + static_cast<long long>(k) * nao * nao;
    for (int a = threadIdx.x; a < nmo; a += blockDim.x) {
        for (int b = 0; b < nmo; ++b) {
            double sr = 0.0, si = 0.0;
            for (int p = 0; p < nao; ++p) {
                double fr = 0.0, fi = 0.0;
                for (int q = 0; q < nao; ++q) {
                    const Z fq = f[p * nao + q];
                    const Z xq = x[q * nao + b];
                    fr += fq.r * xq.r - fq.i * xq.i;
                    fi += fq.r * xq.i + fq.i * xq.r;
                }
                const Z xp = x[p * nao + a];
                // conj(X[p,a]) * (F X)[p,b]
                sr += xp.r * fr + xp.i * fi;
                si += xp.r * fi - xp.i * fr;
            }
            Hp[a * nmo + b] = Z{sr, si};
        }
    }
}

__global__ void k_row_to_col(const Z* row, cuDoubleComplex* col, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n * n)
        return;
    const int r = i / n;
    const int c = i - r * n;
    col[c * n + r] = make_cuDoubleComplex(row[r * n + c].r, row[r * n + c].i);
}

__global__ void k_apply_x(const Z* X, const cuDoubleComplex* Cp, Z* C, int nao,
                          int nmo, int k) {
    // C_ao[:, j] = X * Cp[:, j], Cp column-major.
    const int j = blockIdx.x;
    if (j >= nmo)
        return;
    const Z* x = X + static_cast<long long>(k) * nao * nao;
    Z* c = C + static_cast<long long>(k) * nao * nao;
    for (int i = threadIdx.x; i < nao; i += blockDim.x) {
        double sr = 0.0, si = 0.0;
        for (int a = 0; a < nmo; ++a) {
            const Z xa = x[i * nao + a];
            const cuDoubleComplex cp = Cp[a + j * nmo];
            sr += xa.r * cuCreal(cp) - xa.i * cuCimag(cp);
            si += xa.r * cuCimag(cp) + xa.i * cuCreal(cp);
        }
        c[i * nao + j] = Z{sr, si};
    }
}

__global__ void k_dm_from_c(const Z* C, const double* occ, Z* dm, int nao,
                            int k) {
    const Z* c = C + static_cast<long long>(k) * nao * nao;
    Z* d = dm + static_cast<long long>(k) * nao * nao;
    for (int i = threadIdx.x; i < nao; i += blockDim.x) {
        for (int j = 0; j < nao; ++j) {
            double sr = 0.0, si = 0.0;
            for (int a = 0; a < nao; ++a) {
                if (!(occ[a] > 0.0))
                    continue;
                const Z ci = c[i * nao + a];
                const Z cj = c[j * nao + a];
                sr += occ[a] * (ci.r * cj.r + ci.i * cj.i);
                si += occ[a] * (ci.i * cj.r - ci.r * cj.i);
            }
            d[i * nao + j] = Z{sr, si};
        }
    }
}

int div_up(int n, int b) { return (n + b - 1) / b; }

} // namespace

KRHF::KRHF(cpu::scf::BasisSet basis_in, std::array<int, 3> k_mesh_in,
           int charge_in, double tol_in, int max_iter_in)
    : basis(basis_in), k_mesh(k_mesh_in), charge(charge_in), tol(tol_in),
      max_iter(max_iter_in) {}

void KRHF::solve(const cpu::Frame& frame, bool write, cpu::Frame* out) {
    const int mesh_k[3] = {k_mesh[0], k_mesh[1], k_mesh[2]};
    const cpu::scf::internal::KrhfStatic sys =
        cpu::scf::internal::krhf_assemble(frame, basis, mesh_k, charge);
    const int nk = sys.nk;
    const int nao = sys.nao;
    const int ng = sys.ng;
    const int nx = sys.mesh[0], ny = sys.mesh[1], nz = sys.mesh[2];
    if (nk < 1 || nao < 1 || ng != nx * ny * nz)
        throw std::runtime_error("cuda::scf::KRHF: empty k-point system");

    auto hz = [](const std::complex<double>& z) {
        return Z{z.real(), z.imag()};
    };
    std::vector<Z> h_ao(sys.ao.size()), h_H(sys.H.size()), h_X(sys.X.size());
    for (std::size_t i = 0; i < sys.ao.size(); ++i)
        h_ao[i] = hz(sys.ao[i]);
    for (std::size_t i = 0; i < sys.H.size(); ++i)
        h_H[i] = hz(sys.H[i]);
    for (std::size_t i = 0; i < sys.X.size(); ++i)
        h_X[i] = hz(sys.X[i]);

    const std::size_t nmat = static_cast<std::size_t>(nk) * nao * nao;
    Z *d_ao = nullptr, *d_H = nullptr, *d_X = nullptr, *d_dm = nullptr;
    Z *d_vj = nullptr, *d_vk = nullptr, *d_F = nullptr, *d_vhf = nullptr;
    Z *d_rho = nullptr, *d_tmp = nullptr, *d_rho1 = nullptr, *d_aodm = nullptr;
    Z *d_vrdm = nullptr, *d_Hp = nullptr;
    double *d_coul_j = nullptr, *d_coul_k = nullptr, *d_coords = nullptr;
    double* d_e = nullptr;
    cuDoubleComplex* d_col = nullptr;
    cuDoubleComplex* d_work = nullptr;
    double* d_eig = nullptr;
    int* d_info = nullptr;
    double* d_occ = nullptr;
    cudaStream_t stream = nullptr;
    check_cuda(cudaStreamCreate(&stream), "stream");
    auto guard = [&]() {
        cudaFree(d_ao);
        cudaFree(d_H);
        cudaFree(d_X);
        cudaFree(d_dm);
        cudaFree(d_vj);
        cudaFree(d_vk);
        cudaFree(d_F);
        cudaFree(d_vhf);
        cudaFree(d_rho);
        cudaFree(d_tmp);
        cudaFree(d_rho1);
        cudaFree(d_aodm);
        cudaFree(d_vrdm);
        cudaFree(d_Hp);
        cudaFree(d_coul_j);
        cudaFree(d_coul_k);
        cudaFree(d_coords);
        cudaFree(d_e);
        cudaFree(d_col);
        cudaFree(d_work);
        cudaFree(d_eig);
        cudaFree(d_info);
        cudaFree(d_occ);
        d_ao = d_H = d_X = d_dm = d_vj = d_vk = d_F = d_vhf = nullptr;
        d_rho = d_tmp = d_rho1 = d_aodm = d_vrdm = d_Hp = nullptr;
        d_coul_j = d_coul_k = d_coords = d_e = nullptr;
        d_col = d_work = nullptr;
        d_eig = nullptr;
        d_info = nullptr;
        d_occ = nullptr;
        if (stream) {
            cudaStreamDestroy(stream);
            stream = nullptr;
        }
    };
    try {
        check_cuda(cudaMalloc(&d_ao, h_ao.size() * sizeof(Z)), "ao");
        check_cuda(cudaMalloc(&d_H, nmat * sizeof(Z)), "H");
        check_cuda(cudaMalloc(&d_X, nmat * sizeof(Z)), "X");
        check_cuda(cudaMalloc(&d_dm, nmat * sizeof(Z)), "dm");
        check_cuda(cudaMalloc(&d_vj, nmat * sizeof(Z)), "vj");
        check_cuda(cudaMalloc(&d_vk, nmat * sizeof(Z)), "vk");
        check_cuda(cudaMalloc(&d_F, nmat * sizeof(Z)), "F");
        check_cuda(cudaMalloc(&d_vhf, nmat * sizeof(Z)), "vhf");
        check_cuda(cudaMalloc(&d_rho, static_cast<std::size_t>(ng) * sizeof(Z)),
                   "rho");
        check_cuda(cudaMalloc(&d_tmp, static_cast<std::size_t>(ng) * sizeof(Z)),
                   "tmp");
        check_cuda(cudaMalloc(&d_rho1, static_cast<std::size_t>(nao * nao * ng)
                                           * sizeof(Z)),
                   "rho1");
        check_cuda(
            cudaMalloc(&d_aodm, static_cast<std::size_t>(nao * ng) * sizeof(Z)),
            "aodm");
        check_cuda(
            cudaMalloc(&d_vrdm, static_cast<std::size_t>(nao * ng) * sizeof(Z)),
            "vrdm");
        check_cuda(
            cudaMalloc(&d_Hp, static_cast<std::size_t>(nao * nao) * sizeof(Z)),
            "Hp");
        check_cuda(cudaMalloc(&d_coul_j,
                              static_cast<std::size_t>(ng) * sizeof(double)),
                   "coul_j");
        check_cuda(cudaMalloc(&d_coul_k, sys.coul_k.size() * sizeof(double)),
                   "coul_k");
        check_cuda(cudaMalloc(&d_coords, sys.coords.size() * sizeof(double)),
                   "coords");
        check_cuda(cudaMalloc(&d_e, sizeof(double)), "e");
        check_cuda(cudaMalloc(&d_col, static_cast<std::size_t>(nao * nao)
                                          * sizeof(cuDoubleComplex)),
                   "col");
        check_cuda(
            cudaMalloc(&d_eig, static_cast<std::size_t>(nao) * sizeof(double)),
            "eig");
        check_cuda(cudaMalloc(&d_info, sizeof(int)), "info");
        check_cuda(cudaMemcpyAsync(d_ao, h_ao.data(), h_ao.size() * sizeof(Z),
                                   cudaMemcpyHostToDevice, stream),
                   "H2D ao");
        check_cuda(cudaMemcpyAsync(d_H, h_H.data(), nmat * sizeof(Z),
                                   cudaMemcpyHostToDevice, stream),
                   "H2D H");
        check_cuda(cudaMemcpyAsync(d_X, h_X.data(), nmat * sizeof(Z),
                                   cudaMemcpyHostToDevice, stream),
                   "H2D X");
        check_cuda(
            cudaMemcpyAsync(d_coul_j, sys.coul_j.data(),
                            static_cast<std::size_t>(ng) * sizeof(double),
                            cudaMemcpyHostToDevice, stream),
            "H2D coul");
        check_cuda(cudaMemcpyAsync(d_coul_k, sys.coul_k.data(),
                                   sys.coul_k.size() * sizeof(double),
                                   cudaMemcpyHostToDevice, stream),
                   "H2D coulk");
        check_cuda(cudaMemcpyAsync(d_coords, sys.coords.data(),
                                   sys.coords.size() * sizeof(double),
                                   cudaMemcpyHostToDevice, stream),
                   "H2D coords");

        cusolverDnHandle_t sol;
        check_sol(cusolverDnCreate(&sol), "create");
        check_sol(cusolverDnSetStream(sol, stream), "set stream");
        int lwork = 0;
        check_sol(cusolverDnZheevd_bufferSize(sol, CUSOLVER_EIG_MODE_VECTOR,
                                              CUBLAS_FILL_MODE_LOWER, nao,
                                              d_col, nao, d_eig, &lwork),
                  "buffer");
        check_cuda(cudaMalloc(&d_work, static_cast<std::size_t>(lwork)
                                           * sizeof(cuDoubleComplex)),
                   "work");

        auto diagonalize_k = [&](const Z* Fk, int k, std::vector<double>& eps,
                                 Z* Ck_unused) {
            (void)Ck_unused;
            const int m = sys.nmo[static_cast<std::size_t>(k)];
            k_xform<<<1, 64, 0, stream>>>(d_X, Fk, d_Hp, nao, m, k);
            k_row_to_col<<<(m * m + 127) / 128, 128, 0, stream>>>(d_Hp, d_col,
                                                                  m);
            check_sol(cusolverDnZheevd(sol, CUSOLVER_EIG_MODE_VECTOR,
                                       CUBLAS_FILL_MODE_LOWER, m, d_col, m,
                                       d_eig, d_work, lwork, d_info),
                      "zheevd");
            int info = 0;
            check_cuda(cudaMemcpyAsync(&info, d_info, sizeof(int),
                                       cudaMemcpyDeviceToHost, stream),
                       "D2H info");
            check_cuda(cudaStreamSynchronize(stream), "sync zheevd");
            if (info != 0)
                throw std::runtime_error("cuda::scf::KRHF: zheevd info="
                                         + std::to_string(info));
            eps.assign(static_cast<std::size_t>(nao), 1.0e8);
            check_cuda(
                cudaMemcpyAsync(eps.data(), d_eig,
                                sizeof(double) * static_cast<std::size_t>(m),
                                cudaMemcpyDeviceToHost, stream),
                "D2H eig");
            check_cuda(cudaStreamSynchronize(stream), "sync eig");
            k_apply_x<<<m, 32, 0, stream>>>(d_X, d_col, d_F, nao, m, k);
            // Eigenvectors land in d_F temporarily; caller copies to C.
        };

        // Core-Hamiltonian guess, diagonalized on device.
        std::vector<Z> hC(nmat);
        std::vector<std::vector<double>> eps(static_cast<std::size_t>(nk));
        // Use d_vj as orbital-coefficient scratch (overwritten by J later).
        for (int k = 0; k < nk; ++k) {
            diagonalize_k(d_H, k, eps[static_cast<std::size_t>(k)], nullptr);
            check_cuda(cudaMemcpyAsync(
                           hC.data() + static_cast<std::size_t>(k) * nao * nao,
                           d_F + static_cast<std::size_t>(k) * nao * nao,
                           sizeof(Z) * nao * nao, cudaMemcpyDeviceToHost,
                           stream),
                       "D2H C");
        }
        check_cuda(cudaStreamSynchronize(stream), "sync C");
        // Occupy on the host from the eigenvalues only (khf.get_occ).
        struct Item {
            double e;
            int k, a;
        };
        std::vector<Item> all;
        for (int k = 0; k < nk; ++k)
            for (int a = 0; a < nao; ++a)
                all.push_back({eps[static_cast<std::size_t>(k)]
                                  [static_cast<std::size_t>(a)],
                               k, a});
        std::stable_sort(
            all.begin(), all.end(),
            [](const Item& x, const Item& y) { return x.e < y.e; });
        std::vector<double> hocc(static_cast<std::size_t>(nk * nao), 0.0);
        for (int i = 0; i < sys.nocc; ++i) {
            const Item& it = all[static_cast<std::size_t>(i)];
            hocc[static_cast<std::size_t>(it.k * nao + it.a)] = 2.0;
        }
        // Put C back and build dm on device.
        check_cuda(cudaMemcpyAsync(d_F, hC.data(), nmat * sizeof(Z),
                                   cudaMemcpyHostToDevice, stream),
                   "H2D C");
        check_cuda(cudaMalloc(&d_occ, hocc.size() * sizeof(double)), "occ");
        check_cuda(cudaMemcpyAsync(d_occ, hocc.data(),
                                   hocc.size() * sizeof(double),
                                   cudaMemcpyHostToDevice, stream),
                   "H2D occ");
        for (int k = 0; k < nk; ++k)
            k_dm_from_c<<<1, 32, 0, stream>>>(
                d_F, d_occ + static_cast<std::size_t>(k) * nao, d_dm, nao, k);
        k_hermitize<<<nk, 32, 0, stream>>>(d_dm, nk, nao);

        const int blk = 128;
        const double jw = sys.vol / static_cast<double>(ng);
        const double kw = (1.0 / static_cast<double>(nk))
                          * (sys.vol / static_cast<double>(ng));

        auto build_jk = [&]() {
            k_rho<<<div_up(ng, blk), blk, 0, stream>>>(d_ao, d_dm, d_rho, nk,
                                                       nao, ng);
            fft3_dev(d_rho, d_tmp, nx, ny, nz, false, stream);
            k_scale<<<div_up(ng, blk), blk, 0, stream>>>(d_rho, d_coul_j, ng);
            fft3_dev(d_rho, d_tmp, nx, ny, nz, true, stream);
            k_j_contract<<<div_up(nk * nao * nao, blk), blk, 0, stream>>>(
                d_ao, d_rho, d_vj, nk, nao, ng, jw);
            k_hermitize<<<nk, 32, 0, stream>>>(d_vj, nk, nao);

            k_zero<<<div_up(static_cast<int>(nmat), blk), blk, 0, stream>>>(
                d_vk, static_cast<int>(nmat));
            for (int k2 = 0; k2 < nk; ++k2) {
                k_aodm<<<div_up(ng, blk), blk, 0, stream>>>(d_ao, d_dm, d_aodm,
                                                            nao, ng, k2);
                for (int k1 = 0; k1 < nk; ++k1) {
                    const double kd[3] = {sys.kpts[k2][0] - sys.kpts[k1][0],
                                          sys.kpts[k2][1] - sys.kpts[k1][1],
                                          sys.kpts[k2][2] - sys.kpts[k1][2]};
                    const int zero =
                        (std::abs(kd[0]) + std::abs(kd[1]) + std::abs(kd[2])
                         < 1e-9)
                            ? 1
                            : 0;
                    k_rho1<<<div_up(ng, blk), blk, 0, stream>>>(
                        d_ao, d_ao, d_coords, kd[0], kd[1], kd[2], zero, d_rho1,
                        nao, ng, k1, k2);
                    const double* coul =
                        d_coul_k
                        + static_cast<std::size_t>((k2 * nk + k1) * ng);
                    for (int ij = 0; ij < nao * nao; ++ij) {
                        Z* slice = d_rho1 + static_cast<std::size_t>(ij) * ng;
                        // fft3_dev uses d_tmp of size ng; slice is its own buffer.
                        fft3_dev(slice, d_tmp, nx, ny, nz, false, stream);
                        k_scale<<<div_up(ng, blk), blk, 0, stream>>>(slice,
                                                                     coul, ng);
                        fft3_dev(slice, d_tmp, nx, ny, nz, true, stream);
                    }
                    k_vrdm<<<div_up(ng, blk), blk, 0, stream>>>(
                        d_rho1, d_aodm, d_coords, kd[0], kd[1], kd[2], zero,
                        d_vrdm, nao, ng);
                    k_k_accum<<<div_up(nao * nao, blk), blk, 0, stream>>>(
                        d_vrdm, d_ao, d_vk, nao, ng, k1, kw);
                }
            }
            k_hermitize<<<nk, 32, 0, stream>>>(d_vk, nk, nao);
            k_vhf<<<div_up(static_cast<int>(nmat), blk), blk, 0, stream>>>(
                d_vj, d_vk, d_vhf, static_cast<int>(nmat));
        };

        double e_prev = 0.0;
        bool conv = false;
        int iters = 0;
        double e_ha = 0.0;
        for (int it = 0; it < max_iter; ++it) {
            build_jk();
            k_trace_energy<<<1, 1, 0, stream>>>(d_dm, d_H, d_vhf, nk, nao, d_e);
            check_cuda(cudaMemcpyAsync(&e_ha, d_e, sizeof(double),
                                       cudaMemcpyDeviceToHost, stream),
                       "D2H e");
            check_cuda(cudaStreamSynchronize(stream), "sync e");
            e_ha += sys.enuc;
            iters = it + 1;
            if (it > 0 && std::abs(e_ha - e_prev) < tol) {
                conv = true;
                break;
            }
            e_prev = e_ha;
            k_fock<<<div_up(static_cast<int>(nmat), blk), blk, 0, stream>>>(
                d_H, d_vj, d_vk, d_F, static_cast<int>(nmat));
            // d_F is the Fock. Diagonalize into d_vj as scratch coefficients,
            // then rebuild dm. Save Fock first: diagonalize_k writes C into d_F.
            check_cuda(cudaMemcpyAsync(d_vhf, d_F, nmat * sizeof(Z),
                                       cudaMemcpyDeviceToDevice, stream),
                       "save F");
            for (int k = 0; k < nk; ++k) {
                diagonalize_k(d_vhf, k, eps[static_cast<std::size_t>(k)],
                              nullptr);
                check_cuda(
                    cudaMemcpyAsync(
                        hC.data() + static_cast<std::size_t>(k) * nao * nao,
                        d_F + static_cast<std::size_t>(k) * nao * nao,
                        sizeof(Z) * nao * nao, cudaMemcpyDeviceToHost, stream),
                    "D2H C2");
            }
            check_cuda(cudaStreamSynchronize(stream), "sync C2");
            all.clear();
            for (int k = 0; k < nk; ++k)
                for (int a = 0; a < nao; ++a)
                    all.push_back({eps[static_cast<std::size_t>(k)]
                                      [static_cast<std::size_t>(a)],
                                   k, a});
            std::stable_sort(
                all.begin(), all.end(),
                [](const Item& x, const Item& y) { return x.e < y.e; });
            std::fill(hocc.begin(), hocc.end(), 0.0);
            for (int i = 0; i < sys.nocc; ++i)
                hocc[static_cast<std::size_t>(
                    all[static_cast<std::size_t>(i)].k * nao
                    + all[static_cast<std::size_t>(i)].a)] = 2.0;
            check_cuda(cudaMemcpyAsync(d_F, hC.data(), nmat * sizeof(Z),
                                       cudaMemcpyHostToDevice, stream),
                       "H2D C2");
            check_cuda(cudaMemcpyAsync(d_occ, hocc.data(),
                                       hocc.size() * sizeof(double),
                                       cudaMemcpyHostToDevice, stream),
                       "H2D occ2");
            for (int k = 0; k < nk; ++k)
                k_dm_from_c<<<1, 32, 0, stream>>>(
                    d_F, d_occ + static_cast<std::size_t>(k) * nao, d_dm, nao,
                    k);
            k_hermitize<<<nk, 32, 0, stream>>>(d_dm, nk, nao);
        }
        cusolverDnDestroy(sol);
        if (!conv)
            throw std::runtime_error(
                "cuda::scf::KRHF: SCF did not converge (missing a device "
                "iteration, or krhf_scf_step did not reach tol)");
        energy_ha_ = e_ha;
        converged_ = true;
        iterations_ = iters;
        if (write) {
            const double e_ev = energy_ha_ * units::kHartreeInEv;
            out->meta.get(kEnergy) = e_ev;
            out->meta.get(kEnergyElectronic) = e_ev;
        }
        guard();
    } catch (...) {
        guard();
        throw;
    }
}

void KRHF::run(cpu::Frame& frame) { solve(frame, true, &frame); }
void KRHF::run(const cpu::Frame& frame) { solve(frame, false, nullptr); }

cpu::scf::GradientResult KRHF::forces(cpu::Frame& /*frame*/) {
    throw std::runtime_error(
        "cuda::scf::KRHF::forces: analytic k-point nuclear forces are not "
        "implemented (missing kernel krhf_nuclear_forces)");
}

} // namespace scf
} // namespace cuda
} // namespace ATV
