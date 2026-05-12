/// @file cpu/pw/internal/fft.cpp
/// @brief FFTW3 wrapper implementation.

#include "atomiverse/cpu/pw/internal/fft.h"
#include <cstring>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

FFTPlan3D::FFTPlan3D(int nx, int ny, int nz)
    : nx_(nx), ny_(ny), nz_(nz) {
    nr_ = nx * ny * nz;
    nc_ = (nx / 2 + 1) * ny * nz;
    rbuf_ = fftw_alloc_real(static_cast<size_t>(nr_));
    cbuf_ = fftw_alloc_complex(static_cast<size_t>(nc_));
    if (!rbuf_ || !cbuf_) {
        if (rbuf_) fftw_free(rbuf_);
        if (cbuf_) fftw_free(cbuf_);
        throw std::runtime_error("FFTPlan3D: allocation failed");
    }
    fwd_ = fftw_plan_dft_r2c_3d(nz, ny, nx, rbuf_, cbuf_, FFTW_ESTIMATE);
    bwd_ = fftw_plan_dft_c2r_3d(nz, ny, nx, cbuf_, rbuf_, FFTW_ESTIMATE);
}

FFTPlan3D::~FFTPlan3D() {
    if (fwd_) fftw_destroy_plan(fwd_);
    if (bwd_) fftw_destroy_plan(bwd_);
    if (rbuf_) fftw_free(rbuf_);
    if (cbuf_) fftw_free(cbuf_);
}

void FFTPlan3D::r2c(const Real* in, fftw_complex* out) {
    std::memcpy(rbuf_, in, static_cast<size_t>(nr_) * sizeof(Real));
    fftw_execute_dft_r2c(fwd_, rbuf_, out);
}

void FFTPlan3D::c2r(const fftw_complex* in, Real* out) {
    std::memcpy(cbuf_, in, static_cast<size_t>(nc_) * sizeof(fftw_complex));
    fftw_execute_dft_c2r(bwd_, cbuf_, out);
    std::memcpy(out, rbuf_, static_cast<size_t>(nr_) * sizeof(Real));
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
