/// @file cpu/pw/internal/fft.h
/// @brief RAII wrapper around FFTW3 3D transforms.
#pragma once

#include <fftw3.h>
#include <stdexcept>
#include "atomiverse/core/real.h"

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

class FFTPlan3D {
    fftw_plan fwd_ = nullptr;
    fftw_plan bwd_ = nullptr;
    int nx_, ny_, nz_;
    int nr_, nc_; // real size, complex size (nx/2+1)*ny*nz
    Real* rbuf_ = nullptr;
    fftw_complex* cbuf_ = nullptr;
public:
    FFTPlan3D(int nx, int ny, int nz);
    ~FFTPlan3D();
    FFTPlan3D(const FFTPlan3D&) = delete;
    FFTPlan3D& operator=(const FFTPlan3D&) = delete;

    void r2c(const Real* in, fftw_complex* out);
    void c2r(const fftw_complex* in, Real* out);
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }
    int nr() const { return nr_; }
    int nc() const { return nc_; }
};

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
