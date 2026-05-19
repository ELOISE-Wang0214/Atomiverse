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
    fftw_plan fwd_ = nullptr;   // r2c: real → complex (Hermitian)
    fftw_plan bwd_ = nullptr;   // c2r: complex (Hermitian) → real
    fftw_plan fwd_cc_ = nullptr; // c2c forward: complex → complex
    fftw_plan bwd_cc_ = nullptr; // c2c backward: complex → complex
    int nx_, ny_, nz_;
    int nr_, nc_; // real size, complex (half) size
    Real* rbuf_ = nullptr;
    fftw_complex* cbuf_ = nullptr;
    fftw_complex* cbuf2_ = nullptr; // for c2c
public:
    FFTPlan3D(int nx, int ny, int nz);
    ~FFTPlan3D();
    FFTPlan3D(const FFTPlan3D&) = delete;
    FFTPlan3D& operator=(const FFTPlan3D&) = delete;

    void r2c(const Real* in, fftw_complex* out);
    void c2r(const fftw_complex* in, Real* out);
    void c2c_fwd(fftw_complex* inout);   // in-place c2c forward
    void c2c_bwd(fftw_complex* inout);   // in-place c2c backward
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
