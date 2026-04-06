/// @file torch/transfer.cu
/// @brief Transfer function implementations for torch::Frame.

#ifdef ATV_HAS_TORCH

#include "atomiverse/torch/transfer.cuh"

#include <stdexcept>
#include <torch/torch.h>
#include <c10/cuda/CUDAStream.h>

namespace ATV {
namespace torch {

namespace {

/// @brief Wrap a device Real* as a non-owning 1D Tensor.
Tensor wrap_real(Real* ptr, Int n) {
    if (!ptr || n <= 0) return Tensor();
    auto opts = ::torch::TensorOptions().dtype(kRealDtype).device(::torch::kCUDA);
    return ::torch::from_blob(ptr, {static_cast<int64_t>(n)}, opts);
}

/// @brief Wrap a device Int* as a non-owning 1D Tensor.
Tensor wrap_int(Int* ptr, Int n) {
    if (!ptr || n <= 0) return Tensor();
    auto opts = ::torch::TensorOptions().dtype(kIntDtype).device(::torch::kCUDA);
    return ::torch::from_blob(ptr, {static_cast<int64_t>(n)}, opts);
}

/// @brief Wrap a device Real* scalar (1 element) as a non-owning Tensor.
Tensor wrap_scalar(Real* ptr) {
    if (!ptr) return Tensor();
    auto opts = ::torch::TensorOptions().dtype(kRealDtype).device(::torch::kCUDA);
    return ::torch::from_blob(ptr, {1}, opts);
}

/// @brief Extract Real* from a Tensor. Validates CUDA, contiguous, correct dtype.
Real* unwrap_real(const Tensor& t, const char* field_name) {
    if (!t.defined() || t.numel() == 0) return nullptr;
    if (!t.is_cuda())
        throw std::invalid_argument(std::string(field_name) + ": expected CUDA tensor");
    if (t.scalar_type() != kRealDtype)
        throw std::invalid_argument(std::string(field_name) + ": wrong dtype");
    if (!t.is_contiguous())
        throw std::invalid_argument(std::string(field_name) + ": not contiguous");
    return t.data_ptr<Real>();
}

/// @brief Extract Int* from a Tensor.
Int* unwrap_int(const Tensor& t, const char* field_name) {
    if (!t.defined() || t.numel() == 0) return nullptr;
    if (!t.is_cuda())
        throw std::invalid_argument(std::string(field_name) + ": expected CUDA tensor");
    if (t.scalar_type() != kIntDtype)
        throw std::invalid_argument(std::string(field_name) + ": wrong dtype");
    if (!t.is_contiguous())
        throw std::invalid_argument(std::string(field_name) + ": not contiguous");
    return t.data_ptr<Int>();
}

/// @brief Copy a 1D xtensor<T> to a new CUDA Tensor.
template <typename T, ScalarType dtype>
Tensor xtensor_to_cuda(const xt::xtensor<T, 1>& arr) {
    int64_t n = static_cast<int64_t>(arr.size());
    if (n == 0) return Tensor();
    auto opts = ::torch::TensorOptions().dtype(dtype);
    auto cpu_t = ::torch::from_blob(
        const_cast<T*>(arr.data()), {n}, opts);
    return cpu_t.to(::torch::kCUDA);
}

/// @brief Copy a CUDA Tensor to an xtensor<T>.
template <typename T>
void cuda_to_xtensor(const Tensor& t, xt::xtensor<T, 1>& arr) {
    if (!t.defined() || t.numel() == 0) return;
    auto cpu_t = t.to(::torch::kCPU).contiguous();
    Int n = static_cast<Int>(cpu_t.numel());
    std::memcpy(arr.data(), cpu_t.data_ptr<T>(),
                static_cast<size_t>(n) * sizeof(T));
}

} // anonymous namespace

Frame to_torch(const cuda::Frame& src) {
    Frame dst;
    Int na = src.n_atoms();
    dst.n_atoms_ = na;
    dst.n_bonds_ = src.n_bonds();
    dst.n_angles_ = src.n_angles();
    dst.n_dihedrals_ = src.n_dihedrals();

    // Per-atom arrays
    dst.atoms.x = wrap_real(src.atoms.x, na);
    dst.atoms.y = wrap_real(src.atoms.y, na);
    dst.atoms.z = wrap_real(src.atoms.z, na);
    dst.atoms.vx = wrap_real(src.atoms.vx, na);
    dst.atoms.vy = wrap_real(src.atoms.vy, na);
    dst.atoms.vz = wrap_real(src.atoms.vz, na);
    dst.atoms.fx = wrap_real(src.atoms.fx, na);
    dst.atoms.fy = wrap_real(src.atoms.fy, na);
    dst.atoms.fz = wrap_real(src.atoms.fz, na);
    dst.atoms.mass = wrap_real(src.atoms.mass, na);
    dst.atoms.type = wrap_int(src.atoms.type, na);

    // System scalars
    dst.energy = wrap_scalar(src.energy);
    dst.kinetic = wrap_scalar(src.kinetic);
    if (src.stress) {
        auto opts = ::torch::TensorOptions().dtype(kRealDtype).device(::torch::kCUDA);
        dst.stress = ::torch::from_blob(src.stress, {6}, opts);
    }
    if (src.box) {
        auto opts = ::torch::TensorOptions().dtype(kRealDtype).device(::torch::kCUDA);
        dst.box = ::torch::from_blob(src.box, {9}, opts);
    }

    // Topology
    Int nb = src.n_bonds();
    dst.bonds.atomi = wrap_int(src.bonds.atomi, nb);
    dst.bonds.atomj = wrap_int(src.bonds.atomj, nb);
    dst.bonds.type = wrap_int(src.bonds.type, nb);

    Int nang = src.n_angles();
    dst.angles.atomi = wrap_int(src.angles.atomi, nang);
    dst.angles.atomj = wrap_int(src.angles.atomj, nang);
    dst.angles.atomk = wrap_int(src.angles.atomk, nang);
    dst.angles.type = wrap_int(src.angles.type, nang);

    Int ndih = src.n_dihedrals();
    dst.dihedrals.atomi = wrap_int(src.dihedrals.atomi, ndih);
    dst.dihedrals.atomj = wrap_int(src.dihedrals.atomj, ndih);
    dst.dihedrals.atomk = wrap_int(src.dihedrals.atomk, ndih);
    dst.dihedrals.atoml = wrap_int(src.dihedrals.atoml, ndih);
    dst.dihedrals.type = wrap_int(src.dihedrals.type, ndih);

    // Pair list (may not be populated yet)
    // PairList n_pairs is on device -- we don't D2H here; caller sets n_pairs_ if known.

    return dst;
}

void to_cuda(const Frame& src, cuda::Frame& dst) {
    dst.n_atoms_ = src.n_atoms_;
    dst.n_bonds_ = src.n_bonds_;
    dst.n_angles_ = src.n_angles_;
    dst.n_dihedrals_ = src.n_dihedrals_;

    dst.atoms.x = unwrap_real(src.atoms.x, "x");
    dst.atoms.y = unwrap_real(src.atoms.y, "y");
    dst.atoms.z = unwrap_real(src.atoms.z, "z");
    dst.atoms.vx = unwrap_real(src.atoms.vx, "vx");
    dst.atoms.vy = unwrap_real(src.atoms.vy, "vy");
    dst.atoms.vz = unwrap_real(src.atoms.vz, "vz");
    dst.atoms.fx = unwrap_real(src.atoms.fx, "fx");
    dst.atoms.fy = unwrap_real(src.atoms.fy, "fy");
    dst.atoms.fz = unwrap_real(src.atoms.fz, "fz");
    dst.atoms.mass = unwrap_real(src.atoms.mass, "mass");
    dst.atoms.type = unwrap_int(src.atoms.type, "type");

    dst.energy = unwrap_real(src.energy, "energy");
    dst.kinetic = unwrap_real(src.kinetic, "kinetic");
    dst.stress = unwrap_real(src.stress, "stress");
    dst.box = unwrap_real(src.box, "box");

    dst.bonds.atomi = unwrap_int(src.bonds.atomi, "bonds.atomi");
    dst.bonds.atomj = unwrap_int(src.bonds.atomj, "bonds.atomj");
    dst.bonds.type = unwrap_int(src.bonds.type, "bonds.type");

    dst.angles.atomi = unwrap_int(src.angles.atomi, "angles.atomi");
    dst.angles.atomj = unwrap_int(src.angles.atomj, "angles.atomj");
    dst.angles.atomk = unwrap_int(src.angles.atomk, "angles.atomk");
    dst.angles.type = unwrap_int(src.angles.type, "angles.type");

    dst.dihedrals.atomi = unwrap_int(src.dihedrals.atomi, "dihedrals.atomi");
    dst.dihedrals.atomj = unwrap_int(src.dihedrals.atomj, "dihedrals.atomj");
    dst.dihedrals.atomk = unwrap_int(src.dihedrals.atomk, "dihedrals.atomk");
    dst.dihedrals.atoml = unwrap_int(src.dihedrals.atoml, "dihedrals.atoml");
    dst.dihedrals.type = unwrap_int(src.dihedrals.type, "dihedrals.type");

    dst.pairs.i = unwrap_int(src.pair_i, "pair_i");
    dst.pairs.j = unwrap_int(src.pair_j, "pair_j");
    dst.pairs.dx = unwrap_real(src.pair_dx, "pair_dx");
    dst.pairs.dy = unwrap_real(src.pair_dy, "pair_dy");
    dst.pairs.dz = unwrap_real(src.pair_dz, "pair_dz");
    dst.pairs.dist = unwrap_real(src.pair_dist, "pair_dist");
}

void to_cpu(const Frame& src, cpu::Frame& dst) {
    cuda_to_xtensor<Real>(src.atoms.x, dst.atoms.x);
    cuda_to_xtensor<Real>(src.atoms.y, dst.atoms.y);
    cuda_to_xtensor<Real>(src.atoms.z, dst.atoms.z);
    cuda_to_xtensor<Real>(src.atoms.vx, dst.atoms.vx);
    cuda_to_xtensor<Real>(src.atoms.vy, dst.atoms.vy);
    cuda_to_xtensor<Real>(src.atoms.vz, dst.atoms.vz);
    cuda_to_xtensor<Real>(src.atoms.fx, dst.atoms.fx);
    cuda_to_xtensor<Real>(src.atoms.fy, dst.atoms.fy);
    cuda_to_xtensor<Real>(src.atoms.fz, dst.atoms.fz);
    cuda_to_xtensor<Real>(src.atoms.mass, dst.atoms.mass);
    cuda_to_xtensor<Int>(src.atoms.type, dst.atoms.type);

    if (src.energy.defined() && src.energy.numel() > 0) {
        auto cpu_e = src.energy.to(::torch::kCPU);
        dst.energy = cpu_e.item<Real>();
    }
    if (src.kinetic.defined() && src.kinetic.numel() > 0) {
        auto cpu_k = src.kinetic.to(::torch::kCPU);
        dst.kinetic = cpu_k.item<Real>();
    }
    if (src.stress.defined() && src.stress.numel() == 6) {
        auto cpu_s = src.stress.to(::torch::kCPU).contiguous();
        std::memcpy(dst.stress, cpu_s.data_ptr<Real>(), 6 * sizeof(Real));
    }
    if (src.box.defined() && src.box.numel() == 9) {
        auto cpu_b = src.box.to(::torch::kCPU).contiguous();
        auto* bp = cpu_b.data_ptr<Real>();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                dst.box(i, j) = bp[i * 3 + j];
    }
}

Frame from_cpu(const cpu::Frame& src) {
    Frame dst;
    Int na = src.atoms.size();
    dst.n_atoms_ = na;

    dst.atoms.x = xtensor_to_cuda<Real, kRealDtype>(src.atoms.x);
    dst.atoms.y = xtensor_to_cuda<Real, kRealDtype>(src.atoms.y);
    dst.atoms.z = xtensor_to_cuda<Real, kRealDtype>(src.atoms.z);
    dst.atoms.vx = xtensor_to_cuda<Real, kRealDtype>(src.atoms.vx);
    dst.atoms.vy = xtensor_to_cuda<Real, kRealDtype>(src.atoms.vy);
    dst.atoms.vz = xtensor_to_cuda<Real, kRealDtype>(src.atoms.vz);
    dst.atoms.fx = xtensor_to_cuda<Real, kRealDtype>(src.atoms.fx);
    dst.atoms.fy = xtensor_to_cuda<Real, kRealDtype>(src.atoms.fy);
    dst.atoms.fz = xtensor_to_cuda<Real, kRealDtype>(src.atoms.fz);
    dst.atoms.mass = xtensor_to_cuda<Real, kRealDtype>(src.atoms.mass);
    dst.atoms.type = xtensor_to_cuda<Int, kIntDtype>(src.atoms.type);

    // System scalars
    auto real_opts = ::torch::TensorOptions().dtype(kRealDtype).device(::torch::kCUDA);
    dst.energy = ::torch::zeros({1}, real_opts);
    dst.kinetic = ::torch::zeros({1}, real_opts);
    dst.stress = ::torch::zeros({6}, real_opts);

    // Box: copy Mat3 (9 Reals, row-major)
    Real box_data[9] = {
        src.box(0,0), src.box(0,1), src.box(0,2),
        src.box(1,0), src.box(1,1), src.box(1,2),
        src.box(2,0), src.box(2,1), src.box(2,2)
    };
    auto cpu_opts = ::torch::TensorOptions().dtype(kRealDtype);
    dst.box = ::torch::from_blob(box_data, {9}, cpu_opts).to(::torch::kCUDA);

    // Topology
    if (src.bonds.atomi.size() > 0) {
        dst.n_bonds_ = static_cast<Int>(src.bonds.atomi.size());
        dst.bonds.atomi = xtensor_to_cuda<Int, kIntDtype>(src.bonds.atomi);
        dst.bonds.atomj = xtensor_to_cuda<Int, kIntDtype>(src.bonds.atomj);
        dst.bonds.type = xtensor_to_cuda<Int, kIntDtype>(src.bonds.type);
    }
    if (src.angles.atomi.size() > 0) {
        dst.n_angles_ = static_cast<Int>(src.angles.atomi.size());
        dst.angles.atomi = xtensor_to_cuda<Int, kIntDtype>(src.angles.atomi);
        dst.angles.atomj = xtensor_to_cuda<Int, kIntDtype>(src.angles.atomj);
        dst.angles.atomk = xtensor_to_cuda<Int, kIntDtype>(src.angles.atomk);
        dst.angles.type = xtensor_to_cuda<Int, kIntDtype>(src.angles.type);
    }
    if (src.dihedrals.atomi.size() > 0) {
        dst.n_dihedrals_ = static_cast<Int>(src.dihedrals.atomi.size());
        dst.dihedrals.atomi = xtensor_to_cuda<Int, kIntDtype>(src.dihedrals.atomi);
        dst.dihedrals.atomj = xtensor_to_cuda<Int, kIntDtype>(src.dihedrals.atomj);
        dst.dihedrals.atomk = xtensor_to_cuda<Int, kIntDtype>(src.dihedrals.atomk);
        dst.dihedrals.atoml = xtensor_to_cuda<Int, kIntDtype>(src.dihedrals.atoml);
        dst.dihedrals.type = xtensor_to_cuda<Int, kIntDtype>(src.dihedrals.type);
    }

    return dst;
}

void destroy(Frame& frame) {
    frame = Frame{};
}

} // namespace torch
} // namespace ATV

#endif // ATV_HAS_TORCH
