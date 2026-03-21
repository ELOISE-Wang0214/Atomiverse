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
    dst.x = wrap_real(src.x, na);
    dst.y = wrap_real(src.y, na);
    dst.z = wrap_real(src.z, na);
    dst.vx = wrap_real(src.vx, na);
    dst.vy = wrap_real(src.vy, na);
    dst.vz = wrap_real(src.vz, na);
    dst.fx = wrap_real(src.fx, na);
    dst.fy = wrap_real(src.fy, na);
    dst.fz = wrap_real(src.fz, na);
    dst.mass = wrap_real(src.mass, na);
    dst.type_id = wrap_int(src.type_id, na);

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
    dst.bond_i = wrap_int(src.bonds.i, nb);
    dst.bond_j = wrap_int(src.bonds.j, nb);
    dst.bond_type = wrap_int(src.bonds.type, nb);

    Int nang = src.n_angles();
    dst.angle_i = wrap_int(src.angles.i, nang);
    dst.angle_j = wrap_int(src.angles.j, nang);
    dst.angle_k = wrap_int(src.angles.k, nang);
    dst.angle_type = wrap_int(src.angles.type, nang);

    Int ndih = src.n_dihedrals();
    dst.dihedral_i = wrap_int(src.dihedrals.i, ndih);
    dst.dihedral_j = wrap_int(src.dihedrals.j, ndih);
    dst.dihedral_k = wrap_int(src.dihedrals.k, ndih);
    dst.dihedral_l = wrap_int(src.dihedrals.l, ndih);
    dst.dihedral_type = wrap_int(src.dihedrals.type, ndih);

    // Pair list (may not be populated yet)
    // PairList n_pairs is on device -- we don't D2H here; caller sets n_pairs_ if known.

    return dst;
}

void to_cuda(const Frame& src, cuda::Frame& dst) {
    dst.n_atoms_ = src.n_atoms_;
    dst.n_bonds_ = src.n_bonds_;
    dst.n_angles_ = src.n_angles_;
    dst.n_dihedrals_ = src.n_dihedrals_;

    dst.x = unwrap_real(src.x, "x");
    dst.y = unwrap_real(src.y, "y");
    dst.z = unwrap_real(src.z, "z");
    dst.vx = unwrap_real(src.vx, "vx");
    dst.vy = unwrap_real(src.vy, "vy");
    dst.vz = unwrap_real(src.vz, "vz");
    dst.fx = unwrap_real(src.fx, "fx");
    dst.fy = unwrap_real(src.fy, "fy");
    dst.fz = unwrap_real(src.fz, "fz");
    dst.mass = unwrap_real(src.mass, "mass");
    dst.type_id = unwrap_int(src.type_id, "type_id");

    dst.energy = unwrap_real(src.energy, "energy");
    dst.kinetic = unwrap_real(src.kinetic, "kinetic");
    dst.stress = unwrap_real(src.stress, "stress");
    dst.box = unwrap_real(src.box, "box");

    dst.bonds.i = unwrap_int(src.bond_i, "bond_i");
    dst.bonds.j = unwrap_int(src.bond_j, "bond_j");
    dst.bonds.type = unwrap_int(src.bond_type, "bond_type");

    dst.angles.i = unwrap_int(src.angle_i, "angle_i");
    dst.angles.j = unwrap_int(src.angle_j, "angle_j");
    dst.angles.k = unwrap_int(src.angle_k, "angle_k");
    dst.angles.type = unwrap_int(src.angle_type, "angle_type");

    dst.dihedrals.i = unwrap_int(src.dihedral_i, "dihedral_i");
    dst.dihedrals.j = unwrap_int(src.dihedral_j, "dihedral_j");
    dst.dihedrals.k = unwrap_int(src.dihedral_k, "dihedral_k");
    dst.dihedrals.l = unwrap_int(src.dihedral_l, "dihedral_l");
    dst.dihedrals.type = unwrap_int(src.dihedral_type, "dihedral_type");

    dst.pairs.i = unwrap_int(src.pair_i, "pair_i");
    dst.pairs.j = unwrap_int(src.pair_j, "pair_j");
    dst.pairs.dx = unwrap_real(src.pair_dx, "pair_dx");
    dst.pairs.dy = unwrap_real(src.pair_dy, "pair_dy");
    dst.pairs.dz = unwrap_real(src.pair_dz, "pair_dz");
    dst.pairs.dist = unwrap_real(src.pair_dist, "pair_dist");
}

void to_cpu(const Frame& src, cpu::Frame& dst) {
    cuda_to_xtensor<Real>(src.x, dst.x);
    cuda_to_xtensor<Real>(src.y, dst.y);
    cuda_to_xtensor<Real>(src.z, dst.z);
    cuda_to_xtensor<Real>(src.vx, dst.vx);
    cuda_to_xtensor<Real>(src.vy, dst.vy);
    cuda_to_xtensor<Real>(src.vz, dst.vz);
    cuda_to_xtensor<Real>(src.fx, dst.fx);
    cuda_to_xtensor<Real>(src.fy, dst.fy);
    cuda_to_xtensor<Real>(src.fz, dst.fz);
    cuda_to_xtensor<Real>(src.mass, dst.mass);
    cuda_to_xtensor<Int>(src.type_id, dst.type_id);

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
    Int na = src.n_atoms();
    dst.n_atoms_ = na;

    dst.x = xtensor_to_cuda<Real, kRealDtype>(src.x);
    dst.y = xtensor_to_cuda<Real, kRealDtype>(src.y);
    dst.z = xtensor_to_cuda<Real, kRealDtype>(src.z);
    dst.vx = xtensor_to_cuda<Real, kRealDtype>(src.vx);
    dst.vy = xtensor_to_cuda<Real, kRealDtype>(src.vy);
    dst.vz = xtensor_to_cuda<Real, kRealDtype>(src.vz);
    dst.fx = xtensor_to_cuda<Real, kRealDtype>(src.fx);
    dst.fy = xtensor_to_cuda<Real, kRealDtype>(src.fy);
    dst.fz = xtensor_to_cuda<Real, kRealDtype>(src.fz);
    dst.mass = xtensor_to_cuda<Real, kRealDtype>(src.mass);
    dst.type_id = xtensor_to_cuda<Int, kIntDtype>(src.type_id);

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
    if (src.bond_i.size() > 0) {
        dst.n_bonds_ = static_cast<Int>(src.bond_i.size());
        dst.bond_i = xtensor_to_cuda<Int, kIntDtype>(src.bond_i);
        dst.bond_j = xtensor_to_cuda<Int, kIntDtype>(src.bond_j);
        dst.bond_type = xtensor_to_cuda<Int, kIntDtype>(src.bond_type);
    }
    if (src.angle_i.size() > 0) {
        dst.n_angles_ = static_cast<Int>(src.angle_i.size());
        dst.angle_i = xtensor_to_cuda<Int, kIntDtype>(src.angle_i);
        dst.angle_j = xtensor_to_cuda<Int, kIntDtype>(src.angle_j);
        dst.angle_k = xtensor_to_cuda<Int, kIntDtype>(src.angle_k);
        dst.angle_type = xtensor_to_cuda<Int, kIntDtype>(src.angle_type);
    }
    if (src.dihedral_i.size() > 0) {
        dst.n_dihedrals_ = static_cast<Int>(src.dihedral_i.size());
        dst.dihedral_i = xtensor_to_cuda<Int, kIntDtype>(src.dihedral_i);
        dst.dihedral_j = xtensor_to_cuda<Int, kIntDtype>(src.dihedral_j);
        dst.dihedral_k = xtensor_to_cuda<Int, kIntDtype>(src.dihedral_k);
        dst.dihedral_l = xtensor_to_cuda<Int, kIntDtype>(src.dihedral_l);
        dst.dihedral_type = xtensor_to_cuda<Int, kIntDtype>(src.dihedral_type);
    }

    return dst;
}

void destroy(Frame& frame) {
    frame = Frame{};
}

} // namespace torch
} // namespace ATV

#endif // ATV_HAS_TORCH
