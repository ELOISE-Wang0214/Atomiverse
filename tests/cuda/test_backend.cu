#include <gtest/gtest.h>
#include <vector>

#include "atomiverse/cuda/backend.cuh"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cuda/transfer.cuh"
#include "atomiverse/facade/stage.h"

using namespace ATV;

namespace {

__global__ void add_one_k(double* fx, double* fy, double* fz, const int* n_ptr) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < *n_ptr) {
        atomicAdd(&fx[idx], double(1));
        atomicAdd(&fy[idx], double(1));
        atomicAdd(&fz[idx], double(1));
    }
}

class AddOne : public ComponentBase {
public:
    std::string name() const override { return "add_one"; }
    void launch(const LaunchContext& ctx) override {
        auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
        auto& frame = ctx.frame;
        int block = 256;
        int grid = (static_cast<int>(ctx.frame.n_atoms()) + block - 1) / block;
        if (grid < 1) grid = 1;
        add_one_k<<<grid, block, 0, stream>>>(frame.atoms.fx, frame.atoms.fy,
                                               frame.atoms.fz, frame.atoms.d_size);
    }
};

__global__ void double_vel_k(double* vx, const int* n_ptr) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < *n_ptr) vx[idx] *= double(2);
}

class DoubleVel : public ComponentBase {
public:
    std::string name() const override { return "double_vel"; }
    Scheduling scheduling() const override { return Scheduling::deferrable; }
    void launch(const LaunchContext& ctx) override {
        auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
        auto& frame = ctx.frame;
        int block = 256;
        int grid = (static_cast<int>(ctx.frame.n_atoms()) + block - 1) / block;
        if (grid < 1) grid = 1;
        double_vel_k<<<grid, block, 0, stream>>>(frame.atoms.vx, frame.atoms.d_size);
    }
};

class FreqDoubler : public ComponentBase {
public:
    std::string name() const override { return "freq_doubler"; }
    int frequency() const override { return 3; }
    void launch(const LaunchContext& ctx) override {
        auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
        auto& frame = ctx.frame;
        int block = 256;
        int grid = (static_cast<int>(ctx.frame.n_atoms()) + block - 1) / block;
        if (grid < 1) grid = 1;
        double_vel_k<<<grid, block, 0, stream>>>(frame.atoms.vx, frame.atoms.d_size);
    }
};

} // anonymous namespace

TEST(Backend, KernelDispatch) {
    cpu::Frame frame(32);
    frame.atoms.mass.fill(double(1.0));
    frame.box = Mat3::cubic(10.0);

    Backend backend;
    backend.init();
    backend.materialize(frame);

    AddOne adder;
    PlanBuilder builder;
    builder.add({0, 0, &adder});
    auto plan = builder.build(static_cast<int>(md::Phase::COUNT));
    backend.load(plan);

    BuildContext bctx{backend.context(), backend.frame()};
    adder.build(bctx);

    backend.step(0);
    backend.sync();

    DestroyContext dctx{backend.context(), backend.frame()};
    adder.destroy(dctx);

    cuda::to_cpu(backend.frame(), frame);
    for (int i = 0; i < 32; ++i)
        EXPECT_DOUBLE_EQ(frame.atoms.fx[i], 1.0);

    backend.shutdown();
}

TEST(Backend, BlockingDeferrable) {
    cpu::Frame frame(32);
    frame.atoms.mass.fill(double(1.0));
    frame.atoms.vx.fill(double(1.0));
    frame.box = Mat3::cubic(10.0);

    Backend backend;
    backend.init();
    backend.materialize(frame);

    AddOne adder;
    DoubleVel doubler;

    PlanBuilder builder;
    builder.add({static_cast<int>(md::Phase::force), 0, &adder});
    builder.add({static_cast<int>(md::Phase::end_of_step), 0, &doubler});
    auto plan = builder.build(static_cast<int>(md::Phase::COUNT));
    backend.load(plan);

    BuildContext bctx{backend.context(), backend.frame()};
    adder.build(bctx);
    doubler.build(bctx);

    backend.step(0);
    backend.sync();

    DestroyContext dctx{backend.context(), backend.frame()};
    adder.destroy(dctx);
    doubler.destroy(dctx);

    cuda::to_cpu(backend.frame(), frame);
    for (int i = 0; i < 32; ++i) {
        EXPECT_DOUBLE_EQ(frame.atoms.fx[i], 1.0);
        EXPECT_DOUBLE_EQ(frame.atoms.vx[i], 2.0);
    }

    backend.shutdown();
}

TEST(Backend, ConditionalExecution) {
    cpu::Frame frame(32);
    frame.atoms.mass.fill(double(1.0));
    frame.atoms.vx.fill(double(1.0));
    frame.box = Mat3::cubic(10.0);

    Backend backend;
    backend.init();
    backend.materialize(frame);

    FreqDoubler freq;
    PlanBuilder builder;
    builder.add({static_cast<int>(md::Phase::end_of_step), 0, &freq});
    auto plan = builder.build(static_cast<int>(md::Phase::COUNT));
    backend.load(plan);

    BuildContext bctx{backend.context(), backend.frame()};
    freq.build(bctx);

    backend.run(6);
    backend.sync();

    DestroyContext dctx{backend.context(), backend.frame()};
    freq.destroy(dctx);

    cuda::to_cpu(backend.frame(), frame);
    for (int i = 0; i < 32; ++i)
        EXPECT_DOUBLE_EQ(frame.atoms.vx[i], 4.0);

    backend.shutdown();
}
