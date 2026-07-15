/// @file cpu/task/task_codec.cpp
/// @brief TaskSpec byte-buffer pack/unpack (wire format in task_codec.h).

#include "atomiverse/cpu/task/task_codec.h"

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace ATV {
namespace cpu {
namespace task {

namespace {

constexpr char kMagic[4] = {'A', 'T', 'V', 'T'};

// ── Writer ────────────────────────────────────────────────────────────────

struct Writer {
    std::vector<std::byte>& out;

    template <typename T>
    void raw(const T& v) {
        static_assert(std::is_trivially_copyable_v<T>, "wire scalar");
        const std::size_t at = out.size();
        out.resize(at + sizeof(T));
        std::memcpy(out.data() + at, &v, sizeof(T));
    }
    void u8(std::uint8_t v) { raw(v); }
    void u16(std::uint16_t v) { raw(v); }
    void u32(std::uint32_t v) { raw(v); }
    void i32(std::int32_t v) { raw(v); }
    void u64(std::uint64_t v) { raw(v); }
    void i64(std::int64_t v) { raw(v); }
    void f64(double v) { raw(v); }
    void vec(const std::vector<double>& v) {
        u64(v.size());
        for (double x : v)
            f64(x);
    }
};

// ── Reader (every read is bounds-checked; truncation throws) ─────────────

struct Reader {
    const std::byte* data;
    std::size_t size;
    std::size_t at = 0;

    template <typename T>
    T raw() {
        static_assert(std::is_trivially_copyable_v<T>, "wire scalar");
        if (at + sizeof(T) > size)
            throw std::runtime_error(
                "task_codec: truncated buffer (need " +
                std::to_string(at + sizeof(T)) + " bytes, have " +
                std::to_string(size) + ")");
        T v;
        std::memcpy(&v, data + at, sizeof(T));
        at += sizeof(T);
        return v;
    }
    std::uint8_t u8() { return raw<std::uint8_t>(); }
    std::uint16_t u16() { return raw<std::uint16_t>(); }
    std::uint32_t u32() { return raw<std::uint32_t>(); }
    std::int32_t i32() { return raw<std::int32_t>(); }
    std::uint64_t u64() { return raw<std::uint64_t>(); }
    std::int64_t i64() { return raw<std::int64_t>(); }
    double f64() { return raw<double>(); }
    std::vector<double> vec() {
        const std::uint64_t n = u64();
        std::vector<double> v;
        v.reserve(static_cast<std::size_t>(n));
        for (std::uint64_t i = 0; i < n; ++i)
            v.push_back(f64());
        return v;
    }
};

void write_params(Writer& w, const ComponentParams& params) {
    std::visit(
        [&](const auto& p) {
            using P = std::decay_t<decltype(p)>;
            if constexpr (std::is_same_v<P, PairLJCutParams>) {
                w.f64(p.epsilon);
                w.f64(p.sigma);
                w.f64(p.cutoff);
            } else if constexpr (std::is_same_v<P, NeighborListParams>) {
                w.f64(p.cutoff);
                w.f64(p.skin);
            } else if constexpr (std::is_same_v<P, BondHarmonicParams>) {
                w.vec(p.K);
                w.vec(p.r0);
            } else if constexpr (std::is_same_v<P, AngleHarmonicParams>) {
                w.vec(p.K);
                w.vec(p.theta0);
            } else if constexpr (std::is_same_v<P,
                                                DihedralMultiHarmonicParams>) {
                w.u64(p.Cn.size());
                for (const auto& inner : p.Cn)
                    w.vec(inner);
            } else if constexpr (std::is_same_v<P, VelocityVerletParams>) {
                w.f64(p.dt);
            } else if constexpr (std::is_same_v<P, InitVelocityParams>) {
                w.f64(p.temperature);
                w.u32(p.seed);
            } else if constexpr (std::is_same_v<P, NoseHooverParams>) {
                w.f64(p.T_start);
                w.f64(p.tau_T);
                w.f64(p.dt);
            } else if constexpr (std::is_same_v<P, ThermoParams>) {
                w.i32(p.freq);
            } else if constexpr (std::is_same_v<P, AM1EvaluatorParams>) {
                w.i32(p.charge);
                w.i32(p.multiplicity);
                w.i32(p.max_iter);
                w.f64(p.density_tol);
                w.f64(p.energy_tol);
            }
        },
        params);
}

ComponentParams read_params(Reader& r, std::uint16_t component_id) {
    switch (component_id) {
    case 0: {
        PairLJCutParams p{};
        p.epsilon = r.f64();
        p.sigma = r.f64();
        p.cutoff = r.f64();
        return p;
    }
    case 1: {
        NeighborListParams p{};
        p.cutoff = r.f64();
        p.skin = r.f64();
        return p;
    }
    case 2: {
        BondHarmonicParams p{};
        p.K = r.vec();
        p.r0 = r.vec();
        return p;
    }
    case 3: {
        AngleHarmonicParams p{};
        p.K = r.vec();
        p.theta0 = r.vec();
        return p;
    }
    case 4: {
        DihedralMultiHarmonicParams p{};
        const std::uint64_t n = r.u64();
        p.Cn.reserve(static_cast<std::size_t>(n));
        for (std::uint64_t i = 0; i < n; ++i)
            p.Cn.push_back(r.vec());
        return p;
    }
    case 5: {
        VelocityVerletParams p{};
        p.dt = r.f64();
        return p;
    }
    case 6: {
        InitVelocityParams p{};
        p.temperature = r.f64();
        p.seed = r.u32();
        return p;
    }
    case 7: {
        NoseHooverParams p{};
        p.T_start = r.f64();
        p.tau_T = r.f64();
        p.dt = r.f64();
        return p;
    }
    case 8: {
        ThermoParams p{};
        p.freq = r.i32();
        return p;
    }
    case 9: {
        AM1EvaluatorParams p{};
        p.charge = r.i32();
        p.multiplicity = r.i32();
        p.max_iter = r.i32();
        p.density_tol = r.f64();
        p.energy_tol = r.f64();
        return p;
    }
    default:
        throw std::runtime_error(
            "task_codec: unknown component_id " +
            std::to_string(component_id) + " (enumeration ends at " +
            std::to_string(std::variant_size_v<ComponentParams> - 1) + ")");
    }
}

} // namespace

std::vector<std::byte> pack(const TaskSpec& spec) {
    std::vector<std::byte> out;
    Writer w{out};
    w.raw(kMagic);
    w.u32(TaskSpec::kWireVersion);
    w.u8(static_cast<std::uint8_t>(spec.kind()));
    w.u32(static_cast<std::uint32_t>(spec.entries().size()));
    for (const auto& e : spec.entries()) {
        w.u16(static_cast<std::uint16_t>(e.params.index()));
        w.u8(e.binding ? 1 : 0);
        if (e.binding) {
            w.u32(static_cast<std::uint32_t>(ATV::slot_of(e.binding->phase)));
            w.i32(e.binding->order);
        }
        write_params(w, e.params);
    }
    w.i64(spec.n_steps());
    return out;
}

TaskSpec unpack(const std::byte* data, std::size_t size) {
    Reader r{data, size};

    char magic[4];
    magic[0] = static_cast<char>(r.u8());
    magic[1] = static_cast<char>(r.u8());
    magic[2] = static_cast<char>(r.u8());
    magic[3] = static_cast<char>(r.u8());
    if (std::memcmp(magic, kMagic, sizeof(kMagic)) != 0)
        throw std::runtime_error("task_codec: bad magic (not a TaskSpec "
                                 "buffer)");

    const std::uint32_t version = r.u32();
    if (version != TaskSpec::kWireVersion)
        throw std::runtime_error(
            "task_codec: wire version " + std::to_string(version) +
            " does not match kWireVersion " +
            std::to_string(TaskSpec::kWireVersion));

    const std::uint8_t kind_byte = r.u8();
    if (kind_byte > static_cast<std::uint8_t>(TaskKind::SinglePoint))
        throw std::runtime_error(
            "task_codec: unknown task_kind " +
            std::to_string(static_cast<int>(kind_byte)));
    const TaskKind kind = static_cast<TaskKind>(kind_byte);

    const std::uint32_t n_entries = r.u32();
    std::vector<ComponentEntry> entries;
    entries.reserve(n_entries);
    for (std::uint32_t i = 0; i < n_entries; ++i) {
        const std::uint16_t id = r.u16();
        std::optional<ExplicitBinding> binding;
        if (r.u8() != 0) {
            const std::uint32_t slot = r.u32();
            if (slot >= static_cast<std::uint32_t>(
                            ATV::slot_count<ATV::md::Phase>()))
                throw std::runtime_error(
                    "task_codec: phase slot " + std::to_string(slot) +
                    " out of range for md::Phase");
            const std::int32_t order = r.i32();
            binding = ExplicitBinding{static_cast<ATV::md::Phase>(slot),
                                      order};
        }
        entries.push_back(ComponentEntry{read_params(r, id), binding});
    }
    const std::int64_t n_steps = r.i64();

    // Rebuild through the validating constructors so wire-level task-kind
    // invariant violations cannot yield a half-valid spec.
    if (kind == TaskKind::MD)
        return TaskSpec(std::move(entries), n_steps);

    if (entries.size() != 1)
        throw std::runtime_error(
            "task_codec: SinglePoint task must carry exactly one entry, "
            "got " +
            std::to_string(entries.size()));
    const auto& only = entries.front();
    if (!std::holds_alternative<AM1EvaluatorParams>(only.params))
        throw std::runtime_error(
            std::string("task_codec: SinglePoint task carries ") +
            component_name(only.params.index()) +
            ", which is not an evaluator parameter set");
    if (only.binding)
        throw std::runtime_error(
            "task_codec: SinglePoint entry must not carry an explicit "
            "binding (the materializer pins evaluate::Phase::evaluate)");
    return TaskSpec(std::get<AM1EvaluatorParams>(only.params));
}

} // namespace task
} // namespace cpu
} // namespace ATV
