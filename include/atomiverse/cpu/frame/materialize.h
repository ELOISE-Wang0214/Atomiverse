/// @file cpu/frame/materialize.h
/// @brief cpu::materialize / cpu::write_back -- molrs.Frame <-> cpu::Frame.
///
/// `molrs.Frame` is the canonical exchange layer; `cpu::Frame` is the CPU
/// backend's execution frame. The materializer is shared infrastructure
/// (sibling of `cpu/transfer.h` at the backend root, per the CLAUDE.md
/// backend-layering rule) and is meant to run once per `run()` -- never
/// inside a step loop.
///
/// `materialize()` reads a fixed field set by known block/column name.
/// Required fields (positions, type) fail loud with a named exception;
/// optional fields (velocities, forces, charge, image, tag, topology)
/// default to zero / empty. `write_back()` pushes positions, velocities,
/// forces, charge, and the system scalars back into the same molrs.Frame.
///
/// @see cpu::Frame for the CPU execution frame
/// @see molrs.hpp for the molrs CXX FFI bridge
#pragma once

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/molrs.hpp"

namespace ATV {
namespace cpu {

/// @brief Materialize a molrs.Frame into a new cpu::Frame.
///
/// Reads the `atoms` block (`x, y, z, vx, vy, vz, fx, fy, fz, mass, type,
/// tag, charge`, plus per-axis image columns `ix, iy, iz`), the `bonds`,
/// `angles`, and `dihedrals` topology blocks, and the simulation box.
///
/// Required fields are the atom positions (`x`, `y`, `z`) and `type`;
/// their absence throws `std::runtime_error` naming the missing field.
/// Every other field is optional and defaults to zero / empty when absent.
///
/// @param src molrs.Frame to read (introspected by block/column name).
/// @return Newly constructed cpu::Frame with all available fields filled.
/// @throws std::runtime_error if a required field is missing.
Frame materialize(const molrs::FrameRef& src);

/// @brief Write a cpu::Frame's results back into an existing molrs.Frame.
///
/// Round-trip field contract (arch-refactor-06-frame; the golden test in
/// tests/cpu/test_materialize.cpp asserts both directions):
///
/// | field | materialize (molrs→cpu) | write_back (cpu→molrs) |
/// |---|---|---|
/// | x/y/z [Å]              | required (throws)       | written |
/// | vx/vy/vz, fx/fy/fz     | optional (zero default) | written |
/// | charge [e]             | optional (zero default) | written |
/// | image (ix/iy/iz)       | optional (i32→int8)     | written (int8→i32) |
/// | tag                    | optional (i32→TagInt)   | written (checked i32 narrowing; overflow throws, no silent truncation) |
/// | type                   | required                | NOT written (per-run immutable input) |
/// | mass [amu]             | optional                | NOT written (per-run immutable input) |
/// | bonds/angles/dihedrals | optional → per-atom     | NOT written (per-run static) |
/// | box                    | optional                | NOT written (per-run static by contract) |
/// | scalars/stress/stress_kinetic | NOT read (write-only channel) | written |
///
/// image evolves during MD (unwrapped-coordinate reconstruction) and tag
/// is the identity anchor after multi-rank migration — both round-trip.
/// Per-frame metadata (step, energy annotations) stays outside this
/// contract: it travels through `write_frame_xyz_meta`
/// (cpu/io/frame_writer.h), not through write_back.
///
/// The system scalars go into dedicated blocks: `energy`, `kinetic`,
/// `conserved` as one-row columns of the `scalars` block, and the two
/// 6-component Voigt tensors as a `value` column of their own `stress` /
/// `stress_kinetic` blocks (a molrs block requires uniform column row
/// counts).
///
/// @param src cpu::Frame whose state is exported.
/// @param dst molrs.Frame to update in place (must be the same frame that
///            was materialized, or any frame with a compatible atom count).
/// @throws std::runtime_error when an atom tag does not fit i32.
void write_back(const Frame& src, molrs::FrameRef& dst);

} // namespace cpu
} // namespace ATV
