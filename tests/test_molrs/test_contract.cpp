/// @file tests/test_molrs/test_contract.cpp
/// @brief Drift fence for the molrs CXX API / capability contract ATV requires.

#include <cstdint>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "atomiverse/core/molrs_contract.h"
#include "atomiverse/molrs.hpp"

namespace {

/// Frame-block v2 | Element | AM1-BCC — the bits Atomiverse cannot run without.
constexpr std::uint64_t kRequired = 0x7u;
constexpr std::uint64_t kCapTrajectoryWriter = 1u << 3;
constexpr std::uint64_t kCapRegion = 1u << 4;

TEST(MolrsContract, ApiVersionIsOne) {
    EXPECT_EQ(molrs::cxx_api_version(), 1u);
    EXPECT_EQ(ATV::molrs_contract::expected_api_version(), 1u);
}

TEST(MolrsContract, RequiredCapabilitiesPresent) {
    const auto caps = molrs::cxx_api_capabilities();
    EXPECT_EQ(caps & kRequired, kRequired);
    EXPECT_EQ(ATV::molrs_contract::required_capabilities(), kRequired);
    EXPECT_NO_THROW(ATV::molrs_contract::require_compatible());
}

TEST(MolrsContract, FrameSchemaVersionIsTwo) {
    EXPECT_EQ(molrs::frame_schema_version(), 2u);
}

TEST(MolrsContract, ValidateRejectsWrongVersion) {
    EXPECT_THROW(ATV::molrs_contract::validate(2, kRequired),
                 std::runtime_error);
}

TEST(MolrsContract, ValidateRejectsMissingCaps) {
    EXPECT_THROW(ATV::molrs_contract::validate(1, 0), std::runtime_error);
}

// Optional bits: sibling molrs ≥0.14 ships TrajectoryWriter; ≥0.15 adds Region.
// ATV must keep working when they are absent (submodule pin) or present.
TEST(MolrsContract, OptionalCapabilityBitsAreKnownMasks) {
    const auto caps = molrs::cxx_api_capabilities();
    if (caps & kCapTrajectoryWriter) {
        EXPECT_NE(caps & kCapTrajectoryWriter, 0u);
    }
    if (caps & kCapRegion) {
        EXPECT_NE(caps & kCapRegion, 0u);
    }
}

} // namespace
