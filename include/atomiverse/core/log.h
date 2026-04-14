/// @file core/log.h
/// @brief Project-wide logging via spdlog.
///
/// CUDA-free header -- do NOT include from .cuh or .cu files.
/// Components in .cu files that need logging should delegate to .cpp helpers.
///
/// Usage:
/// @code
///   #include "atomiverse/core/log.h"
///   ATV_INFO("Simulation started: {} atoms, {} steps", n_atoms, n_steps);
///   ATV_DEBUG("Component '{}' built", name);
/// @endcode
///
/// Optional initialization (call once at program start):
/// @code
///   ATV::log_init(spdlog::level::debug);
/// @endcode
#pragma once

#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#endif

#include <chrono>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <memory>

namespace ATV {

using LogClock = std::chrono::steady_clock;

/// @brief Elapsed wall time in milliseconds.
inline double elapsed_ms(
    LogClock::time_point start,
    LogClock::time_point end = LogClock::now())
{
    return std::chrono::duration<double, std::milli>(end - start).count();
}

/// @brief Get the project-wide logger instance.
/// @return Reference to the "atv" named spdlog logger (created on first call).
inline spdlog::logger& logger() {
    static std::shared_ptr<spdlog::logger> instance = spdlog::stdout_color_mt("atv");
    return *instance;
}

/// @brief Initialize logging with Atomiverse defaults.
/// @param level Log level (default: info). Call once at program start; safe to skip.
inline void log_init(spdlog::level::level_enum level = spdlog::level::info) {
    logger().set_level(level);
    logger().set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
    logger().flush_on(spdlog::level::info);
}

} // namespace ATV

// -- Convenience macros (compile to nothing when spdlog strips the level) ----

#define ATV_TRACE(...)    SPDLOG_LOGGER_TRACE(&ATV::logger(), __VA_ARGS__)
#define ATV_DEBUG(...)    SPDLOG_LOGGER_DEBUG(&ATV::logger(), __VA_ARGS__)
#define ATV_INFO(...)     SPDLOG_LOGGER_INFO(&ATV::logger(), __VA_ARGS__)
#define ATV_WARN(...)     SPDLOG_LOGGER_WARN(&ATV::logger(), __VA_ARGS__)
#define ATV_ERROR(...)    SPDLOG_LOGGER_ERROR(&ATV::logger(), __VA_ARGS__)
#define ATV_CRITICAL(...) SPDLOG_LOGGER_CRITICAL(&ATV::logger(), __VA_ARGS__)
