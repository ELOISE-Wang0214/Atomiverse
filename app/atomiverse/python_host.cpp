/// @file app/atomiverse/python_host.cpp
/// @brief Embedded-CPython host (contract in python_host.h).

#include "python_host.h"

#include <Python.h>

// POSIX getpid/setenv: the launcher process model is Linux-only by charter.
#include <unistd.h>

#include <cstdlib>
#include <stdexcept>
#include <string>
#include <utility>

#ifndef ATV_LAUNCHER_BUILD_PYTHON
#error "ATV_LAUNCHER_BUILD_PYTHON must carry the build-time interpreter path"
#endif

namespace ATV {
namespace app {

PythonHost::PythonHost(std::string script, int argc, char** argv)
    : script_(std::move(script)), argc_(argc), argv_(argv) {
    // The marker MUST land before Py_Initialize: the bindings' LauncherGate
    // (in the extension module's copy of the library) reads the environment
    // lazily at the first run() — the value is OUR pid, so a child process
    // inheriting the environment stays silently serial (PID mismatch).
    const std::string pid = std::to_string(static_cast<long>(::getpid()));
    ::setenv("ATV_LAUNCHER", pid.c_str(), /*overwrite=*/1);
}

int PythonHost::run() {
    PyConfig config;
    PyConfig_InitPythonConfig(&config);

    // 04's build-time interpreter identity: PyConfig's getpath computes
    // sys.prefix / site-packages from this executable's environment
    // (pyvenv.cfg honoured), which is exactly the venv the wheel was
    // installed into. Runtime discovery replaces this in chain link 05.
    PyStatus status = PyConfig_SetBytesString(
        &config, &config.executable, ATV_LAUNCHER_BUILD_PYTHON);
    if (PyStatus_Exception(status)) {
        PyConfig_Clear(&config);
        throw std::runtime_error("atomiverse: PyConfig executable setup "
                                 "failed");
    }
    status = PyConfig_SetBytesString(&config, &config.run_filename,
                                     script_.c_str());
    if (PyStatus_Exception(status)) {
        PyConfig_Clear(&config);
        throw std::runtime_error("atomiverse: PyConfig script setup failed");
    }
    // The script's own argv (sys.argv = [script, args...]); parse_argv
    // stays off — the launcher consumed its argv already.
    config.parse_argv = 0;
    status = PyConfig_SetBytesArgv(&config, argc_, argv_);
    if (PyStatus_Exception(status)) {
        PyConfig_Clear(&config);
        throw std::runtime_error("atomiverse: PyConfig argv setup failed");
    }

    status = Py_InitializeFromConfig(&config);
    PyConfig_Clear(&config);
    if (PyStatus_Exception(status)) {
        if (PyStatus_IsExit(status))
            return status.exitcode;
        throw std::runtime_error(
            std::string("atomiverse: Python initialization failed: ") +
            (status.err_msg != nullptr ? status.err_msg : "unknown"));
    }

    // Runs the script as __main__, prints any traceback, finalizes the
    // interpreter, and returns the process exit code. Driver destructors
    // inside the script run during interpreter teardown here — before the
    // caller broadcasts SHUTDOWN and finalizes MPI (the Backend::shutdown-
    // before-MPI_Finalize ordering holds by construction).
    return Py_RunMain();
}

} // namespace app
} // namespace ATV
