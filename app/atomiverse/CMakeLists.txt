# The `atomiverse` MPI launcher: every rank native, rank 0 embeds CPython.
# ONLY this executable links libpython — libatomiverse and _atomiverse.so
# stay untouched (the wheel's serial `python atv.py` path is unaffected).
#
# The CMake target is `atomiverse_launcher` (the library target already owns
# the name `atomiverse`); the on-disk binary is still `atomiverse`.

find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter Development.Embed)

add_executable(atomiverse_launcher
    main.cpp
    python_host.cpp
)
set_target_properties(atomiverse_launcher PROPERTIES
    OUTPUT_NAME atomiverse
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    # 05-packaging: conda/venv layouts keep libpython next to ../lib.
    INSTALL_RPATH "$ORIGIN/../lib")
target_link_libraries(atomiverse_launcher PRIVATE
    atomiverse
    Python3::Python
    MPI::MPI_CXX)
# 04's build-time interpreter identity (runtime discovery is chain link 05):
# PyConfig.executable is set to the interpreter this build was configured
# against, so its environment (venv/conda) resolves stdlib + site-packages.
target_compile_definitions(atomiverse_launcher PRIVATE
    ATV_LAUNCHER_BUILD_PYTHON="${Python3_EXECUTABLE}")
