# Exact Atomiverse -> molrs CXX consumer contract.
#
# This is Atomiverse's single CMake-readable source of truth. The actual molrs
# version lives in molrs-cxxapi/CXX_API_VERSION and is checked even when a
# developer supplies ATV_MOLRS_PATH.
#
# Pin policy (cpu-breadth-01 ac-002):
#   - ATV_MOLRS_CXX_API_VERSION must match molrs-cxxapi/CXX_API_VERSION exactly.
#   - ATV_MOLRS_GIT_REF is the CI checkout ref of MolCrafts/molrs (tag or SHA).
#     Bump both together when promoting a new molrs release. Local ATV_MOLRS_PATH
#     still cannot bypass the CXX API version check.
set(ATV_MOLRS_CXX_API_VERSION 1)
set(ATV_MOLRS_CXX_REQUIRED_CAPABILITIES 7) # Frame-block v2 | Element | AM1-BCC
set(ATV_MOLRS_GIT_REF "v0.8.2" CACHE STRING
    "Pinned MolCrafts/molrs git tag or SHA for CI checkout")

function(atv_validate_molrs_source_contract molrs_root)
    set(version_file "${molrs_root}/molrs-cxxapi/CXX_API_VERSION")
    if(NOT CMAKE_SCRIPT_MODE_FILE)
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                     "${version_file}")
    endif()
    if(NOT EXISTS "${version_file}")
        message(FATAL_ERROR
            "Atomiverse requires molrs CXX API v${ATV_MOLRS_CXX_API_VERSION}, "
            "but ${version_file} does not exist. ATV_MOLRS_PATH='${molrs_root}' "
            "is not a compatible molrs source tree.")
    endif()

    file(READ "${version_file}" actual_version)
    string(STRIP "${actual_version}" actual_version)
    if(NOT actual_version MATCHES "^[0-9]+$")
        message(FATAL_ERROR
            "Invalid molrs CXX API version in ${version_file}: "
            "expected one unsigned decimal integer, found '${actual_version}'.")
    endif()
    if(NOT actual_version STREQUAL "${ATV_MOLRS_CXX_API_VERSION}")
        message(FATAL_ERROR
            "molrs CXX API mismatch: Atomiverse requires v"
            "${ATV_MOLRS_CXX_API_VERSION}; found v${actual_version} at "
            "ATV_MOLRS_PATH='${molrs_root}'. Local path overrides do not "
            "bypass the contract.")
    endif()
endfunction()
