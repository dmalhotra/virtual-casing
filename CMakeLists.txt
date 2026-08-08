cmake_minimum_required(VERSION 3.16)
project(virtual_casing)
#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

OPTION (USE_PROFILER "Activate code profiling." OFF)
if (USE_PROFILER)
    add_definitions (-DSCTL_PROFILE=5 -DSCTL_VERBOSE)
endif ()
set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON")

if(DEFINED ENV{CONDA_PREFIX})
    include_directories($ENV{CONDA_PREFIX}/include)
    link_directories("$ENV{CONDA_PREFIX}/lib")
    message(STATUS "Conda prefix is $ENV{CONDA_PREFIX}")
endif()

# Add linker and incldue directories for conda build before finding any library
if (DEFINED ENV{CONDA_BUILD})
  link_directories("$ENV{BUILD_PREFIX}/lib")
  include_directories("$ENV{BUILD_PREFIX}/include")
  message(STATUS "Prefix environment variable is $ENV{BUILD_PREFIX}")
endif ()


# Wheels turn this off (pyproject.toml): one N=1 dgemm, no LAPACK, ~80% of wheel size
option (VC_USE_BLAS "Link BLAS/LAPACK; SCTL falls back to its own routines without them" ON)
if (VC_USE_BLAS)
    find_package (BLAS)
endif ()
if (BLAS_FOUND)
    add_definitions (-DSCTL_HAVE_BLAS)
    if (NOT ${BLA_VENDOR} MATCHES "^Intel")
        find_package (LAPACK)
        if(LAPACK_FOUND)
            add_definitions (-DSCTL_HAVE_LAPACK)
        endif()
    else()
        add_definitions (-DSCTL_HAVE_LAPACK)
    endif()
endif ()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(FFTW COMPONENTS DOUBLE_LIB)
if(FFTW_DOUBLE_LIB_FOUND)
    message(STATUS "FFTW_INCLUDE_DIR is ${FFTW_INCLUDE_DIRS}")
    add_definitions (-DSCTL_HAVE_FFTW)
    include_directories("${FFTW_INCLUDE_DIRS}")
endif()
# float/long-double FFTW is never called; linking it only bloats the wheel

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

include(CheckCXXCompilerFlag)
IF(DEFINED ENV{VC_ARCH_FLAG})
    message(STATUS "User provided architecture flag: $ENV{VC_ARCH_FLAG}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 $ENV{VC_ARCH_FLAG} -ffp-contract=fast")
elseif(DEFINED ENV{CONDA_BUILD})
    message(STATUS "conda build environment detected. Let conda set compilation flags accordingly.")
else()
    # Temporary fix till clang 15 is default on macs with Apple silicon
    message(STATUS "Local build detected. Set compilation flags accordingly.")
    unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
    CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
    if(COMPILER_SUPPORTS_MARCH_NATIVE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffp-contract=fast")
    elseif(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "arm64")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mcpu=native -ffp-contract=fast")
    endif()
endif()

find_package(OpenMP REQUIRED)

if (SKBUILD)
    find_package(pybind11 CONFIG REQUIRED)

    # SCTL selects its vector ISA at compile time and has no runtime dispatch, so
    # reaching AVX2/AVX-512 means compiling the extension once per level and
    # choosing in __init__.py. AVX-512 cannot be a single-wheel target: Intel
    # fused it off on client CPUs from Alder Lake onwards.
    option (VC_DISPATCH "Build one extension per x86-64 instruction-set level" ON)
    set(vc_variants baseline)
    if (VC_DISPATCH AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
        list(APPEND vc_variants v3 v4)
        set(vc_march_v3 "-march=x86-64-v3")
        set(vc_march_v4 "-march=x86-64-v4")
    endif()

    foreach (variant IN LISTS vc_variants)
        set(tgt _vc_${variant})
        pybind11_add_module(${tgt} src/python.cpp)
        target_compile_definitions(${tgt} PRIVATE VC_MODULE_NAME=${tgt})
        # after CMAKE_CXX_FLAGS, so this -march wins over VC_ARCH_FLAG
        target_compile_options(${tgt} PRIVATE ${vc_march_${variant}}
            "$<$<CXX_COMPILER_ID:PGI>:-Wc,--pending_instantiations=1000>")
        target_include_directories(${tgt} PRIVATE "$<$<BOOL:${FFTW_DOUBLE_LIB_FOUND}>:${FFTW_INCLUDE_DIRS}>"
                                                  "include" "extern/SCTL/include" "extern/BIEST/include")
        target_link_libraries(${tgt} PRIVATE "$<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>"
            "$<$<BOOL:${FFTW_DOUBLE_LIB_FOUND}>:FFTW::Double>"
            "$<$<BOOL:${BLAS_FOUND}>:${BLAS_LIBRARIES}>" "$<$<BOOL:${LAPACK_FOUND}>:${LAPACK_LIBRARIES}>")
        install(TARGETS ${tgt} LIBRARY DESTINATION ${PROJECT_NAME})
    endforeach()

    pybind11_add_module(_cpu src/cpu_features.cpp) # baseline: must run anywhere
    install(TARGETS _cpu LIBRARY DESTINATION ${PROJECT_NAME})
    install(FILES src/__init__.py DESTINATION ${PROJECT_NAME})
endif()

add_executable(vc_testing EXCLUDE_FROM_ALL test/virtual-casing.cpp)
# Fix for PowerPC with PGI compiler
target_compile_options(vc_testing PUBLIC "$<$<CXX_COMPILER_ID:PGI>:-Wc,--pending_instantiations=1000>")
target_include_directories(vc_testing PRIVATE "$<$<BOOL:${FFTW_DOUBLE_LIB_FOUND}>:${FFTW_INCLUDE_DIRS}>"
					      "include" "extern/SCTL/include" "extern/BIEST/include")
target_link_libraries(vc_testing PUBLIC $<$<BOOL:${FFTW_DOUBLE_LIB_FOUND}>:FFTW::Double>
    "$<$<BOOL:${BLAS_FOUND}>:${BLAS_LIBRARIES}>"
    "$<$<BOOL:${LAPACK_FOUND}>:${LAPACK_LIBRARIES}>" $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
