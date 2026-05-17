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


find_package (BLAS)
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
find_package(FFTW COMPONENTS FLOAT_LIB DOUBLE_LIB LONGDOUBLE_LIB DOUBLE_OPENMP_LIB)
if(FFTW_DOUBLE_LIB_FOUND)
    message(STATUS "FFTW_INCLUDE_DIR is ${FFTW_INCLUDE_DIRS}")
    add_definitions (-DSCTL_HAVE_FFTW)
    include_directories("${FFTW_INCLUDE_DIRS}")
endif()
if (FFTW_FLOAT_LIB_FOUND)
    add_definitions (-DSCTL_HAVE_FFTWF)
endif()
if (FFTW_LONGDOUBLE_LIB_FOUND)
    add_definitions (-DSCTL_HAVE_FFTWL)
endif()
if (FFTW_FLOAT_OPENMP_FOUND)
    add_definitions (-DSCTL_FFTWF_THREADS)
endif()

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
    pybind11_add_module(${PROJECT_NAME} src/python.cpp)
    target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CXX_COMPILER_ID:PGI>:-Wc,--pending_instantiations=1000>")
    target_include_directories(${PROJECT_NAME} PRIVATE "$<$<BOOL:${FFTW_DOUBLE_LIB_FOUND}>:${FFTW_INCLUDE_DIRS}>"
    						       "include" "extern/SCTL/include" "extern/BIEST/include")
    target_link_libraries(${PROJECT_NAME} PRIVATE "$<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>"
        "$<$<BOOL:${FFTW_DOUBLE_LIB_FOUND}>:FFTW::Double>" "$<$<BOOL:${FFTW_FLOAT_LIB_FOUND}>:FFTW::Float>"
        "$<$<BOOL:${FFTW_LONGDOUBLE_LIB_FOUND}>:FFTW::LongDouble>" "$<$<BOOL:${FFTW_FLOAT_OPENMP_FOUND}>:FFTW::DoubleOpenMP>"
        "$<$<BOOL:${BLAS_FOUND}>:${BLAS_LIBRARIES}>" "$<$<BOOL:${LAPACK_FOUND}>:${LAPACK_LIBRARIES}>")

    install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION .)
endif()

add_executable(vc_testing EXCLUDE_FROM_ALL test/virtual-casing.cpp)
# Fix for PowerPC with PGI compiler
target_compile_options(vc_testing PUBLIC "$<$<CXX_COMPILER_ID:PGI>:-Wc,--pending_instantiations=1000>")
target_include_directories(vc_testing PRIVATE "$<$<BOOL:${FFTW_DOUBLE_LIB_FOUND}>:${FFTW_INCLUDE_DIRS}>"
					      "include" "extern/SCTL/include" "extern/BIEST/include")
target_link_libraries(vc_testing PUBLIC $<$<BOOL:${FFTW_DOUBLE_LIB_FOUND}>:FFTW::Double>
    $<$<BOOL:${FFTW_FLOAT_LIB_FOUND}>:FFTW::Float> $<$<BOOL:${FFTW_LONGDOUBLE_LIB_FOUND}>:FFTW::LongDouble>
    $<$<BOOL:${FFTW_FLOAT_OPENMP_FOUND}>:FFTW::DoubleOpenMP> "$<$<BOOL:${BLAS_FOUND}>:${BLAS_LIBRARIES}>"
    "$<$<BOOL:${LAPACK_FOUND}>:${LAPACK_LIBRARIES}>" $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
