# The spot finder and the image reader, included from the top-level
# CMakeLists, which requires HDF5 and the bitshuffle submodule.
#
# One project, one C++ standard, one build type. Kept here: the version check against
# vcpkg.json, because vcpkg reads that file before CMake runs and the two would
# otherwise drift -- "which build is this?" is the first question asked in an
# incident, and a quietly wrong answer is worse than none.

enable_language(C)  # bitshuffle and lz4

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg.json" mxi_manifest)
string(REGEX MATCH "\"version-string\"[ \t]*:[ \t]*\"([^\"]+)\""
       mxi_manifest_version "${mxi_manifest}")
if(NOT CMAKE_MATCH_1 STREQUAL PROJECT_VERSION)
    message(FATAL_ERROR
            "version-string in vcpkg.json is ${CMAKE_MATCH_1} but project() says "
            "${PROJECT_VERSION}; make them agree")
endif()


# Fast math on the device, for measuring what it costs and what it buys. Off,
# because the whole claim the device backends make is that they agree with
# dext.cc bit for bit, and fast math is licence to reassociate the very
# arithmetic that claim is about -- the Metal compiler has it on unless told
# otherwise, which is why -fno-fast-math is passed explicitly below.
#
# It applies to the device only. dext.cc is the reference the device is compared
# against, so compiling that with fast math would move the thing being measured
# and leave nothing to measure it with.
#
# With it on, ctest -R dext_gpu tolerates a background differing in its last
# bits and says by how much, but still fails if the set of pixels found has
# changed; bench_dext_gpu reports the difference and prints its timings anyway,
# since the timings are the reason for the build. Named without the project's
# usual prefix because that is what was asked for.
# The spot finder's options were SPOTFINDER_*, from when it was a project of
# its own; mxi is one project, and they are MXI_* (Graeme, 1.0). The old names
# are taken, with a warning, through 1.0.x.
foreach(name CUDA METAL AVX2 TESTS SANITIZE)
    if(DEFINED SPOTFINDER_${name})
        message(DEPRECATION
                "SPOTFINDER_${name} is MXI_${name} now; taking it as that")
        set(MXI_${name} "${SPOTFINDER_${name}}" CACHE STRING "" FORCE)
        unset(SPOTFINDER_${name} CACHE)  # said once, not at every configure
    endif()
endforeach()

option(ENABLE_FAST_MATH "Compile the device kernels with fast math" OFF)

option(MXI_TESTS "Build the tests" ON)
option(MXI_CUDA "Build the CUDA threshold kernels" OFF)
option(MXI_METAL "Build the Metal threshold kernels" OFF)

# One backend at a time. gpu::find has one definition, and a build with both
# would either fail at the linker or -- worse -- link whichever the linker saw
# first, so that --gpu ran a device nobody chose.
if(MXI_CUDA AND MXI_METAL)
    message(FATAL_ERROR
            "MXI_CUDA and MXI_METAL are both on, and they implement "
            "the same interface. Pick one.")
endif()

# CUDA is enabled as a language only when asked for, so that a machine without a
# toolkit still configures. check_language reports its absence as a plain message
# rather than the hard failure enable_language would give.
if(MXI_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(NOT CMAKE_CUDA_COMPILER)
        message(FATAL_ERROR
                "MXI_CUDA is on but no CUDA compiler was found. Put nvcc "
                "on PATH, or set CMAKE_CUDA_COMPILER, or configure with "
                "-DMXI_CUDA=OFF.")
    endif()
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    # Left to the user, since the right answer is whichever cards this will run
    # on. "native" needs CMake 3.24; before that, name them, e.g. "70;80;90".
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
            set(CMAKE_CUDA_ARCHITECTURES native)
        else()
            set(CMAKE_CUDA_ARCHITECTURES 70)
        endif()
    endif()
    message(STATUS "CUDA: building for architectures ${CMAKE_CUDA_ARCHITECTURES}")
endif()

# Metal, likewise, is off by default and everything it needs is looked for only
# when it is asked for, so that a Linux machine still configures.
#
# metal-cpp is a header-only wrapper Apple distributes as a zip rather than
# through any package manager, so it is pointed at rather than found:
#
#   cmake -S . -B build -DMXI_METAL=ON -DMETAL_CPP_DIR=/path/to/metal-cpp
if(MXI_METAL)
    if(NOT APPLE)
        message(FATAL_ERROR
                "MXI_METAL is on but this is not macOS. The CPU path runs "
                "everywhere; -DMXI_CUDA=ON is the other device.")
    endif()

    find_path(METAL_CPP_INCLUDE_DIR
              NAMES Metal/Metal.hpp
              HINTS ${METAL_CPP_DIR} ENV METAL_CPP_DIR
              PATH_SUFFIXES metal-cpp
              DOC "Directory containing Metal/Metal.hpp from Apple's metal-cpp")
    if(NOT METAL_CPP_INCLUDE_DIR)
        message(FATAL_ERROR
                "metal-cpp was not found. Download it from "
                "https://developer.apple.com/metal/cpp/ and configure with "
                "-DMETAL_CPP_DIR=<the unpacked directory>.")
    endif()

    find_program(MXI_XCRUN xcrun)
    if(NOT MXI_XCRUN)
        message(FATAL_ERROR
                "xcrun was not found, so the shader cannot be compiled. Install "
                "the Xcode command line tools.")
    endif()

    find_library(METAL_FRAMEWORK Metal REQUIRED)
    find_library(FOUNDATION_FRAMEWORK Foundation REQUIRED)

    message(STATUS "Metal: metal-cpp at ${METAL_CPP_INCLUDE_DIR}")

    # The shader is compiled at build time rather than from source at run time.
    # newLibraryWithSource needs the Metal compiler present on the machine that
    # runs the binary, which is a beamline machine, and a shader that only fails
    # to compile at the first frame of an acquisition is a bad trade for the
    # second or so it saves here.
    #
    # -fno-fast-math because the point of that file is that it agrees with
    # dext.cc: fast math is on by default and would reassociate the dispersion
    # test out from under the comparison. -DENABLE_FAST_MATH asks for the other
    # one, deliberately.
    if(ENABLE_FAST_MATH)
        set(mxi_metal_math -ffast-math)
    else()
        set(mxi_metal_math -fno-fast-math)
    endif()

    set(mxi_metal_air "${CMAKE_CURRENT_BINARY_DIR}/dext_metal.air")
    set(mxi_metallib "${CMAKE_CURRENT_BINARY_DIR}/dext_metal.metallib")
    set(mxi_metallib_cc "${CMAKE_CURRENT_BINARY_DIR}/dext_metal_library.cc")

    add_custom_command(
        OUTPUT ${mxi_metal_air}
        COMMAND ${MXI_XCRUN} -sdk macosx metal
                -std=metal3.0 ${mxi_metal_math} -Wall -Werror
                -c "${CMAKE_CURRENT_SOURCE_DIR}/src/dext_metal.metal"
                -o ${mxi_metal_air}
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/dext_metal.metal"
        COMMENT "Compiling src/dext_metal.metal")

    add_custom_command(
        OUTPUT ${mxi_metallib}
        COMMAND ${MXI_XCRUN} -sdk macosx metallib ${mxi_metal_air}
                -o ${mxi_metallib}
        DEPENDS ${mxi_metal_air}
        COMMENT "Linking dext_metal.metallib")

    add_custom_command(
        OUTPUT ${mxi_metallib_cc}
        COMMAND ${CMAKE_COMMAND} -DINPUT=${mxi_metallib}
                -DOUTPUT=${mxi_metallib_cc}
                -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_metallib.cmake"
        DEPENDS ${mxi_metallib}
                "${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_metallib.cmake"
        COMMENT "Embedding dext_metal.metallib")
endif()

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/bitshuffle/src/bitshuffle.h")
    message(FATAL_ERROR
            "third_party/bitshuffle is missing: run 'git submodule update --init'")
endif()

# Sanitizers, for chasing a fault that only appears under load. Applied with
# add_compile_options rather than CMAKE_CXX_FLAGS so that they reach everything
# built in this tree -- bitshuffle included, which matters: a sanitizer only
# reports an access it instrumented, so a read off the end of a buffer inside
# bitshuffle is invisible unless bitshuffle itself was built this way.
set(MXI_SANITIZE "" CACHE STRING
    "Sanitizers to build with: ON (address and undefined), address, undefined, address,undefined, or thread")
set_property(CACHE MXI_SANITIZE PROPERTY STRINGS
             "" ON address undefined "address,undefined" thread)
# ON, as the one switch this build had before the spot finder's list joined it:
# address and undefined behaviour.
if(MXI_SANITIZE STREQUAL "ON" OR MXI_SANITIZE STREQUAL "TRUE")
    set(MXI_SANITIZE "address,undefined")
endif()

if(MXI_SANITIZE)
    if(MSVC)
        message(FATAL_ERROR "MXI_SANITIZE is not supported with MSVC here")
    endif()
    # nvcc does not take -fsanitize itself; it has to be handed to the host
    # compiler, one flag per sanitizer, because nvcc splits -Xcompiler arguments
    # on commas. A device fault needs compute-sanitizer instead.
    string(REPLACE "," ";" mxi_sanitizers "${MXI_SANITIZE}")
    set(mxi_sanitize_host "")
    set(mxi_sanitize_cuda "")
    foreach(sanitizer IN LISTS mxi_sanitizers)
        list(APPEND mxi_sanitize_host "-fsanitize=${sanitizer}")
        list(APPEND mxi_sanitize_cuda "-Xcompiler=-fsanitize=${sanitizer}")
    endforeach()
    # Undefined behaviour ABORTS rather than prints and carries on: otherwise a
    # sanitizer run reports and the test still passes, which cannot fail.
    if("undefined" IN_LIST mxi_sanitizers)
        list(APPEND mxi_sanitize_host "-fno-sanitize-recover=undefined")
        list(APPEND mxi_sanitize_cuda "-Xcompiler=-fno-sanitize-recover=undefined")
    endif()

    add_compile_options(
        "$<$<COMPILE_LANGUAGE:C,CXX>:${mxi_sanitize_host};-fno-omit-frame-pointer>"
        "$<$<COMPILE_LANGUAGE:CUDA>:${mxi_sanitize_cuda};-Xcompiler=-fno-omit-frame-pointer>")
    add_link_options(${mxi_sanitize_host})
    message(STATUS "sanitizers: ${MXI_SANITIZE}, applied to everything built here")
endif()

# Warning flags as usage requirements, so that every target picks them up by
# linking this and the compiler check lives in one place. An INTERFACE library
# with no location contributes nothing to the link line.
add_library(mxi_warnings INTERFACE)
if(NOT MSVC)
    target_compile_options(mxi_warnings INTERFACE -Wall -Wextra)
endif()

# bitshuffle has no CMake build of its own -- it is packaged for Python -- so the
# C files that make up the library are compiled here. bshuf_h5filter.c with
# them: the .rflx writer compresses its integer datasets with bitshuffle-LZ4,
# HDF5 filter 32008, through HDF5's filter pipeline (docs/rflx.md). The plugin
# sources, for loading the filter into another program's HDF5, are not needed.
add_library(bitshuffle STATIC
    third_party/bitshuffle/src/bitshuffle.c
    third_party/bitshuffle/src/bitshuffle_core.c
    third_party/bitshuffle/src/iochain.c
    third_party/bitshuffle/src/bshuf_h5filter.c
    third_party/bitshuffle/lz4/lz4.c)
set_target_properties(bitshuffle PROPERTIES C_STANDARD 99)
target_include_directories(bitshuffle SYSTEM PUBLIC
    third_party/bitshuffle/src third_party/bitshuffle/lz4)
target_include_directories(bitshuffle SYSTEM PRIVATE ${HDF5_C_INCLUDE_DIRS})

# bitshuffle chooses its SIMD path at compile time, from what the compiler says
# the target has: __SSE2__, __AVX2__, or __ARM_NEON with __aarch64__. Two of those
# need no help, since SSE2 is baseline on x86-64 and NEON is baseline on aarch64
# -- a plain build is already vectorised on both. AVX2 is the tier that has to be
# asked for, and asking for it is only meaningful on x86.
option(MXI_AVX2
       "Compile bitshuffle's AVX2 path; the binary then requires an AVX2 CPU" ON)

set(mxi_target_arch "${CMAKE_SYSTEM_PROCESSOR}")
if(CMAKE_OSX_ARCHITECTURES)
    # A universal or cross build overrides the host processor and may name more
    # than one architecture, so "arm64" and "arm64;x86_64" both fall through to
    # the else branch below -- which is what we want, since -mavx2 would fail the
    # arm64 slice.
    set(mxi_target_arch "${CMAKE_OSX_ARCHITECTURES}")
endif()

if(NOT MXI_AVX2)
    message(STATUS "bitshuffle: AVX2 off by request")
elseif(NOT mxi_target_arch MATCHES "^(x86_64|amd64|AMD64|i[3-6]86)$")
    message(STATUS
            "bitshuffle: ${mxi_target_arch} is not x86, so AVX2 does not "
            "apply; NEON is baseline on aarch64")
else()
    include(CheckCCompilerFlag)
    if(MSVC)
        set(mxi_avx2_flag /arch:AVX2)
    else()
        set(mxi_avx2_flag -mavx2)
    endif()
    # Ask the compiler rather than assuming: a clang targeting arm64 rejects
    # -mavx2 outright, and this is the last line of defence if the architecture
    # test above is ever fooled.
    check_c_compiler_flag(${mxi_avx2_flag} MXI_COMPILER_TAKES_AVX2)
    if(MXI_COMPILER_TAKES_AVX2)
        target_compile_options(bitshuffle PRIVATE ${mxi_avx2_flag})
        message(STATUS
                "bitshuffle: AVX2 path enabled, so this build will not run on a "
                "CPU without AVX2 -- configure with -DMXI_AVX2=OFF for a "
                "portable SSE2 build")
    else()
        message(STATUS
                "bitshuffle: ${mxi_avx2_flag} rejected by the C compiler, "
                "falling back to SSE2")
    endif()
endif()

# What every target carrying device kernels is given. Collected here so that the
# executable, the device test and the benchmark cannot end up compiled three
# different ways -- which would make the test's verdict about a different binary
# from the one being timed.
set(mxi_gpu_definitions "")
set(mxi_cuda_options "")
if(ENABLE_FAST_MATH)
    if(MXI_CUDA OR MXI_METAL)
        list(APPEND mxi_gpu_definitions MXI_FAST_MATH)
        list(APPEND mxi_cuda_options
             "$<$<COMPILE_LANGUAGE:CUDA>:-use_fast_math>")
        message(STATUS
                "fast math: ON for the device kernels. They no longer agree "
                "with dext.cc bit for bit, so run ctest -R dext_gpu to see "
                "what changed before believing a .refl from this build.")
    else()
        # Rather than silently doing nothing: fast math here means the device
        # kernels, and dext.cc is left alone on purpose because it is the
        # reference they are measured against.
        message(WARNING
                "ENABLE_FAST_MATH is on but neither device backend is being "
                "built. It applies to the device kernels only -- the CPU "
                "threshold is the reference and is never compiled with it -- so "
                "this build is unaffected.")
    endif()
endif()

find_package(Threads REQUIRED)

# HDF5, the C library only, and only H5D chunk reads plus the virtual-dataset
# property list calls. Both need 1.10.3 or newer.
find_package(HDF5 REQUIRED COMPONENTS C)

# HDF5 2.0 gave H5Dread_chunk a buffer-size argument, and its versioned-API macros
# point the old name at the new function, so the call has two shapes and no macro
# tells them apart. Ask the compiler which one is in front of it rather than
# guessing from H5_VERSION_GE: vcpkg and the distributions are on different sides
# of this. Compiled as a static library so the probe does not need to link.
include(CheckCXXSourceCompiles)
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
set(CMAKE_REQUIRED_INCLUDES ${HDF5_C_INCLUDE_DIRS})
check_cxx_source_compiles("
#include <hdf5.h>
#include <cstddef>
#include <cstdint>
herr_t probe(hid_t dataset, const hsize_t *offset, std::uint32_t *filters,
             void *buffer, std::size_t *size) {
  return H5Dread_chunk(dataset, H5P_DEFAULT, offset, filters, buffer, size);
}
" MXI_H5DREAD_CHUNK_TAKES_SIZE)
unset(CMAKE_REQUIRED_INCLUDES)
unset(CMAKE_TRY_COMPILE_TARGET_TYPE)
if(MXI_H5DREAD_CHUNK_TAKES_SIZE)
    message(STATUS "HDF5: H5Dread_chunk takes a buffer size (2.0 and later)")
else()
    message(STATUS "HDF5: H5Dread_chunk takes no buffer size (before 2.0)")
endif()

# Where frames come from: an NXmx master file whose virtual dataset is unpacked
# so its chunks can be read directly. HDF5 stays behind this library rather than
# reaching the tool.
add_library(mxi_image_series STATIC src/nxmx_reader.cc)
target_include_directories(mxi_image_series PUBLIC src)
target_include_directories(mxi_image_series SYSTEM PRIVATE ${HDF5_C_INCLUDE_DIRS})
target_link_libraries(mxi_image_series
    PUBLIC mxi PRIVATE ${HDF5_C_LIBRARIES} mxi_warnings)
target_compile_definitions(mxi_image_series PRIVATE
    $<$<BOOL:${MXI_H5DREAD_CHUNK_TAKES_SIZE}>:MXI_H5DREAD_CHUNK_TAKES_SIZE>)

# bitshuffle+LZ4, plain LZ4 and uncompressed chunks.
add_library(mxi_decompress STATIC src/decompress.cc)
target_include_directories(mxi_decompress PUBLIC src)
target_link_libraries(mxi_decompress PRIVATE bitshuffle mxi_warnings)

# The extended dispersion threshold, explicitly instantiated for 16 and 32 bit
# pixels.
add_library(mxi_dext STATIC src/dext.cc)
target_include_directories(mxi_dext PUBLIC src)
target_link_libraries(mxi_dext PRIVATE mxi_warnings)

# Ordering a device's signal pixels by index. Host code, shared by both device
# backends, and its own library because it is worth testing on a machine with no
# GPU at all -- which is where it was written and where its bug was found.
add_library(mxi_signal_order STATIC src/signal_order.cc)
target_include_directories(mxi_signal_order PUBLIC src)
target_link_libraries(mxi_signal_order PRIVATE mxi_warnings)

# The backend-agnostic half of gpu::: which window each stage uses, whether to
# profile, and the last frame's split. Plain C++ with no toolkit in it, so it
# builds whether the backend is CUDA or Metal and is compiled once either way.
add_library(mxi_dext_gpu STATIC src/dext_gpu.cc)
target_include_directories(mxi_dext_gpu PUBLIC src)
target_link_libraries(mxi_dext_gpu PRIVATE mxi_warnings)

# Grouping signal pixels six-connected in three dimensions, as DIALS does, with
# its centroids and its filters. No Boost: the streaming flush renumbers
# components as it compacts, and a union-find it owns is both smaller than
# adjacency_list and the reason the whole sweep need not be held in memory.
add_library(mxi_spots STATIC src/spots.cc)
target_include_directories(mxi_spots PUBLIC src)
target_link_libraries(mxi_spots PRIVATE mxi_warnings)

# Writing a DIALS reflection table. A .refl is msgpack around raw column dumps,
# so this needs neither DIALS nor a msgpack library; what it does need is the
# column type names to be exactly right, since DIALS refuses a name it does not
# know.


add_executable(mxi_find apps/mxi_find.cc)
target_compile_definitions(mxi_find PRIVATE
                           MXI_VERSION="${PROJECT_VERSION}")
target_link_libraries(mxi_find
    PRIVATE mxi mxi_image_series mxi_decompress mxi_dext
            mxi_spots
            Threads::Threads mxi_warnings)

# The kernels go into the executable rather than a library of their own, so that
# adding more .cu files is a one-line change here.
if(MXI_CUDA)
    target_sources(mxi_find PRIVATE src/dext_cuda.cu)
    target_compile_definitions(mxi_find PRIVATE
                               MXI_CUDA MXI_GPU
                               ${mxi_gpu_definitions})
    target_link_libraries(mxi_find
        PRIVATE mxi_signal_order mxi_dext_gpu)
    target_compile_options(mxi_find PRIVATE ${mxi_cuda_options})
    # No separable compilation: there is one .cu file and nothing calls across
    # translation units on the device, so the extra nvcc device-link step buys
    # nothing -- and it would be handed the host linker's flags, which it does
    # not understand.
    set_target_properties(mxi_find
                          PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)
endif()

# Metal, unlike CUDA, goes into a library of its own rather than sources hung on
# the executable, because tests/spots/test_dext_gpu.cc has to link it too: a backend
# only reachable through the tool cannot be diffed against the CPU without a
# detector and a series.
if(MXI_METAL)
    add_library(mxi_dext_metal STATIC
                src/dext_metal.cc ${mxi_metallib_cc})
    target_include_directories(mxi_dext_metal PUBLIC src)
    # SYSTEM: metal-cpp is not clean under -Wall -Wextra and it is not ours to
    # fix. The warnings that matter are mxi's own.
    target_include_directories(mxi_dext_metal SYSTEM PUBLIC
                               ${METAL_CPP_INCLUDE_DIR})
    target_link_libraries(mxi_dext_metal
        PUBLIC ${METAL_FRAMEWORK} ${FOUNDATION_FRAMEWORK}
        PRIVATE mxi_signal_order mxi_dext_gpu mxi_warnings)

    target_link_libraries(mxi_find PRIVATE mxi_dext_metal)
    target_compile_definitions(mxi_find PRIVATE
                               MXI_METAL MXI_GPU
                               ${mxi_gpu_definitions})
endif()

install(TARGETS mxi_find RUNTIME DESTINATION bin)

if(MXI_TESTS)
    enable_testing()

    # Not a regression: it pins the domain the 32-bit window sums are exact over,
    # and fails if a change to the kernel size or to the masking eats the
    # headroom.
    add_executable(test_dext_squares tests/spots/test_dext_squares.cc)
    target_link_libraries(test_dext_squares
        PRIVATE mxi_dext mxi_warnings)
    add_test(NAME dext_squares COMMAND test_dext_squares)

    # No GPU needed: the fused kernel's phases, run on the CPU tile by tile,
    # against dext() -- the tile arithmetic, checked where it can be.
    add_executable(test_dext_fused tests/spots/test_dext_fused.cc)
    target_link_libraries(test_dext_fused
        PRIVATE mxi_dext mxi_warnings)
    add_test(NAME dext_fused COMMAND test_dext_fused)

    # No GPU needed: the ordering is host code, and this is where its heap
    # overflow was caught.
    add_executable(test_signal_order tests/spots/test_signal_order.cc)
    target_link_libraries(test_signal_order
        PRIVATE mxi_signal_order mxi_warnings)
    add_test(NAME signal_order COMMAND test_signal_order)

    # The DIALS side: the grouping, the file, the experiment list, and the three
    # of them together over planted frames. None of these needs a GPU, a server
    # or a detector, which is deliberate -- they are the tests that say whether
    # what dials.index is handed means what it says.
    add_executable(test_dials_spots tests/spots/test_dials_spots.cc)
    target_link_libraries(test_dials_spots
        PRIVATE mxi_spots mxi_warnings)
    add_test(NAME dials_spots COMMAND test_dials_spots)

    # Not a ctest: it is a measurement to read, not a pass or a fail. It is
    # here because the grouping looks like the expensive stage -- it is the one
    # that cannot be spread over threads -- and at the density a real frame has
    # it is a tenth of a per cent of one. The numbers are in its header.
    add_executable(bench_dials_spots tests/spots/bench_dials_spots.cc)
    target_link_libraries(bench_dials_spots
        PRIVATE mxi_spots mxi_warnings)

    # mxi_integrate needs pixels, so it lives here with the HDF5 guard rather
    # than in the top-level list: everything else in apps/ builds without HDF5.
    # Every program mirrors its output to mxi_<program>.log where it runs;
    # checked through real binaries, one of them the spot finder, whose main is
    # its own, so it lives inside this guard with it.
    add_test(NAME log_mirror
             COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/tests/log_mirror.sh
                     $<TARGET_FILE:mxi_background> $<TARGET_FILE:mxi_index>
                     $<TARGET_FILE:mxi_find>)

    add_executable(mxi_integrate ${CMAKE_CURRENT_SOURCE_DIR}/apps/mxi_integrate.cc)
    target_link_libraries(mxi_integrate
                          PRIVATE mxi mxi_image_series mxi_decompress
                                  mxi_warnings)
    target_include_directories(mxi_integrate PRIVATE src apps)

    # An experiment list from an NXmx master, as dials.import writes it.
    add_executable(mxi_import ${CMAKE_CURRENT_SOURCE_DIR}/apps/mxi_import.cc
                              ${CMAKE_CURRENT_SOURCE_DIR}/src/nxmx_import.cc)
    target_link_libraries(mxi_import PRIVATE mxi mxi_warnings ${HDF5_C_LIBRARIES})
    target_include_directories(mxi_import PRIVATE src apps)
    target_include_directories(mxi_import SYSTEM PRIVATE ${HDF5_C_INCLUDE_DIRS})

    # The largest count on a valid pixel, and whether it fits in 16 bits.
    add_executable(mxi_max ${CMAKE_CURRENT_SOURCE_DIR}/apps/mxi_max.cc)
    target_link_libraries(mxi_max
                          PRIVATE mxi mxi_image_series mxi_decompress
                                  mxi_warnings)
    target_include_directories(mxi_max PRIVATE src apps)

    # How fast the images can be read and decompressed, apart from the
    # programs: H5Dread_chunk against pread. Not used by the pipeline.
    add_executable(mxi_readtest ${CMAKE_CURRENT_SOURCE_DIR}/apps/mxi_readtest.cc)
    target_link_libraries(mxi_readtest
                          PRIVATE mxi mxi_image_series mxi_decompress
                                  mxi_warnings)
    target_include_directories(mxi_readtest PRIVATE src apps)
    target_include_directories(mxi_readtest SYSTEM PRIVATE ${HDF5_C_INCLUDE_DIRS})
    target_link_libraries(mxi_readtest PRIVATE ${HDF5_C_LIBRARIES})

    add_executable(test_refl tests/spots/test_refl.cc)
    target_link_libraries(test_refl PRIVATE mxi mxi_spots mxi_warnings)
    add_test(NAME refl COMMAND test_refl)

    # Holds the writer's output to what it produces today, so that replacing it
    # with the general writer in ../src/refl.cc can be shown to change nothing.
    # See docs/spotfinder.md: this must not be weakened to pass.
    add_executable(test_refl_golden tests/spots/test_refl_golden.cc)
    # mxi links the table writer; mxi_spots the spots themselves
    # too puts the archive on the link line twice and ld warns about it.
    target_link_libraries(test_refl_golden
                          PRIVATE mxi mxi_spots mxi_warnings)
    add_test(NAME refl_golden COMMAND test_refl_golden)

    add_executable(test_expt tests/spots/test_expt.cc)
    target_link_libraries(test_expt PRIVATE mxi mxi_warnings)
    add_test(NAME expt COMMAND test_expt)

    add_executable(test_dials_end_to_end tests/spots/test_dials_end_to_end.cc)
    target_link_libraries(test_dials_end_to_end
        PRIVATE mxi_dext mxi mxi_spots mxi_warnings)
    add_test(NAME dials_end_to_end COMMAND test_dials_end_to_end)

    # The CPU and the device must agree exactly, so this is only built when there
    # is a device to disagree with. It skips, rather than fails, when the backend
    # is compiled in but the machine has none. Either backend: the CUDA kernels
    # are compiled into the executable rather than a library, so the test needs
    # the .cu too.
    if(MXI_CUDA OR MXI_METAL)
        set(mxi_gpu_test_libs mxi_dext mxi_signal_order
                                     mxi_dext_gpu mxi_warnings)
        if(MXI_METAL)
            list(APPEND mxi_gpu_test_libs mxi_dext_metal)
            set(mxi_gpu_test_sources "")
        else()
            set(mxi_gpu_test_sources src/dext_cuda.cu)
        endif()

        add_executable(test_dext_gpu tests/spots/test_dext_gpu.cc
                       ${mxi_gpu_test_sources})
        target_include_directories(test_dext_gpu PRIVATE tests/spots)
        target_link_libraries(test_dext_gpu PRIVATE ${mxi_gpu_test_libs})
        target_compile_definitions(test_dext_gpu PRIVATE
                                   ${mxi_gpu_definitions})
        target_compile_options(test_dext_gpu PRIVATE ${mxi_cuda_options})
        add_test(NAME dext_gpu COMMAND test_dext_gpu)
        set_tests_properties(dext_gpu PROPERTIES SKIP_RETURN_CODE 77)

        # Which window is faster is a question about the hardware, so it is
        # answered by measurement. Not a ctest: it takes minutes and its output
        # is a number to read, not a pass or a fail.
        add_executable(bench_dext_gpu tests/spots/bench_dext_gpu.cc
                       ${mxi_gpu_test_sources})
        target_include_directories(bench_dext_gpu PRIVATE tests/spots)
        target_link_libraries(bench_dext_gpu
            PRIVATE ${mxi_gpu_test_libs} Threads::Threads)
        target_compile_definitions(bench_dext_gpu PRIVATE
                                   ${mxi_gpu_definitions})
        target_compile_options(bench_dext_gpu PRIVATE ${mxi_cuda_options})
    endif()
endif()
