cmake_minimum_required(VERSION 3.20)
# The version is here and in vcpkg.json, and cmake/spots.cmake checks the two
# against each other: vcpkg reads its manifest before CMake runs, so nothing
# else would catch them drifting. "Which build is this?" is the first question
# asked in an incident and a quietly wrong answer is worse than none.
project(mxi VERSION 1.1.2 LANGUAGES CXX)

# Installed programs keep the paths to the libraries they were linked against
# outside the system's own -- HDF5 from conda or Homebrew -- or they would not
# find them; and the install has nothing to rewrite.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# No third-party dependencies, deliberately. The spot finder hand-rolled its
# msgpack writer for the same reason: this has to build on a laptop and on a
# beamline machine with nothing installed on it.
# FFTW for the transform, which is the largest phase of indexing on a fast
# machine. The built-in radix-2 in src/fft.cc is the default and needs nothing;
# this is faster and is GPL, which is why it is opt-in rather than assumed.
#
#   sudo apt-get install libfftw3-dev      # Debian and Ubuntu
#   brew install fftw                      # macOS
#   cmake -S . -B build -DMXI_FFTW=ON
option(MXI_FFTW "Use FFTW for the transform instead of the built-in radix-2" OFF)
if(MXI_FFTW)
    # pkg-config first, through an imported target rather than by hand.
    # pkg_check_modules sets FFTW3_LIBRARIES to the link NAME -- "fftw3" -- and
    # keeps the directory in FFTW3_LIBRARY_DIRS, so linking the first without
    # the second works only where the library is already on the linker's
    # default path. It is on Linux, and it is not under Homebrew's
    # /opt/homebrew, which is where this failed. The imported target carries
    # both.
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(FFTW3 QUIET IMPORTED_TARGET fftw3)
    endif()
    if(TARGET PkgConfig::FFTW3)
        set(MXI_FFTW_TARGET PkgConfig::FFTW3)
    else()
        # No pkg-config, or no .pc file. find_library returns a full path, so
        # this needs no directory of its own -- but it does need to be told
        # where Homebrew puts things, since /opt/homebrew is not searched by
        # default.
        find_library(FFTW3_LIBRARY NAMES fftw3
                     HINTS ENV FFTW_ROOT ENV HOMEBREW_PREFIX
                     PATHS /opt/homebrew /usr/local /opt/local
                     PATH_SUFFIXES lib)
        find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h
                  HINTS ENV FFTW_ROOT ENV HOMEBREW_PREFIX
                  PATHS /opt/homebrew /usr/local /opt/local
                  PATH_SUFFIXES include)
        if(FFTW3_LIBRARY AND FFTW3_INCLUDE_DIR)
            add_library(mxi_fftw INTERFACE)
            target_link_libraries(mxi_fftw INTERFACE ${FFTW3_LIBRARY})
            target_include_directories(mxi_fftw INTERFACE ${FFTW3_INCLUDE_DIR})
            set(MXI_FFTW_TARGET mxi_fftw)
        endif()
    endif()
    if(NOT MXI_FFTW_TARGET)
        message(FATAL_ERROR
                "MXI_FFTW is on but FFTW was not found. Install it with "
                "'apt-get install libfftw3-dev' or 'brew install fftw', or "
                "point at it with -DFFTW3_LIBRARY=<path to libfftw3> and "
                "-DFFTW3_INCLUDE_DIR=<path to fftw3.h>. Configure without "
                "-DMXI_FFTW=ON to use the built-in transform.")
    endif()
    # FFTW's own threading, if it was built with any. Two flavours exist and
    # either will do: fftw3_threads is pthreads, fftw3_omp is OpenMP. Debian
    # ships the second, Homebrew the first. Absent both, the transform runs on
    # one thread and everything still works.
    find_library(FFTW3_THREADS_LIBRARY NAMES fftw3_threads fftw3_omp
                 HINTS ENV FFTW_ROOT ENV HOMEBREW_PREFIX
                 PATHS /opt/homebrew /usr/local /opt/local
                 PATH_SUFFIXES lib)
    if(FFTW3_THREADS_LIBRARY)
        message(STATUS "FFT: FFTW (${MXI_FFTW_TARGET}), threaded via "
                       "${FFTW3_THREADS_LIBRARY}")
    else()
        message(STATUS "FFT: FFTW (${MXI_FFTW_TARGET}), single threaded -- "
                       "no fftw3_threads or fftw3_omp found")
    endif()
else()
    message(STATUS "FFT: the built-in radix-2; -DMXI_FFTW=ON for FFTW")
endif()

# gemmi, for space groups: operators, the reciprocal asymmetric unit, absences.
# A submodule pinned at a release, as bitshuffle is. Only its symmetry source
# is compiled -- the rest used so far is headers -- and with its own warnings
# off: the no-warnings rule is about this code, and gemmi's headers come in as
# system headers so that it stays so.
if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/gemmi/include/gemmi/symmetry.hpp)
    message(FATAL_ERROR "third_party/gemmi is empty: run git submodule update --init")
endif()
add_library(mxi_gemmi STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/gemmi/src/symmetry.cpp)
target_include_directories(mxi_gemmi SYSTEM PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/gemmi/include)
target_compile_options(mxi_gemmi PRIVATE -w)
set_target_properties(mxi_gemmi PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(mxi
  src/linalg.cc
  src/geometry.cc
  src/predict.cc
  src/json.cc
  src/expt.cc
  src/refl.cc
  src/rflx.cc
  src/cosym.cc
  src/views.cc
  src/processing.cc
  src/two_theta.cc
  src/correlation.cc
  src/spot_table.cc
  src/untrusted.cc
  src/experiment_images.cc
  src/hdf5_handle.cc
  src/combine.cc
  src/fft.cc
  $<$<BOOL:${MXI_FFTW}>:src/fft_fftw.cc>
  src/index.cc
  src/refine.cc
  src/args.cc
  src/derivatives.cc
  src/shoebox.cc
  src/profile_model.cc
  src/profile_grid.cc
  src/mask.cc
  src/forward.cc
  src/background.cc
  src/integrate.cc
  src/reference.cc
  src/summary.cc
  src/postrefine.cc
  src/symmetry.cc
  src/laue.cc
  src/scale_model.cc
  src/parallel.cc
  src/fit_batch.cc
  src/scale.cc
  src/resolution.cc
)

target_link_libraries(mxi PUBLIC mxi_gemmi)

if(MXI_FFTW)
    # PUBLIC: fft3d is called from mxi and the symbol has to resolve in every
    # executable that links mxi, not only inside the archive.
    target_compile_definitions(mxi PUBLIC MXI_USE_FFTW)
    target_link_libraries(mxi PUBLIC ${MXI_FFTW_TARGET})
    if(FFTW3_THREADS_LIBRARY)
        # Before the main library on the link line, as FFTW asks.
        target_link_libraries(mxi PUBLIC ${FFTW3_THREADS_LIBRARY})
        target_compile_definitions(mxi PRIVATE MXI_FFTW_THREADS)
        find_package(Threads REQUIRED)
        target_link_libraries(mxi PUBLIC Threads::Threads)
    endif()
endif()
# Always: src/parallel.cc runs std::threads, FFTW's threads or not.
find_package(Threads REQUIRED)
target_link_libraries(mxi PUBLIC Threads::Threads)
target_include_directories(mxi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_options(mxi PRIVATE -Wall -Wextra -Wpedantic)

add_library(mxi_index OBJECT apps/mxi_index.cc)
target_link_libraries(mxi_index PRIVATE mxi)
target_compile_options(mxi_index PRIVATE -Wall -Wextra)

add_library(mxi_refine OBJECT apps/mxi_refine.cc)
target_link_libraries(mxi_refine PRIVATE mxi)
target_compile_options(mxi_refine PRIVATE -Wall -Wextra)

add_library(mxi_scale OBJECT apps/mxi_scale.cc)
target_link_libraries(mxi_scale PRIVATE mxi)
target_compile_options(mxi_scale PRIVATE -Wall -Wextra)

# gemmi's MTZ writer, for mxi_export: a library of its own, as it wants zlib
# and nothing else here does.
find_package(ZLIB REQUIRED)
add_library(mxi_gemmi_mtz STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/gemmi/src/mtz.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/gemmi/src/gz.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/gemmi/src/sprintf.cpp)
target_include_directories(mxi_gemmi_mtz SYSTEM PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/gemmi/include
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/gemmi/third_party)
target_link_libraries(mxi_gemmi_mtz PUBLIC ZLIB::ZLIB mxi_gemmi)
target_compile_options(mxi_gemmi_mtz PRIVATE -w)

add_library(mxi_export OBJECT apps/mxi_export.cc)
target_link_libraries(mxi_export PRIVATE mxi mxi_gemmi_mtz)
target_compile_options(mxi_export PRIVATE -Wall -Wextra)

add_library(mxi_symmetry OBJECT apps/mxi_symmetry.cc)
target_link_libraries(mxi_symmetry PRIVATE mxi)
target_compile_options(mxi_symmetry PRIVATE -Wall -Wextra)
add_library(mxi_show OBJECT apps/mxi_show.cc)
target_link_libraries(mxi_show PRIVATE mxi)
target_compile_options(mxi_show PRIVATE -Wall -Wextra)

add_library(mxi_correlation OBJECT apps/mxi_correlation.cc)
target_link_libraries(mxi_correlation PRIVATE mxi)
target_compile_options(mxi_correlation PRIVATE -Wall -Wextra)
add_library(mxi_report OBJECT apps/mxi_report.cc)
target_link_libraries(mxi_report PRIVATE mxi)
target_compile_options(mxi_report PRIVATE -Wall -Wextra)

add_library(mxi_residuals OBJECT apps/mxi_residuals.cc)
target_link_libraries(mxi_residuals PRIVATE mxi)
target_compile_options(mxi_residuals PRIVATE -Wall -Wextra)

add_library(mxi_profile OBJECT apps/mxi_profile.cc)
target_link_libraries(mxi_profile PRIVATE mxi)
target_compile_options(mxi_profile PRIVATE -Wall -Wextra)

# Each predicted reflection's foreground region, for dials.image_viewer (was
# mxi_mask, whose name the mask editor took, 22 October).
add_library(mxi_foreground OBJECT apps/mxi_foreground.cc)
target_link_libraries(mxi_foreground PRIVATE mxi)
target_compile_options(mxi_foreground PRIVATE -Wall -Wextra)

add_library(mxi_background OBJECT apps/mxi_background.cc)
target_link_libraries(mxi_background PRIVATE mxi)
target_compile_options(mxi_background PRIVATE -Wall -Wextra)

add_library(mxi_forward OBJECT apps/mxi_forward.cc)
target_link_libraries(mxi_forward PRIVATE mxi)
target_compile_options(mxi_forward PRIVATE -Wall -Wextra)

add_library(mxi_convert OBJECT apps/mxi_convert.cc)
target_link_libraries(mxi_convert PRIVATE mxi)
target_compile_options(mxi_convert PRIVATE -Wall -Wextra)
add_library(mxi_combine OBJECT apps/mxi_combine.cc)
target_link_libraries(mxi_combine PRIVATE mxi)
target_compile_options(mxi_combine PRIVATE -Wall -Wextra)
add_library(mxi_split OBJECT apps/mxi_split.cc)
target_link_libraries(mxi_split PRIVATE mxi)
target_compile_options(mxi_split PRIVATE -Wall -Wextra)
# mxi_view, in Python (napari, optional): copied beside the programs, its
# version filled in, with the module mxi's Python tools share.
configure_file(tools/mxi_view.py ${CMAKE_CURRENT_BINARY_DIR}/mxi_view @ONLY
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
configure_file(tools/mxi_tools.py ${CMAKE_CURRENT_BINARY_DIR}/mxi_tools.py
  COPYONLY)
add_library(mxi_untrusted OBJECT apps/mxi_untrusted.cc)
target_link_libraries(mxi_untrusted PRIVATE mxi)
target_compile_options(mxi_untrusted PRIVATE -Wall -Wextra)

add_library(mxi_grid OBJECT apps/mxi_grid.cc)
target_link_libraries(mxi_grid PRIVATE mxi)
target_compile_options(mxi_grid PRIVATE -Wall -Wextra)

# The spot finder and the image reader, which need HDF5 and the bitshuffle
# submodule -- both required (below): without them nothing reads an image, and
# a build that quietly leaves out spot finding and integration is not one
# worth having (Graeme, v1.0).
#
# enable_testing() is called BEFORE it, not after: a subdirectory's add_test
# calls are only registered with the top-level CTest if testing was enabled
# before the subdirectory was entered. It sat after for a long time, which
# meant `ctest` from the top ran this project's tests and silently none of the
# spot finder's.
# SSE4.1 on x86, for the rounding instruction. Baseline x86-64 predates it, so
# std::rint compiles to a sequence rather than roundsd and the scoring loop in
# indexing runs three times slower: 3.5 nanoseconds a value against 1.2, and
# there are three values per reflection per triple of candidates scored.
# Anything from 2009 onwards has it; -DMXI_SSE41=OFF for a machine that does
# not.
#
# ARM needs nothing: the rounding instruction has always been there, which is
# why the same scoring loop was twenty-two times faster on an Apple machine
# than in this container, and why that difference first looked like the data.
option(MXI_SSE41 "Use SSE4.1 on x86, for the rounding instruction" ON)
if(MXI_SSE41 AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i[36]86")
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-msse4.1 MXI_COMPILER_TAKES_SSE41)
    if(MXI_COMPILER_TAKES_SSE41)
        add_compile_options(-msse4.1)
        message(STATUS "x86: SSE4.1 enabled for the rounding instruction")
    endif()
endif()

# AddressSanitizer and UndefinedBehaviorSanitizer, for running the tests under:
#
#     cmake -S . -B build-san -DMXI_SANITIZE=ON && cmake --build build-san
#     ctest --test-dir build-san
#
# -fno-sanitize-recover makes undefined behaviour ABORT rather than print and
# carry on. Without it UBSan reports and the test still passes, which is a
# sanitizer run that cannot fail. Global, so the spot finder and the libraries
# are instrumented as well as the tests: memory misused in a library and
# noticed in a test is only caught if both are built this way.
# MXI_SANITIZE itself is in cmake/spots.cmake: ON, or a list of sanitizers.


enable_testing()

# HDF5 is required (docs/rflx.md, Graeme's decision, October 2026): mxi's own
# files, .rflx, are HDF5, and so every program needs it, not only the spot
# finder and the importer.
# Every program's --version (src/args.cc).
target_compile_definitions(mxi PRIVATE MXI_VERSION="${PROJECT_VERSION}")
find_package(HDF5 REQUIRED COMPONENTS C)
target_include_directories(mxi SYSTEM PUBLIC ${HDF5_C_INCLUDE_DIRS})
target_link_libraries(mxi PUBLIC ${HDF5_C_LIBRARIES})
# bitshuffle (with LZ4) decompresses the images and compresses the .rflx's
# integer datasets: required, as HDF5 is.
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/bitshuffle/LICENSE")
    message(FATAL_ERROR
            "The bitshuffle submodule is not checked out, and without it nothing "
            "reads an image: git submodule update --init --recursive")
endif()
include(cmake/spots.cmake)
target_link_libraries(mxi PUBLIC bitshuffle)
target_compile_definitions(mxi PRIVATE MXI_HAVE_BITSHUFFLE)

# Profile fitting on a device, for mxi_integrate --gpu: CUDA's kernel where the
# spot finder's CUDA is on, and otherwise the stubs, which fit in single
# precision on the CPU instead.
if(MXI_CUDA)
    target_sources(mxi_integrate PRIVATE src/fit_device_cuda.cc src/fit_cuda.cu)
    target_link_libraries(mxi_integrate PRIVATE mxi_cuda_device)
elseif(MXI_METAL)
    # The shader compiled, linked and embedded as the spot finder's is, with the
    # same flags, under a name of its own.
    set(fit_metal_air "${CMAKE_CURRENT_BINARY_DIR}/fit_metal.air")
    set(fit_metallib "${CMAKE_CURRENT_BINARY_DIR}/fit_metal.metallib")
    set(fit_metallib_cc "${CMAKE_CURRENT_BINARY_DIR}/fit_metal_library.cc")
    add_custom_command(
        OUTPUT ${fit_metal_air}
        COMMAND ${MXI_XCRUN} -sdk macosx metal
                -std=metal3.0 ${mxi_metal_math} -Wall -Werror
                -c "${CMAKE_CURRENT_SOURCE_DIR}/src/fit_metal.metal"
                -o ${fit_metal_air}
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/fit_metal.metal"
        COMMENT "Compiling src/fit_metal.metal")
    add_custom_command(
        OUTPUT ${fit_metallib}
        COMMAND ${MXI_XCRUN} -sdk macosx metallib ${fit_metal_air} -o ${fit_metallib}
        DEPENDS ${fit_metal_air}
        COMMENT "Linking fit_metal.metallib")
    add_custom_command(
        OUTPUT ${fit_metallib_cc}
        COMMAND ${CMAKE_COMMAND} -DINPUT=${fit_metallib} -DOUTPUT=${fit_metallib_cc}
                -DSYMBOL=mxi_fit_metallib
                -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_metallib.cmake"
        DEPENDS ${fit_metallib} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_metallib.cmake"
        COMMENT "Embedding fit_metal.metallib")
    target_sources(mxi_integrate PRIVATE src/fit_device_metal.cc ${fit_metallib_cc})
    target_link_libraries(mxi_integrate PRIVATE mxi_metal)
    # SYSTEM: metal-cpp is not clean under -Wall -Wextra.
    target_include_directories(mxi_integrate SYSTEM PRIVATE ${METAL_CPP_INCLUDE_DIR})
    target_link_libraries(mxi_integrate PRIVATE ${METAL_FRAMEWORK} ${FOUNDATION_FRAMEWORK})
else()
    target_sources(mxi_integrate PRIVATE src/fit_device_none.cc)
endif()

add_executable(mxi_tests
  tests/main.cc
  tests/test_combine.cc
  tests/spots/test_dext_squares.cc
  tests/spots/test_dext_fused.cc
  tests/spots/test_signal_order.cc
  tests/spots/test_dials_spots.cc
  tests/spots/test_refl.cc
  tests/spots/test_refl_golden.cc
  tests/spots/test_expt.cc
  tests/spots/test_dials_end_to_end.cc
  tests/test_geometry.cc
  tests/test_real_geometry.cc
  tests/test_multi_axis.cc
  tests/test_index.cc
  tests/test_refine.cc
  tests/test_args.cc
  tests/test_derivatives.cc
  tests/test_precision.cc
  tests/test_expt_format.cc
  tests/test_background.cc
  tests/test_integrate.cc
  tests/test_reference.cc
  tests/test_refl.cc
  tests/test_shoebox.cc
  tests/test_summary.cc
  tests/test_postrefine.cc
  tests/test_symmetry.cc
  tests/test_lattice.cc
  tests/test_laue.cc
  tests/test_rflx.cc
  tests/test_cosym.cc
  tests/test_views.cc
  tests/test_processing.cc
  tests/test_two_theta.cc
  tests/test_correlation.cc
  tests/test_space_group.cc
  tests/test_resolution.cc
  tests/test_timing.cc
  tests/test_parallel.cc
  tests/test_fill_row.cc
  tests/test_scale_model.cc
  tests/test_scale.cc
  tests/test_mask.cc
  tests/test_profile_model.cc
)
target_link_libraries(mxi_tests PRIVATE mxi mxi_dext mxi_signal_order mxi_spots)
target_compile_options(mxi_tests PRIVATE -Wall -Wextra)
add_test(NAME mxi_tests COMMAND mxi_tests)

# One program, mxi, every other a subcommand of it: mxi find, mxi index ...
# (apps/mxi.cc). Each program above is an object library, its entry
# mxi_<name>_main, its own sources and settings kept; mxi links them all. One
# binary rather than 24, each of which carried its own HDF5 and library -- 156
# MB where one is some 8 -- and one name to install (Graeme, 28 October).
set(MXI_PROGRAMS
  mxi_background mxi_combine mxi_convert mxi_correlation mxi_export mxi_find
  mxi_foreground mxi_forward mxi_grid mxi_import mxi_index mxi_integrate
  mxi_max mxi_profile mxi_readtest mxi_refine mxi_report mxi_residuals
  mxi_scale mxi_show mxi_split mxi_symmetry mxi_untrusted)
# No source in two commands: they are one binary, and a source compiled into
# two is defined twice -- src/cuda_device.cu was, in mxi_find and
# mxi_integrate, and the CUDA release build failed to link (1.1.0). Shared code
# goes in a library (mxi_cuda_device, mxi_metal), linked once. Checked at
# configure, for whatever backend is configured.
set(mxi_seen_sources "")
foreach(program ${MXI_PROGRAMS})
  get_target_property(program_sources ${program} SOURCES)
  foreach(source ${program_sources})
    get_filename_component(source "${source}" ABSOLUTE)
    list(FIND mxi_seen_sources "${source}" seen)
    if(NOT seen EQUAL -1)
      list(GET mxi_seen_owners ${seen} owner)
      message(FATAL_ERROR "${source} is compiled into both ${owner} and "
              "${program}, which are one binary: put it in a library both link")
    endif()
    list(APPEND mxi_seen_sources "${source}")
    list(APPEND mxi_seen_owners ${program})
  endforeach()
endforeach()

# The target mxi_program, its file mxi: the target mxi is the library's.
add_executable(mxi_program apps/mxi.cc)
set_target_properties(mxi_program PROPERTIES OUTPUT_NAME mxi)
target_link_libraries(mxi_program PRIVATE ${MXI_PROGRAMS} mxi)
target_compile_definitions(mxi_program PRIVATE MXI_VERSION="${PROJECT_VERSION}")
target_compile_options(mxi_program PRIVATE -Wall -Wextra)
# In the build directory, each program's old name a link to mxi, which runs
# the program it is called by: the tests run them so. Not installed.
foreach(program ${MXI_PROGRAMS})
  file(CREATE_LINK mxi ${CMAKE_CURRENT_BINARY_DIR}/${program} SYMBOLIC)
endforeach()

# Installing: mxi, with mxi_view and the module it shares beside it --
#   cmake --install build --prefix ~/.local
# -- which then wants to be on the PATH. The library is linked in; nothing
# else is needed beside it (the Metal libraries are compiled in).
install(TARGETS mxi_program RUNTIME DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mxi_view DESTINATION bin)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mxi_tools.py DESTINATION bin)
