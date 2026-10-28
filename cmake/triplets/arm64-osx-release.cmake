# vcpkg's arm64-osx-release, with a minimum macOS: the release's static HDF5
# and zlib built for macOS 13 and later, as mxi is (the release workflow's
# CMAKE_OSX_DEPLOYMENT_TARGET). Without it they took the build machine's own
# version, and the binaries would have refused anything older than the
# runner's macOS 15. An overlay: -DVCPKG_OVERLAY_TRIPLETS=cmake/triplets.
set(VCPKG_TARGET_ARCHITECTURE arm64)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE static)
set(VCPKG_CMAKE_SYSTEM_NAME Darwin)
set(VCPKG_OSX_ARCHITECTURES arm64)
set(VCPKG_BUILD_TYPE release)
set(VCPKG_OSX_DEPLOYMENT_TARGET 13.0)
