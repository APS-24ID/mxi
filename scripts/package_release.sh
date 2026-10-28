#!/bin/sh
# A release's binaries from a build: every program, mxi_view and the module
# it shares with them, installed into mxi-VERSION-PLATFORM/bin, with the
# licence, the README and the release notes; checked to need nothing but the
# system's own libraries -- HDF5 linked in, from vcpkg -- and archived.
#
#   scripts/package_release.sh build v1.0.6 linux-x86_64-cuda
set -eu
build=${1:?usage: package_release.sh BUILD TAG PLATFORM}
tag=${2:?}
platform=${3:?}
name="mxi-${tag#v}-${platform}"
root=$(cd "$(dirname "$0")/.." && pwd)
rm -rf "$name"
# Stripped: each program carries HDF5 within it, and its symbols are not needed.
cmake --install "$build" --prefix "$name" --strip
cp "$root/LICENSE" "$root/README.md" "$root/RELEASE_NOTES.md" "$name/"

# The notices of the code built into the programs, which their licences ask
# to go with binary copies: gemmi (MPL-2.0: its source is in the release's
# source archive), bitshuffle (MIT) and the LZ4 within it (BSD-2), from the
# submodules; HDF5, libaec and zlib from vcpkg's copyright files beside the
# libraries it built; NVIDIA's CUDA runtime, linked statically, where built
# with CUDA. A notice missing fails the package.
lic="$name/licenses"
mkdir -p "$lic"
cp "$root/LICENSE" "$lic/mxi.txt"
cp "$root/third_party/gemmi/LICENSE.txt" "$lic/gemmi.txt"
cp "$root/third_party/bitshuffle/LICENSE" "$lic/bitshuffle.txt"
cp "$root/third_party/bitshuffle/lz4/LICENSE" "$lic/lz4.txt"
installed=${VCPKG_INSTALLED_DIR:-$build/vcpkg_installed}
for package in hdf5 libaec zlib; do
    found=$(find "$installed" -path "*/share/$package/copyright" 2>/dev/null | head -1)
    if [ -z "$found" ]; then
        echo "no copyright for $package under $installed" >&2
        exit 1
    fi
    cp "$found" "$lic/$package.txt"
done
if grep -q '^MXI_CUDA:BOOL=ON' "$build/CMakeCache.txt" 2>/dev/null; then
    cat > "$lic/nvidia-cuda-runtime.txt" <<'NOTE'
The programs are linked statically with NVIDIA's CUDA runtime library
(cudart), which is redistributed under the NVIDIA CUDA Toolkit End User
License Agreement: https://docs.nvidia.com/cuda/eula/
NOTE
fi
cat > "$lic/README" <<'NOTE'
mxi is BSD-3-Clause (mxi.txt). The programs also contain, linked in: gemmi
(MPL-2.0; its source is in the release's mxi-<version>-source.tar.gz),
bitshuffle (MIT) and the LZ4 within it (BSD-2), HDF5, libaec and zlib, each
under the licence here; and, in the CUDA build, NVIDIA's CUDA runtime. FFTW,
which is GPL, is not used: the built-in transform is.
NOTE
# Self-contained: no HDF5, zlib or vcpkg library to be found at run time. On
# Linux libcuda (the driver) is the system's, as it must be.
bad=0
for p in "$name"/bin/mxi_*; do
    [ -f "$p" ] && [ -x "$p" ] || continue
    head -c 4 "$p" | grep -q '#!' && continue   # mxi_view, Python
    if [ "$(uname)" = Darwin ]; then
        libs=$(otool -L "$p" | tail -n +2 | awk '{print $1}' | grep -v -E '^(/usr/lib/|/System/)' || true)
    else
        libs=$(ldd "$p" | awk '{print $1}' | grep -i -E 'hdf5|libz\.|libsz|vcpkg' || true)
    fi
    if [ -n "$libs" ]; then
        echo "$p needs: $libs" >&2
        bad=1
    fi
done
[ "$bad" = 0 ] || { echo "not self-contained" >&2; exit 1; }

# No FFTW, which is GPL: the release builds are BSD with the notices above.
for p in "$name"/bin/mxi_*; do
    if [ -x "$p" ] && grep -a -q -i 'fftw_plan' "$p"; then
        echo "$p contains FFTW, which is GPL: build without -DMXI_FFTW=ON" >&2
        exit 1
    fi
done
count=0
for p in "$name"/bin/mxi_*; do
    [ -x "$p" ] && count=$((count + 1))
done
tar -czf "$name.tar.gz" "$name"
echo "$name.tar.gz: $count programs"
