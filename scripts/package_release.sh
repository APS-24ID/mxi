#!/bin/sh
# A release's binaries from a build: mxi, with mxi_view and the module it
# shares beside it, installed into mxi-VERSION-PLATFORM/bin, with the
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
# FFTW, where the build has it (MXI_FFTW): its GPL, from vcpkg's copyright
# file -- or a system's, building by hand -- and the README saying what that
# makes of the binary.
fftw=no
if grep -q '^MXI_FFTW:BOOL=ON' "$build/CMakeCache.txt" 2>/dev/null; then
    fftw=yes
    found=$(find "$installed" -path "*/share/fftw3/copyright" 2>/dev/null | head -1)
    [ -n "$found" ] || found=$(ls /usr/share/doc/libfftw3-dev/copyright 2>/dev/null | head -1)
    if [ -z "$found" ]; then
        echo "the build has FFTW but no licence for it was found" >&2
        exit 1
    fi
    cp "$found" "$lic/fftw.txt"
fi
cat > "$lic/README" <<'NOTE'
mxi is BSD-3-Clause (mxi.txt). The program also contains, linked in: gemmi
(MPL-2.0; its source is in the release's mxi-<version>-source.tar.gz),
bitshuffle (MIT) and the LZ4 within it (BSD-2), HDF5, libaec and zlib, each
under the licence here; and, in the CUDA build, NVIDIA's CUDA runtime.
NOTE
if [ "$fftw" = yes ]; then
    cat >> "$lic/README" <<'NOTE'

It also contains FFTW (fftw.txt), for indexing's transform, which is under the
GNU General Public License, version 2 or later. The program as distributed
here is therefore under the GPL's terms as a whole; mxi's own source remains
BSD-3-Clause, and every part's source is open: mxi's and its submodules' in
the release's source archive, FFTW's at https://www.fftw.org.
NOTE
fi
# Self-contained: no HDF5, zlib or vcpkg library to be found at run time. On
# Linux libcuda (the driver) is the system's, as it must be.
bad=0
for p in "$name"/bin/mxi "$name"/bin/mxi_*; do
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

# FFTW in the program exactly when the build says so, so that its licence
# goes with it whenever it is there.
if grep -a -q 'fftw_plan' "$name/bin/mxi"; then
    [ "$fftw" = yes ] || { echo "bin/mxi contains FFTW the build does not record" >&2; exit 1; }
fi
# One program, mxi, every other its subcommand.
[ -x "$name/bin/mxi" ] || { echo "no bin/mxi" >&2; exit 1; }
count=$("$name/bin/mxi" help | grep -c '^  [a-z]')
tar -czf "$name.tar.gz" "$name"
echo "$name.tar.gz: mxi and its $count commands"
