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
count=0
for p in "$name"/bin/mxi_*; do
    [ -x "$p" ] && count=$((count + 1))
done
tar -czf "$name.tar.gz" "$name"
echo "$name.tar.gz: $count programs"
