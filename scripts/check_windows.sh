#!/bin/sh
# What can be checked of a Windows build without Windows (docs/windows.md):
# every source compiled for Windows by MinGW-w64, syntax only, warnings on --
# HDF5's headers for Linux standing in, a stub for glibc's features.h -- and
# the platform layer's Windows half built and run under Wine, where Wine is.
#
#   scripts/check_windows.sh     (apt install g++-mingw-w64-x86-64-posix wine64)
set -eu
cd "$(dirname "$0")/.."
cxx=x86_64-w64-mingw32-g++
command -v "$cxx" >/dev/null || { echo "no $cxx: apt install g++-mingw-w64-x86-64-posix" >&2; exit 2; }
h5=$(dirname "$(find /usr/include -name hdf5.h | head -1)")
stub=$(mktemp -d)
trap 'rm -rf "$stub"' EXIT
touch "$stub/features.h"
flags="-std=c++20 -fsyntax-only -Wall -Wextra -D__USE_MINGW_ANSI_STDIO=1 -DMXI_VERSION=\"check\" -Isrc -Iapps \
  -Ithird_party/gemmi/include -Ithird_party/bitshuffle/src -Ithird_party/bitshuffle/lz4 \
  -I$stub -I$h5"
files=0 failed=0 warned=0
for f in src/*.cc apps/*.cc; do
    case $f in *metal*|*cuda*|*fft_fftw*|*readtest*) continue ;; esac
    files=$((files + 1))
    # shellcheck disable=SC2086
    if ! out=$("$cxx" $flags "$f" 2>&1); then
        failed=$((failed + 1)); echo "$out" | grep -m3 error >&2
    fi
    n=$(echo "$out" | grep "warning:" | grep -v -c "third_party\|/usr/include\|$stub" || true)
    warned=$((warned + n))
done
echo "compiled for Windows: $files sources, $failed failed, $warned warnings in mxi's code"
exe="$stub/platform_windows.exe"
"$cxx" -std=c++20 -O1 -Wall -Wextra -D__USE_MINGW_ANSI_STDIO=1 -Isrc \
    tests/windows/platform_windows.cc -static -o "$exe"
wine=$(command -v wine64 || command -v wine || ls /usr/lib/wine/wine64 2>/dev/null || true)
if [ -n "$wine" ]; then
    (cd "$stub" && WINEDEBUG=-all WINEPREFIX="$stub/prefix" "$wine" "$exe" 2>/dev/null | tr -d '\r' | grep -v "^wine:")
else
    echo "no Wine: the platform test built, not run"
fi
[ "$failed" -eq 0 ] && [ "$warned" -eq 0 ]
