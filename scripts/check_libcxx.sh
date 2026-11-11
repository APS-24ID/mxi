#!/bin/sh
# Every source and test compiled against libc++, Apple's C++ library, on
# Linux (clang++ -stdlib=libc++, syntax only): GCC's library includes more
# than the standard asks, so what compiles with it may not on a Mac. Not
# what an "#if defined(__linux__)" hides, which is true here whichever
# library -- <numbers>, put under one in tests/test_combine.cc (7 November),
# compiled here and not on Graeme's Mac; python/tests/test_includes.py
# guards that. The Metal, CUDA and FFTW sources left out.
#
#   scripts/check_libcxx.sh     (apt install clang libc++-dev libc++abi-dev)
set -eu
cd "$(dirname "$0")/.."
command -v clang++ >/dev/null || { echo "no clang++: apt install clang libc++-dev" >&2; exit 2; }
h5=$(dirname "$(find /usr/include -name hdf5.h | head -1)")
files=0 failed=0
for f in src/*.cc apps/*.cc tests/*.cc tests/spots/*.cc tests/tools/*.cc; do
    [ -f "$f" ] || continue
    case $f in *metal*|*cuda*|*fft_fftw*|*dext_gpu*) continue ;; esac
    files=$((files + 1))
    if ! out=$(clang++ -std=c++20 -stdlib=libc++ -fsyntax-only -w -DMXI_VERSION='"check"' \
            -Isrc -Iapps -Itests -Ithird_party/gemmi/include -Ithird_party/bitshuffle/src \
            -Ithird_party/bitshuffle/lz4 -I"$h5" "$f" 2>&1); then
        failed=$((failed + 1))
        echo "$f: $(echo "$out" | grep -m1 error)" >&2
    fi
done
echo "compiled against libc++: $files files, $failed failed"
[ "$failed" -eq 0 ]
