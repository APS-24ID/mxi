# Making a release

A release is made by pushing a tag: GitHub builds it, tests it and publishes
it (`.github/workflows/release.yml`). Nothing is published unless every build
and its tests pass. Actions minutes are free for a public repository.

```sh
# The version in its five places (CMakeLists.txt, vcpkg.json, python/
# pyproject.toml, python/src/mxeq/__init__.py, docs/spots.md's example), and
# RELEASE_NOTES.md's section for it, then:
git tag -a v1.1.1 -m "mxi 1.1.1"
git push origin v1.1.1
```

It can be run by hand too, for a tag pushed already -- to make a release again
after a failure, say -- from the Actions page (release, Run workflow) or
`gh workflow run release.yml -f tag=v1.1.1`. A tag whose release exists must
have it deleted first (`gh release delete v1.1.1`): GitHub keeps one release a
tag.

Push a candidate first -- `v1.1.1-rc1` -- to see the whole of it work: a tag
with `-rc`, `-alpha` or `-beta` makes a pre-release. Its notes are 1.1.1's
section where there is one yet, and otherwise the "Unreleased" section, said
to be a candidate's. The release itself must have its section: without one
nothing is published.

## What it makes

| asset | what |
|---|---|
| `mxi-1.1.1-source.tar.gz` | the complete source: the repository and its submodules, gemmi and bitshuffle, at the commits the tag pins. GitHub's own source archives leave the submodules out, and without them nothing builds; this is also the source for a conda-forge recipe |
| `mxi-1.1.1-linux-x86_64.tar.gz` | `mxi`, built on Ubuntu 22.04, for the CPU: it runs anywhere, no CUDA driver needed |
| `mxi-1.1.1-linux-x86_64-cuda.tar.gz` | `mxi`, built on Ubuntu 22.04 with CUDA for Volta to Hopper (V100 to H100 and GH200) |
| `mxi-1.1.1-linux-arm64.tar.gz` | `mxi` for ARM64 Linux, built on GitHub's ARM runner with Ubuntu 22.04, for the CPU: a Raspberry Pi 4 or 5 with a 64-bit OS, ARM servers |
| `mxi-1.1.1-linux-arm64-cuda.tar.gz` | the same with CUDA, NVIDIA's toolkit for ARM servers (sbsa), for Hopper -- the GH200 -- with its PTX for later GPUs; not Jetson, whose toolkit is another |
| `mxi-1.1.1-macos-arm64.tar.gz` | `mxi`, for Apple silicon, with Metal: macOS 13 or later, its static libraries too (an overlay triplet, `cmake/triplets`); any API newer than 13 used without a check fails the build, so that the binaries cannot fail on an older macOS |
| `mxeq-1.1.1-py3-none-any.whl` | the referee, `pip install` it |

Each binary archive is `bin/` -- `mxi`, every step a command of it, with `mxi
view`'s Python and its module `mxi_tools.py` beside it -- with the licence,
the README and the release notes. Unpack it and put its `bin` on the PATH. The
program needs nothing installed: HDF5 and zlib are linked in, from vcpkg
(`vcpkg.json`), and on Linux the C++ runtime too, so that it runs on
distributions older than the build's, RHEL among them; the CUDA driver is the
system's. `scripts/package_release.sh` refuses to package a build that would
need anything else. `mxi view` needs napari, as always (`docs/mask.md`).

One program, not 24 each with its own HDF5: the archive is a few megabytes, where 24 programs had made it some 50.

## Licences

mxi's source is BSD-3-Clause. The release builds include FFTW, for indexing's
transform: a vcpkg feature, `fftw` (`vcpkg.json`), asked for with
`-DVCPKG_MANIFEST_FEATURES=fftw -DMXI_FFTW=ON`. FFTW is GPL, version 2 or
later, so the program as released is distributed under the GPL's terms as a
whole -- a modest commitment, every part's source being open. A build without
FFTW (the default, `-DMXI_FFTW=OFF`) uses the built-in transform and is BSD
with the notices below.

The code linked into the program comes with its notices, in each archive's
`licenses/`: gemmi (MPL-2.0, its source in the source archive), bitshuffle
(MIT) and the LZ4 within it (BSD-2), from the submodules; HDF5, libaec, zlib
and FFTW, from vcpkg's copyright files; NVIDIA's CUDA runtime in the CUDA
build; and a README saying which is which, and what FFTW makes of the whole.
`scripts/package_release.sh` fails if a notice is missing, or if the program
contains FFTW the build does not record -- a new library linked in needs its
own added.

## macOS

The binaries are not signed, so macOS blocks them once downloaded. Clear the
quarantine flag once, after unpacking -- the release notes say so:

```sh
xattr -dr com.apple.quarantine mxi-1.1.1-macos-arm64
```

Signing and notarising would need an Apple Developer ID and its credentials
as the repository's secrets.

## What the builds test, and what they cannot

Each build runs ctest and the Python suite, the tests that need no external
data, before its archive is kept. Hosted runners have no GPU to rely on, so
the device tests (ctest's `dext_gpu`) are left out: the CUDA and Metal builds
are released tested on their CPU paths, their device paths tested on real
machines as before; the CPU build is tested whole.

The Metal build downloads Apple's metal-cpp headers from the URL at the top of
the workflow (`METAL_CPP_URL`); if Apple renames the zip, the current name is
on https://developer.apple.com/metal/cpp/.

## The scripts, by hand

```sh
scripts/source_archive.sh v1.1.1         # mxi-1.1.1-source.tar.gz
scripts/release_notes.py v1.1.1 > notes.md
scripts/package_release.sh build v1.1.1 linux-x86_64
```

A static build like the release's, here, needs vcpkg:

```sh
git clone https://github.com/microsoft/vcpkg && vcpkg/bootstrap-vcpkg.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake \
      -DVCPKG_TARGET_TRIPLET=x64-linux-release
```
