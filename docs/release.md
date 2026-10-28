# Making a release

A release is made by pushing a tag: GitHub builds it, tests it and publishes
it (`.github/workflows/release.yml`). Nothing is published unless every build
and its tests pass. Actions minutes are free for a public repository.

```sh
# The version in its five places (CMakeLists.txt, vcpkg.json, python/
# pyproject.toml, python/src/mxeq/__init__.py, docs/spots.md's example), and
# RELEASE_NOTES.md's section for it, then:
git tag -a v1.0.7 -m "mxi 1.0.7"
git push origin v1.0.7
```

It can be run by hand too, for a tag pushed already -- to make a release again
after a failure, say -- from the Actions page (release, Run workflow) or
`gh workflow run release.yml -f tag=v1.0.7`. A tag whose release exists must
have it deleted first (`gh release delete v1.0.7`): GitHub keeps one release a
tag.

Push a candidate first -- `v1.0.7-rc1` -- to see the whole of it work: a tag
with `-rc`, `-alpha` or `-beta` makes a pre-release. Its notes are 1.0.7's
section where there is one yet, and otherwise the "Unreleased" section, said
to be a candidate's. The release itself must have its section: without one
nothing is published.

## What it makes

| asset | what |
|---|---|
| `mxi-1.0.7-source.tar.gz` | the complete source: the repository and its submodules, gemmi and bitshuffle, at the commits the tag pins. GitHub's own source archives leave the submodules out, and without them nothing builds; this is also the source for a conda-forge recipe |
| `mxi-1.0.7-linux-x86_64.tar.gz` | every program, built on Ubuntu 22.04, for the CPU: it runs anywhere, no CUDA driver needed |
| `mxi-1.0.7-linux-x86_64-cuda.tar.gz` | every program, built on Ubuntu 22.04 with CUDA for Volta to Hopper (V100 to H100 and GH200) |
| `mxi-1.0.7-macos-arm64.tar.gz` | every program, for Apple silicon, with Metal: macOS 13 or later, its static libraries too (an overlay triplet, `cmake/triplets`); any API newer than 13 used without a check fails the build, so that the binaries cannot fail on an older macOS |
| `mxeq-1.0.7-py3-none-any.whl` | the referee, `pip install` it |

Each binary archive is `bin/` -- the programs, `mxi_view` and its module
`mxi_tools.py` -- with the licence, the README and the release notes. Unpack
it and put its `bin` on the PATH. The programs need nothing installed: HDF5
and zlib are linked in, from vcpkg (`vcpkg.json`), and on Linux the C++
runtime too, so that they run on distributions older than the build's, RHEL
among them; the CUDA driver is the system's. `scripts/package_release.sh`
refuses to package a build that would need anything else. `mxi_view` needs
napari, as always (`docs/mask.md`).

The archive is some 50 MB: each program carries its own HDF5.

## Licences

mxi is BSD-3-Clause, and the release builds stay so: FFTW, which is GPL, is
opt-in (`-DMXI_FFTW=ON`), the workflow says `-DMXI_FFTW=OFF`, and
`scripts/package_release.sh` refuses a program containing it. The code linked
into the programs comes with its notices, in each archive's `licenses/`:
gemmi (MPL-2.0, its source in the source archive), bitshuffle (MIT) and the
LZ4 within it (BSD-2), from the submodules; HDF5, libaec and zlib, from
vcpkg's copyright files; and NVIDIA's CUDA runtime in the CUDA build. A notice
missing fails the package -- a new library linked in needs its own added.

## macOS

The binaries are not signed, so macOS blocks them once downloaded. Clear the
quarantine flag once, after unpacking -- the release notes say so:

```sh
xattr -dr com.apple.quarantine mxi-1.0.7-macos-arm64
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
scripts/source_archive.sh v1.0.7         # mxi-1.0.7-source.tar.gz
scripts/release_notes.py v1.0.7 > notes.md
scripts/package_release.sh build v1.0.7 linux-x86_64
```

A static build like the release's, here, needs vcpkg:

```sh
git clone https://github.com/microsoft/vcpkg && vcpkg/bootstrap-vcpkg.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake \
      -DVCPKG_TARGET_TRIPLET=x64-linux-release
```
