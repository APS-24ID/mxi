# Making a release

A release is made by pushing a tag: GitHub builds it, tests it and publishes
it (`.github/workflows/release.yml`). Nothing is published unless every build
and its tests pass.

```sh
# The version in its five places (CMakeLists.txt, vcpkg.json, python/
# pyproject.toml, python/src/mxeq/__init__.py, docs/spots.md's example), and
# RELEASE_NOTES.md's section for it, then:
git tag -a v1.0.7 -m "mxi 1.0.7"
git push origin v1.0.7
```

Push a candidate first -- `v1.0.7-rc1` -- to see the whole of it work: a tag
with `-rc`, `-alpha` or `-beta` makes a pre-release, its notes 1.0.7's.

## What it makes

| asset | what |
|---|---|
| `mxi-1.0.7-source.tar.gz` | the complete source: the repository and its submodules, gemmi and bitshuffle, at the commits the tag pins. GitHub's own source archives leave the submodules out, and without them nothing builds; this is also the source for a conda-forge recipe |
| `mxi-1.0.7-linux-x86_64-cuda.tar.gz` | every program, built on Ubuntu 22.04 with CUDA for Volta to Hopper (V100 to H100 and GH200) |
| `mxi-1.0.7-macos-arm64.tar.gz` | every program, for Apple silicon, with Metal |
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
machines as before.

The Metal build downloads Apple's metal-cpp headers from the URL at the top of
the workflow (`METAL_CPP_URL`); if Apple renames the zip, the current name is
on https://developer.apple.com/metal/cpp/.

## The scripts, by hand

```sh
scripts/source_archive.sh v1.0.7         # mxi-1.0.7-source.tar.gz
scripts/release_notes.py v1.0.7 > notes.md
scripts/package_release.sh build v1.0.7 linux-x86_64-cuda
```

A static build like the release's, here, needs vcpkg:

```sh
git clone https://github.com/microsoft/vcpkg && vcpkg/bootstrap-vcpkg.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake \
      -DVCPKG_TARGET_TRIPLET=x64-linux-release
```
