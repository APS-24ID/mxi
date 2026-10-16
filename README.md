# mxi

```
░▒▓ █▀▄▀█ ▀▄▀ █
░▒▓ █ ▀ █ ▄▀▄ █
```

Macromolecular crystallography data reduction for rotation data: from NXmx
images to scaled, unmerged intensities -- spot finding, indexing,
refinement, integration, symmetry, scaling and clustering of many crystals --
in C++, with spot finding and profile fitting on the GPU through CUDA or Metal.
It follows DIALS's methods, reads and writes DIALS's files, and any step can be
exchanged for DIALS's own.

* **How it works:** `docs/algorithms.md`, each stage's methods and sources.
* **How it differs from DIALS, and why:** `docs/review.md`.
* **Every document:** `docs/README.md`. **What is open:** `docs/outstanding.md`.

## Building

Requires CMake 3.20 or later, a C++20 compiler, and HDF5's C library; and two
submodules, gemmi (space groups) and bitshuffle (the images' compression). The
configure stops if any is missing.

```sh
sudo apt-get install libhdf5-dev        # Debian and Ubuntu
brew install hdf5                       # macOS
# or: conda install -c conda-forge hdf5

git submodule update --init --recursive
cmake -S . -B build && cmake --build build
ctest --test-dir build
```

Options:

| option | |
| --- | --- |
| `-DSPOTFINDER_CUDA=ON` | spot finding and profile fitting on an NVIDIA GPU; `-DCMAKE_CUDA_ARCHITECTURES="80;90"` for cards other than this machine's |
| `-DSPOTFINDER_METAL=ON -DMETAL_CPP_DIR=DIR` | the same on Apple silicon; `DIR` holds `Metal/Metal.hpp`, Apple's metal-cpp from <https://developer.apple.com/metal/cpp/>. One device backend at a time. |
| `-DMXI_FFTW=ON` | indexing's transform by FFTW, about three times the built-in's speed; GPL, so opt-in. `-DFFTW3_LIBRARY=... -DFFTW3_INCLUDE_DIR=...` where it is not found. |
| `-DSPOTFINDER_AVX2=OFF` | a binary for x86 CPUs without AVX2 |
| `-DMXI_SSE41=OFF` | for x86 CPUs without SSE4.1 (before 2009) |
| `-DMXI_SANITIZE=ON` | AddressSanitizer and UndefinedBehaviorSanitizer, for the tests |

## Usage

### The main sequence

```sh
mxi_import    master.nxs              # the experiment list, a sweep each master  -> imported.rflx
mxi_find      imported.rflx           # strong spots on every image                -> strong.rflx
mxi_index     strong.rflx             # the lattice, by 3D FFT                     -> indexed.rflx
mxi_refine    indexed.rflx            # beam, detector, scan-varying crystal       -> refined.rflx
mxi_integrate refined.rflx            # summation and profile fitting              -> integrated.rflx
mxi_symmetry  integrated.rflx         # the Laue group and space group, reindexed  -> symmetrized.rflx
mxi_scale     symmetrized.rflx        # scaled, with merging statistics            -> scaled.rflx
mxi_export    scaled.rflx             # an unmerged MTZ                            -> scaled.mtz
```

Each step reads the last one's file and writes its own, under the name shown
unless `-o` names another: a `.rflx`, mxi's own format -- HDF5 holding the
experiment list, the reflections and the processing's summary
(`docs/rflx.md`). Every program also writes what it prints to
`mxi_<program>.log` and takes `--help` and `--version`; the main sequence's
programs take `--timing`, where the time went, and use every core;
`-j N` or `--threads N` sets fewer, where a program has the choice. `mxi_find --gpu` and `mxi_integrate --gpu`
use a GPU where the build has one.

**Several sweeps of one crystal:** name every master to `mxi_import`; each
program takes them all, the sweeps indexed in one basis, each its own crystal
and scale model (`docs/multi-sweep.md`).

**Several crystals in one beam:** `mxi_index --max-lattices N`; each lattice
is refined, integrated and scaled as its own experiment, the detector shared,
cosym bringing them to one setting in `mxi_symmetry`
(`docs/multi-lattice.md`).

**Many crystals, a data set each:** after scaling them together,
`mxi_correlation scaled.rflx --write-clusters` finds the groups that belong
together by their intensities, writes each as its own `.rflx` to scale alone,
and an HTML report (`docs/correlation.md`).

The options most often wanted:

| program | options |
| --- | --- |
| `mxi_find` | `-g`/`--gpu`; `--min-spot-size N`; `--no-shoeboxes` |
| `mxi_index` | `--d-min D`; `--max-cell A`; `--max-lattices N`; `--no-joint` sweeps indexed apart |
| `mxi_refine` | `--static` one crystal setting for the scan; `--scan-varying N` N control points; `--shared-crystal` |
| `mxi_integrate` | `-g`/`--gpu`; `--d-min D`, `--d-max D`; `--postrefine`; `--summation-only` |
| `mxi_symmetry` | `--max-delta D` the lattice's tolerance; `--no-cosym` |
| `mxi_scale` | `--d-min-auto`, `--d-min D`, `--d-max D`; `--space-group NAME`; `--anomalous`; `--absorption-level low|medium|high`; `--error-model-fit regression|individual` |
| `mxi_correlation` | `--write-clusters`, `--prefix P`; `--html FILE`; `--min-samples N`; `--dimensions N` |
| `mxi_export` | `--partiality-threshold P`; `--min-isigi S`; `--d-min D` |

### The other tools

| program | |
| --- | --- |
| `mxi_show` | what a file holds: the models, the reflections by what has been done to them, each program's summary (`docs/output.md`) |
| `mxi_report` | an HTML report of a file's processing, every step's figures and plots (`docs/output.md`) |
| `mxi_convert` | between `.rflx` and DIALS's `.expt` and `.refl` (below) |
| `mxi_max` | the largest count on a valid pixel, and whether the images fit in 16 bits (for `mxi_find --gpu-force`) |
| `mxi_profile` | the profile model estimated from strong spots, and how much of each spot it holds |
| `mxi_residuals` | refinement's residuals by resolution, detector module, or anything else |
| `mxi_mask` | shoeboxes marking the integration region, for the image viewer |
| `mxi_grid` | spot density in Kabsch's frame, and spot widths across the detector |
| `mxi_forward` | the model rendered onto the pixels, against the data |
| `mxi_background` | the robust background fitted to pixel values |
| `mxi_readtest` | how fast the images read and decompress |

## Working with DIALS

**The files.** Every program reads DIALS's `.expt` and `.refl` as well as a
`.rflx`, and writes them on request (`--output-expt`, `--output-refl`, or `-o`
naming a `.refl`). `mxi_convert` goes either way: `mxi_convert indexed.rflx`
writes `indexed.expt` and `indexed.refl`; `mxi_convert indexed.expt
indexed.refl` the `.rflx`. Two things DIALS cannot read stay behind: each
sweep's own crystal-less experiment, which mxi carries through indexing
(dials/dials#1029; `--keep-sweeps` keeps them), and mxi's scale model.

**Swapping steps.** Any step can be DIALS's: `dials.symmetry` and `dials.scale`
take mxi's integrated reflections, mxi's programs take DIALS's at any stage,
and `dials.image_viewer`, `dials.export` and `dials.merge` read mxi's output.

**The equivalents:**

| mxi | DIALS |
| --- | --- |
| `mxi_import` | `dials.import` |
| `mxi_find` | `dials.find_spots` |
| `mxi_index` | `dials.index` |
| `mxi_refine` | `dials.refine` |
| `mxi_integrate` | `dials.integrate` |
| `mxi_symmetry` | `dials.symmetry`, with `dials.cosym` for several crystals |
| `mxi_scale` | `dials.scale`, with `dials.estimate_resolution` and `dials.two_theta_refine` |
| `mxi_correlation` | `dials.correlation_matrix` |
| `mxi_export` | `dials.export` |
| `mxi_show`, `mxi_report` | `dials.show`, `dials.report` |

Where they differ, and by how much, is in `docs/review.md`.

## mxeq: comparisons and regression runs

`mxeq` is a Python package that judges mxi's output, independently of it: it
reads files and imports nothing from `src/`.

```sh
pip install -e python                    # mxeq
pip install -e "python[plots,fixtures]"  # with its plots and test-data makers
```

**Comparing with DIALS** (`python/README.md`): `mxeq check` compares two
pipelines at a boundary; `compare-expt` two experiment lists; `trend`,
`explain` and `disagree` where two integrations differ, and why; `unique` two
scaled data sets a reflection at a time; `equivalents` integration's bias
against the data's own symmetry equivalents, with no second program.

**Regression runs** (`docs/regression.md`): `mxeq regress` builds a commit,
runs the chain over a set of real data sets, and keeps every figure in a
history, so that runs, commits and machines can be compared:

```sh
mxeq regress init                        # datasets.toml and build.toml, to fill in
mxeq regress run                         # the current commit over every data set
mxeq regress compare --against reference # what changed, past its tolerance
```

## Tests

`ctest --test-dir build` runs the C++ tests. `pytest python/tests` runs
mxeq's, and those of the programs that need real data, given the programs and
data through `MXI_*` environment variables; without them those tests are
skipped (`python/README.md`).
