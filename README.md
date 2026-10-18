# mxi

```
░▒▓ █▀▄▀█ ▀▄▀ █
░▒▓ █ ▀ █ ▄▀▄ █
```

mxi is a data processing package for NXmx format macromolecular crystallography
(MX) data: from the diffraction images to scaled, unmerged intensities ready for
structure solution -- spot finding, indexing, refinement, integration, symmetry,
scaling and the clustering of data from many crystals.

**It is a new, standalone package.** mxi is a completely new project, written from
scratch in C++. It contains no DIALS code and does not need DIALS installed.
It is compatible with DIALS: it follows DIALS' methods, reads and writes
DIALS' files, and any step can be exchanged for DIALS' own (see [Working with
DIALS](#working-with-dials)).

**Its purpose is almost instant data processing.** mxi aims to process data as
DIALS does, with results of similar quality, an order of magnitude more quickly
-- fast enough to keep pace with modern detectors, and to give results while the
crystal is still mounted. Spot finding and integration can run on a GPU,
through CUDA or Apple's Metal.

**It reads NXmx format data only.**
[NXmx](https://journals.iucr.org/m/issues/2020/05/00/ti5018/index.html)
is the community standard format for X-ray
diffraction images: a layout, built on the HDF5 file format, that keeps the
images together with everything needed to interpret them -- the X-ray beam,
the detector's geometry, how the crystal was rotated -- so that any program
can read data from any beamline the same way. Modern detectors, such as
DECTRIS's EIGER, write it directly; older formats can be converted, for example
CBF by [cbf2nxmx](https://github.com/graeme-winter/cbf2nxmx).

* **How it works:** [`docs/algorithms.md`](docs/algorithms.md), each stage's methods and sources.
* **How it differs from DIALS, and why:** [`docs/review.md`](docs/review.md).
* **Open issues:** [`docs/outstanding.md`](docs/outstanding.md).

## Building

Requires CMake 3.20 or later, a C++20 compiler, and HDF5's C library; and two
submodules, gemmi (space groups) and bitshuffle (the images' compression). The
configure stops if any is missing.

```sh
# Debian and Ubuntu
sudo apt-get install libhdf5-dev
# or macOS
brew install hdf5
# or: conda install -c conda-forge hdf5

git submodule update --init --recursive
cmake -S . -B build && cmake --build build
ctest --test-dir build
```

Build options:

| option | |
| --- | --- |
| `-DMXI_CUDA=ON` | spot finding and profile fitting on an NVIDIA GPU; `-DCMAKE_CUDA_ARCHITECTURES="80;90"` for cards other than this machine's |
| `-DMXI_METAL=ON -DMETAL_CPP_DIR=DIR` | the same on Apple silicon; `DIR` holds `Metal/Metal.hpp`, Apple's metal-cpp from <https://developer.apple.com/metal/cpp/>. One device backend at a time. |
| `-DMXI_FFTW=ON` | indexing's transform by FFTW, about three times the built-in's speed; GPL, so opt-in. `-DFFTW3_LIBRARY=... -DFFTW3_INCLUDE_DIR=...` where it is not found. |
| `-DMXI_AVX2=OFF` | a binary for x86 CPUs without AVX2 |
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
([`docs/rflx.md`](docs/rflx.md)). Every program also writes what it prints to
`mxi_<program>.log` and takes `--help` and `--version`; the main sequence's
programs take `--timing`, where the time went, and use every core;
`-j N` or `--threads N` sets fewer, where a program has the choice.
`mxi_find --gpu` and `mxi_integrate --gpu` use a GPU where the build has one.

**Several sweeps of one crystal:** name every master to `mxi_import`; each
program takes them all, the sweeps indexed in one basis, each its own crystal
and scale model ([`docs/multi-sweep.md`](docs/multi-sweep.md)).

**Several crystals in one beam:** `mxi_index --max-lattices N`; each lattice
is refined, integrated and scaled as its own experiment, the detector shared,
cosym bringing them to one setting in `mxi_symmetry`
([`docs/multi-lattice.md`](docs/multi-lattice.md)).

**Many crystals, a data set each:** after scaling them together,
`mxi_correlation scaled.rflx --write-clusters` finds the groups that belong
together by their intensities, writes each as its own `.rflx` to scale alone,
and an HTML report ([`docs/correlation.md`](docs/correlation.md)).

**Data sets processed separately** -- each integrated on its own, say -- are
joined by `mxi_combine a.rflx b.rflx ... -o combined.rflx`, then go through
`mxi_symmetry` and `mxi_scale` together; `mxi_split` takes data sets out
again, each alone or a chosen few together.

The options most often wanted:

| program | options |
| --- | --- |
| `mxi_find` | `-g`/`--gpu`; `--min-spot-size N`; `--no-shoeboxes` |
| `mxi_index` | `--d-min D`; `--max-cell A`; `--max-lattices N`; `--no-joint` sweeps indexed apart |
| `mxi_refine` | `--static` one crystal setting for the scan; `--scan-varying N` N control points; `--shared-crystal`; `--joint` several sweeps in one fit, where by default each sweep is refined on its own |
| `mxi_integrate` | `-g`/`--gpu`; `--d-min D`, `--d-max D`; `--postrefine`; `--summation-only` |
| `mxi_symmetry` | `--max-delta D` the lattice's tolerance; `--no-cosym` |
| `mxi_scale` | `--d-min-auto`, `--d-min D`, `--d-max D`; `--space-group NAME`; `--anomalous`; `--absorption-level low|medium|high`; `--error-model-fit regression|individual`; `--scale-restraint W` each sweep's scale tied to its level (1) |
| `mxi_correlation` | `--write-clusters`, `--prefix P`; `--html FILE`; `--min-samples N`; `--dimensions N` |
| `mxi_export` | `--partiality-threshold P`; `--min-isigi S`; `--d-min D` |

### The other tools

| program | |
| --- | --- |
| `mxi_show` | what a file holds: the models, the reflections by what has been done to them, each program's summary ([`docs/output.md`](docs/output.md)) |
| `mxi_report` | an HTML report of a file's processing, every step's figures and plots ([`docs/output.md`](docs/output.md)) |
| `mxi_convert` | between `.rflx` and DIALS' `.expt` and `.refl` (below) |
| `mxi_combine` | several runs' files joined into one: data sets processed separately, to go through symmetry and scaling together |
| `mxi_split` | a file's data sets apart, each into its own file, or the ones named into one |
| `mxi_max` | the largest count on a valid pixel, and whether the images fit in 16 bits (for `mxi_find --gpu-force`) |
| `mxi_profile` | the profile model estimated from strong spots, and how much of each spot it holds |
| `mxi_residuals` | refinement's residuals by resolution, detector module, or anything else |
| `mxi_mask` | shoeboxes marking the integration region, for the image viewer |
| `mxi_grid` | spot density in Kabsch's frame, and spot widths across the detector |
| `mxi_forward` | the model rendered onto the pixels, against the data |
| `mxi_background` | the robust background fitted to pixel values |
| `mxi_readtest` | how fast the images read and decompress |

## Working with DIALS

**The files.** Every program reads DIALS' `.expt` and `.refl` as well as a
`.rflx`, and writes them on request (`--output-expt`, `--output-refl`, or `-o`
naming a `.refl`). `mxi_convert` goes either way: `mxi_convert indexed.rflx`
writes `indexed.expt` and `indexed.refl`; `mxi_convert indexed.expt
indexed.refl` the `.rflx`. Two things DIALS cannot read stay behind: each
sweep's own crystal-less experiment, which mxi carries through indexing
(dials/dials#1029; `--keep-sweeps` keeps them), and mxi's scale model.

**Swapping steps.** Any step can be DIALS': `dials.symmetry` and `dials.scale`
take mxi's integrated reflections, mxi's programs take DIALS' at any stage,
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

Where they differ, and by how much, is in [`docs/review.md`](docs/review.md).

## mxeq: comparisons and regression runs

`mxeq` is a Python package that judges mxi's output, independently of it: it
reads files and imports nothing from `src/`.

```sh
pip install -e python                    # mxeq
pip install -e "python[plots,fixtures]"  # with its plots and test-data makers
```

**Comparing with DIALS** ([`python/README.md`](python/README.md)): `mxeq check` compares two
pipelines at a boundary; `compare-expt` two experiment lists; `trend`,
`explain` and `disagree` where two integrations differ, and why; `unique` two
scaled data sets a reflection at a time; `equivalents` integration's bias
against the data's own symmetry equivalents, with no second program.

**Regression runs** ([`docs/regression.md`](docs/regression.md)): `mxeq regress` builds a commit,
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
skipped ([`python/README.md`](python/README.md)).

## Acknowledgements

mxi implements methods described in the following papers, and is tested
against the software they describe; please cite them alongside mxi:

* **DIALS:** Winter, G., Waterman, D. G., Parkhurst, J. M., Brewster, A. S.,
  Gildea, R. J., Gerstel, M., Fuentes-Montero, L., Vollmar, M.,
  Michels-Clark, T., Young, I. D., Sauter, N. K. & Evans, G. (2018).
  DIALS: implementation and evaluation of a new integration package.
  *Acta Cryst.* D74, 85–97. https://doi.org/10.1107/S2059798317017235
* **XDS:** Kabsch, W. (2010). XDS. *Acta Cryst.* D66, 125–132.
  https://doi.org/10.1107/S0907444909047337; and Kabsch, W. (2010).
  Integration, scaling, space-group assignment and post-refinement.
  *Acta Cryst.* D66, 133–144. https://doi.org/10.1107/S0907444909047374
* **cosym:** Gildea, R. J. & Winter, G. (2018). Determination of Patterson
  group symmetry from sparse multi-crystal data sets in the presence of an
  indexing ambiguity. *Acta Cryst.* D74, 405–410.
  https://doi.org/10.1107/S2059798318002978
* **Scaling:** Beilsten-Edmands, J., Winter, G., Gildea, R., Parkhurst, J.,
  Waterman, D. & Evans, G. (2020). Scaling diffraction data in the DIALS
  software package: algorithms and new approaches for multi-crystal scaling.
  *Acta Cryst.* D76, 385–399. https://doi.org/10.1107/S2059798320003198
* **Clustering:** Thompson, A. J., Beilsten-Edmands, J., Tam, C.,
  Sanchez-Weatherby, J., Sandy, J., Mikolajek, H., Axford, D., Jaho, S.,
  Hough, M. A. & Winter, G. (2025). Enhanced intensity-based clustering of
  isomorphous multi-crystal data sets in the presence of subtle variations.
  *Acta Cryst.* D81, 278–290. https://doi.org/10.1107/S2059798325004589

The other methods mxi uses, and their references, are in
[`docs/algorithms.md`](docs/algorithms.md).

mxi was developed at the Northeastern Collaborative Access Team (NE-CAT)
beamlines, which are funded by the National Institute of General Medical
Sciences from the National Institutes of Health (P30 GM124165). NE-CAT is
at the Advanced Photon Source, a U.S. Department of Energy (DOE) Office of
Science User Facility operated for the DOE Office of Science by Argonne
National Laboratory under Contract No. DE-AC02-06CH11357.
