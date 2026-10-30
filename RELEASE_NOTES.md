# Release notes

Newest first. Each release lists first what changes results or commands -- a
default changed, an option renamed or gone -- then what is new, then what is
fixed. The tag of each (`git tag -n99 v1.0.5`) carries the same summary.

## Unreleased -- since 1.1.2

**Changes commands**
- Every command takes `-j N` or `--threads N`, the same option with the same
  meaning: the threads for everything the command does in parallel, 0 (the
  default) every core the job may use. Four options for one part's threads are
  gone, `-j` covering them: `--fft-threads` and `--jacobian-threads` of
  indexing, `--jacobian-threads` and `--normal-threads` of refinement
  (APS-24ID/mxi#1).

**New**
- `mxi verify`: how much of what each lattice predicts was seen -- every
  reflection predicted to the resolution the indexed spots reach, a prediction
  observed where a strong spot of its sweep has fractional indices within 0.3
  of its own, as indexing judges one, and how many were also indexed to the
  lattice; by resolution shell, predictions on masked pixels counted apart.
  Indexing counts the spots a lattice explains; this counts what it should
  explain and does not, which a supercell or a second crystal leaves high in
  the first and low here. After indexing or any later step
  (docs/multi-lattice.md).

**Fixed**
- "Every core" means every core the job may use: on Linux the process's CPU
  affinity, which a batch system's allocation, `taskset` or a cpuset sets, not
  the machine's count, which on a shared node had every job take every core.
- Refinement's result no longer depends on the number of threads: its normal
  equations are summed in blocks fixed by the data, in block order, the same
  bit for bit on one thread or many. They had differed in the last bits, which
  is why `--normal-threads` existed.
- `mxi correlation`, `mxi show` and `mxi report` say when an option is
  refused; they had ignored an option they did not know whenever the files
  came first.

**Changed**
- `mxi find` shows its progress through a series, as `mxi integrate` does: a
  bar on a terminal, a line each tenth of the way otherwise, from 100 images
  (`--progress-from N`). On a slow machine -- a Raspberry Pi -- a long series
  had run in silence.

## 1.1.2 -- 29 October 2026

**Fixed**
- The Linux releases run on RHEL 9 and its clones (Rocky, Alma): they are
  built in manylinux_2_28, needing glibc 2.28 or later, where 1.1.1's, built
  on Ubuntu 22.04, needed 2.35 and RHEL 9 has 2.34. Packaging now refuses a
  binary needing anything newer.
- mxi builds from source on RHEL 9: eleven sources used `std::sort`,
  `std::unique` and the like without including `<algorithm>`, which newer
  compilers' headers happen to bring in and RHEL 9's GCC 11 does not. A test
  checks every source for the headers it uses.

**New**
- Releases for ARM64 Linux as well: on the CPU, for a Raspberry Pi 4 or 5 with
  a 64-bit OS and ARM servers, and with CUDA for the GH200 and GB200 (from
  NVIDIA's ARM-server toolkit). The CUDA builds use CUDA 12.8, adding
  Blackwell.

**Changed**
- The CUDA builds take the toolkit from NVIDIA's redistributable archives,
  each checked against its SHA-256, rather than NVIDIA's RHEL 8 package
  repository, whose index named files it could not supply and stopped the
  first two candidates of this release.

## 1.1.1 -- 28 October 2026

1.1.0 was tagged but never published: its CUDA build failed to link. 1.1.1 is
1.1.0 with that fixed, and the first release of everything listed under 1.1.0
below -- one program, `mxi`, every step a command of it, and the releases
with FFTW and every licence.

**Fixed**
- The CUDA build links: the CUDA layer the GPU code shares
  (`src/cuda_device.cu`) had been compiled into both `mxi find` and `mxi
  integrate`, and so defined twice once they were one program. It is a
  library of its own now, as Metal's was, and the build refuses at configure
  any source compiled into two commands.
- The documents brought up to date with the single program throughout.

## 1.1.0 -- not published (tagged 28 October 2026)

**Changes commands**
- **One program, `mxi`, every step a command of it:** `mxi import`, `mxi find`,
  `mxi index`, `mxi refine`, `mxi integrate`, `mxi symmetry`, `mxi scale`,
  and the rest, where they had been 24 programs, `mxi_find` and so on. `mxi
  help` lists the commands, `mxi COMMAND --help` each one's options, `mxi
  --version` the version; `mxi view` runs the image viewer. A script that ran
  `mxi_find strong.rflx` runs `mxi find strong.rflx`. The documents and the
  tutorial say so throughout.
- What files record is unchanged, so that files written by 1.0 read as before:
  each step is named in a file's summary as it was (`mxi_refine`), and each
  command's log is `mxi_<command>.log`.
- One program rather than 24, each of which carried its own HDF5: a release
  archive is some 5 MB, where 1.0.7's was some 54.

**New**
- The releases include FFTW for indexing's transform (a vcpkg feature,
  `fftw`). FFTW is GPL, so the released program is distributed under the GPL's
  terms as a whole; mxi's source remains BSD-3-Clause, and a build without
  FFTW -- the default -- uses the built-in transform.

**Fixed**
- The release archives carry the licences of everything linked into the
  program -- gemmi, bitshuffle and its LZ4, HDF5, libaec, zlib, FFTW, and
  NVIDIA's CUDA runtime in the CUDA build -- in `licenses/`, as those licences
  ask of binary copies, with a README saying which is which.

## 1.0.7 -- 28 October 2026

**New**
- Releases on GitHub, built when a tag is pushed (or by hand, `gh workflow run
  release.yml -f tag=...`): Linux x86_64 on the CPU and with CUDA, and macOS
  arm64 with Metal (macOS 13 and later), each self-contained -- HDF5 linked in
  -- and tested before it is published; the complete source with its
  submodules; mxeq's wheel (`docs/release.md`). The README says how to install
  them.

**Fixed**
- Building with `-DMXI_TESTS=OFF` makes every program: `mxi_import`,
  `mxi_integrate`, `mxi_max` and `mxi_readtest` had been defined only with the
  tests.
- `mxi_view` marks a spot on exactly the images its box spans. Its extent had
  come from the centroid's variance -- its uncertainty, not the spot's spread
  -- with an image either side at least, so that a spot on one image was
  marked on its neighbours too: against the spot finder's boxes on a 16M
  sweep, only 15% of spots had been marked on just their own images.
- `mxi_view` starts in about a second on an integrated 16M sweep of millions
  of reflections, where it took some 15 s: its markers are picked out for the
  image shown as the slider moves -- a few hundred -- rather than built for
  every image of the sweep at once, 6 million points for the spots and 26
  million for the predictions.
- `mxi_report`'s cell over the scan is drawn against the rotation angle, a
  point a degree over the whole scan. It had read the scan-varying model's
  control points -- a few dozen -- as image numbers, so that a 3600-image
  sweep's cell was drawn over its first 37 images; `mxi_refine`'s printed
  ranges came from the same control points and are now the model's values.
  Reports of files refined before still draw correctly, from the file's
  crystals.

## 1.0.6 -- 26 October 2026

**Changes results**
- `mxi_integrate --min-zeta` now also keeps reflections below the cut out of
  the profile model's estimate, which had used a fixed 0.05.

**Changes commands**
- The diagnostic that marked each reflection's foreground for
  dials.image_viewer, `mxi_mask`, is now `mxi_foreground`; there is no program
  called `mxi_mask` any more.

**New**
- `mxi_view`: the images in napari (an optional dependency), a sweep at a
  time, each shown summed over as many images as asked, with the strong spots
  and predictions on every image they are recorded on; and the detector's
  untrusted regions -- a backstop, its arm, a shadow -- drawn as rectangles,
  ellipses and polygons, with Undo, Clear and Same as last, and saved to
  `masked.rflx`, the input left as it was. `mxi_untrusted` lists, sets and
  clears them with no window. Spot finding and integration leave them out as
  they do the detector's bad pixels (`docs/mask.md`).
- `mxi_scale --exclude E:FIRST:LAST`: an experiment's images left out of
  scaling, as `dials.scale exclude_images=`; as many ranges as wanted,
  numbered as the file's experiments and images.
- `mxi_refine --sample-offset`: each lattice's centre of diffraction off the
  rotation axis, turning with the crystal -- two parameters, amplitude and
  phase -- used by prediction and so by integration.
- `mxeq truth`: a scaled data set against a refined structure's Fc^2 (a SHELXL
  `.fcf`): calibration -- against Fc^2, and against each reflection's own
  equivalents -- bins by what processing could get wrong, the model's part
  against the processing's, figures, the worst observations' shoeboxes.
- `cmake --install` puts every program, `mxi_view` with them, in one `bin`.

**Changed**
- `mxi_integrate`'s progress is a bar, as tqdm draws one, on a terminal, as
  wide as the terminal allows -- written to the terminal alone, so the log does
  not keep it -- and the plain line each tenth of the way, as before, where the
  output is not a terminal.

## 1.0.5 -- 22 October 2026

**Changes results**
- `mxi_symmetry` puts an orthorhombic cell's axes in order of length,
  a ≤ b ≤ c, as DIALS reports (and any other non-monoclinic lattice's, where
  its symmetry leaves a choice). Before, every permutation tied and the first
  found won, so a cell already in order could be reindexed `c,a,b`. Cell
  lengths within 0.5% count as equal, so a cubic cell's operator is
  unchanged.

**Fixed**
- Regression tests for the above: the setting rule, and an orthorhombic
  crystal through the whole symmetry analysis.

**Documents**
- `docs/outstanding.md` checked against the code and its history: four items
  closed, nine rewritten to what is now true.

## 1.0.4 -- 20 October 2026

**Changes results**
- `mxi_correlation` normalises each data set by its own Wilson fall-off before
  correlating them, as DIALS does (`--no-normalise` for the old behaviour). A
  crystal whose intensities fall off differently had looked less like the rest
  than it is, and could be cut off as an outlier.

**New**
- `mxi_index --known-cell a,b,c,α,β,γ` with `--known-space-group NAME`: the
  basis chosen to match a known cell, as `dials.index unit_cell=`. Refused
  where nothing matches; a warning where the match indexes under half what the
  data's own best basis does.
- `mxi_index --real-space-grid`: the basis vectors by the real-space grid
  search of Gildea et al. (2014), given a known cell (`--grid-step DEG`).
- A tutorial, [`docs/tutorials/cows-pigs-people.md`](docs/tutorials/cows-pigs-people.md):
  36 insulin crystals from three species, processed together and sorted into
  species by `mxi_correlation`, with a script that runs it.

**Fixed**
- `mxi_integrate` on several sweeps no longer prints the names of its internal
  scratch files.

## 1.0.3 -- 18 October 2026

**Changes results**
- `mxi_refine` refines each sweep on its own by default, even several sweeps of
  one crystal that share a beam and detector (`--joint` for one fit of all).
- The scan-varying pass starts from the static pass's reflections, its outliers
  kept out. It had re-admitted them -- on a small molecule a third of the
  reflections -- and could end worse than the static fit, the cell wandering
  within a sweep; RMSDs are now DIALS's.
- Each scan-varying control point is tied weakly to the static crystal
  (`--crystal-restraint W`, 10⁶ by default; 0 for none).

**New**
- `mxi_refine --interval DEG`: a scan-varying control point every DEG degrees of
  each sweep (10 by default; 20–30 suit sparse data).
- A warning when a crystal changes within a sweep by more than a crystal can.

## 1.0.2 -- 18 October 2026

**Changes results**
- `mxi_index --no-joint` keeps each sweep's own beam and detector. It had
  merged identical beams into one, which `mxi_refine` then fitted for all of
  them: data sets processed together scaled worse than processed apart. Re-index
  such runs.
- `mxi_scale` ties each sweep's scale weakly to its own level
  (`--scale-restraint W`, 1 by default), where an unrestrained short wedge had
  wandered through zero; a scale at or below zero is a warning.

**New**
- `mxi_combine`: several runs' files joined into one, to be taken through
  symmetry and scaling together. `mxi_split`: a file's data sets apart, each
  alone or a chosen few together.

**Fixed**
- A file holding several data sets -- from `mxi_split`, `mxi_combine` or
  `mxi_correlation --write-clusters` -- integrated every data set but the first
  from the first data set's images.
- `--help` was an error in seven programs.
- The documents' commands are checked against the programs' own options by a
  test; obsolete running instructions corrected.

## 1.0.1 -- 17 October 2026

**Changes commands**
- `mxi_find` takes its input as every program does: an experiment list or an
  NXmx master, positionally. `-e` and `-gpu` are gone.
- Build options are `MXI_*`: `MXI_CUDA`, `MXI_METAL`, `MXI_AVX2`, ... The old
  `SPOTFINDER_*` names are accepted, with a warning, through 1.0.x.
  `MXI_SANITIZE` takes `ON` or a list: address, undefined, thread.
- `mxeq.checks.spot_table` and `mxeq.checks.planted_spots` (were
  `spotfinder_refl` and `spotfinder_spots`).

**Changed**
- The spot finder is one stage of mxi: one source tree and one build, the
  library's readers and writers, one Metal and one CUDA layer shared with
  profile fitting, its tests in `mxi_tests`. Strong spots are byte for byte the
  same.
- `.refl` files list their columns in DIALS's order, shoeboxes included.

**Documents**
- The README says what mxi is and what NXmx is, gives the acknowledgements, and
  its links work on GitHub (and are tested to keep working).

## 1.0.0 -- 16 October 2026

The first release. MX data processing for NXmx rotation data, from images to
scaled, unmerged intensities: `mxi_import`, `mxi_find`, `mxi_index`,
`mxi_refine`, `mxi_integrate`, `mxi_symmetry`, `mxi_scale` and `mxi_export`,
with `mxi_correlation` for many crystals and tools around them. A new,
standalone package written in C++, compatible with DIALS: it follows DIALS's
methods, reads and writes its files, and any step can be exchanged for
DIALS's own. Spot finding and profile fitting on a GPU, through CUDA or Metal.
`mxeq`, in Python, compares mxi with DIALS and runs regression tests over real
data. The build requires HDF5 and the bitshuffle submodule; every program
answers `--version`, and takes `-j` or `--threads`.
