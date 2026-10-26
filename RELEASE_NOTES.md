# Release notes

Newest first. Each release lists first what changes results or commands -- a
default changed, an option renamed or gone -- then what is new, then what is
fixed. The tag of each (`git tag -n99 v1.0.5`) carries the same summary.

## Unreleased -- since 1.0.5

**New**
- `mxi_mask`: draw the detector's untrusted regions -- a backstop, its arm, a
  shadow -- over a sum of each sweep's images, as circles, rectangles and
  polygons; writes `masked.rflx`, the input left as it was. `mxi_untrusted`
  lists, sets and clears them with no window. Spot finding and integration
  leave them out as they do the detector's bad pixels (`docs/mask.md`).
- `mxi_refine --sample-offset`: each lattice's centre of diffraction off the
  rotation axis, turning with the crystal -- two parameters, amplitude and
  phase -- used by prediction and so by integration.
- `mxeq truth`: a scaled data set against a refined structure's Fc^2 (a SHELXL
  `.fcf`): calibration, bins by what processing could get wrong, the model's
  part against the processing's, figures, the worst observations' shoeboxes.
- `cmake --install` puts every program, `mxi_mask` with them, in one `bin`.
- `mxi_scale --exclude E:FIRST:LAST`: an experiment's images left out of
  scaling, as `dials.scale exclude_images=`; as many ranges as wanted,
  numbered as the file's experiments and images.

**Changes results**
- `mxi_integrate --min-zeta` now also keeps reflections below the cut out of
  the profile model's estimate, which had used a fixed 0.05.

**Changes commands**
- The diagnostic that marked each reflection's foreground for
  dials.image_viewer, `mxi_mask`, is now `mxi_foreground`: the name went to the
  mask editor.

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
