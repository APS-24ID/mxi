# mxi -- working notes

## What this is

MX data reduction for NXmx rotation data, from the images to scaled, unmerged
intensities: spot finding, indexing, refinement, integration, symmetry and
scaling, for one sweep or many, one crystal or many -- with `mxi_correlation`
to find which crystals belong together, and `mxi_combine` and `mxi_split` to
join and part runs. It reads and writes DIALS' formats, so any stage can be
swapped for DIALS' own. Version 1.0.5; `RELEASE_NOTES.md` says what each
release changed. `docs/review.md` is the guide for someone meeting the code
for the first time; `docs/algorithms.md` lists the methods.

Two pieces share this repository and are not one program:

| | where | its notes |
| --- | --- | --- |
| the pipeline, spot finding among its stages | `src/`, `apps/` | this file; `docs/spots.md`, `docs/integration.md`, `docs/symmetry.md`, `docs/scaling.md`, `docs/correlation.md` |
| `mxeq`, the referee | `python/` | `python/CLAUDE.md` |

`mxeq` must stay independent of everything else here. It is the referee, and a
referee that depends on the thing it judges is not one. Nothing in `src/` or
`apps/` may import it, and it must never import them; it reads files, which is
the whole point.

The build requires HDF5 and the gemmi and bitshuffle submodules, and stops at
configure without them (`git submodule update --init --recursive`): nothing is
left out silently. The spot finder was a project of its own until October; it
now uses the library's table writer (`src/spot_table.cc` builds a library
table) and the library's experiment lists (`src/experiment_images.cc`, the one
reader tolerant of stills, kept deliberately).

## The documents, and what is open

`docs/README.md` says what each document is for -- reference, notes or history
-- and `docs/outstanding.md` is the one list of open work. When something on
that list is done, it leaves the list, and the document it concerns says what
was found.

## Handover: where the work stands

STATUS: written 4 October 2026, 325 commits; several sweeps, the cell
regularised and mxi_export added 6 October, 349 commits; several crystals,
sweeps carried through and the variances 12 October, 414 commits; the
programs' output, mxi_show and mxi_report, compression and the fit's weights
13 October, 434 commits; the error model, the regression runs, the cell
from two-theta and memory and time on long scans 14 October, 473 commits;
mxi_correlation, partials, the error model's fits, centred input and
integration's speed 15 October, 494 commits.
1.0.0 to 1.0.5, 16 to 22 October, 530 commits: the spot finder made one
stage of mxi; joint and isolated processing made to agree (--no-joint's
beams, refinement sweep by sweep); refinement's scan-varying pass and its
restraints, from a small molecule; mxi_combine and mxi_split; --known-cell and
the real-space grid search; mxi_correlation's normalisation; the
cows-pigs-people tutorial; the orthorhombic setting; release notes.
22 and 23 October, to 547 commits, unreleased: `mxeq truth`, a scaled data set
judged against a refined structure's Fc^2 (figures, zeta, the outliers'
shoeboxes); `mxi_integrate --min-zeta` reaching the profile model; untrusted
regions -- `mxi_view` to draw them, `mxi_untrusted`, applied in spot finding and
integration; `cmake --install`; `mxi_refine --sample-offset`, each lattice's
centre of diffraction off the axis.
For a new session to
pick up
from. `docs/outstanding.md` is the list of open work; this section is what is in
flight and how the work is done here. People follow this project: keep the
README, `docs/review.md` and this handover true as things change.

**The project is mxi.** It began as dials-metal, and its history was rewritten
with git filter-repo to the name every program already had; the old history is
kept apart. The chain is complete for one sweep, for several of one crystal
indexed together (`docs/multi-sweep.md`), and for several crystals in one beam
(`docs/multi-lattice.md`), and needs nothing of DIALS:
`mxi_import`, `mxi_find`, `mxi_index`, `mxi_refine`, `mxi_integrate`,
`mxi_symmetry`, `mxi_scale`, each interchangeable with DIALS at its boundary,
every core by default.

**mxi_export** (6 October, `docs/export.md`): an unmerged MTZ as dials.export
writes it -- the same intensities and corrections, columns, batch numbering
and headers in the Cambridge frame -- through gemmi's MTZ writer. AIMLESS reads
its output and understands it (Graeme); no column-by-column comparison with
dials.export has been made.

**Several sweeps** (5 and 6 October, `docs/multi-sweep.md`): every program from
`mxi_import` to `mxi_scale` takes several sweeps of one crystal, as DIALS does,
and works well -- Graeme -- on the four cubic insulin sweeps at their own
orientations of https://zenodo.org/records/8376818 (ins10_1 to ins10_4; ins10_1 is
the insulin sweep used throughout).
The protocol, Graeme's from DIALS: index with one matrix, which puts the sweeps
in one basis; the last indexing cycle and `mxi_refine` then refine a crystal for
each sweep apart (`--shared-crystal` keeps one); integration takes each sweep
alone and joins them; symmetry pools them; scaling gives each a model of its
own, the overall scale and B fixed over all of them together. On the way: a
scan not starting at image one had every angle wrong by its start (z anchored
at image one, where dxtbx anchors it at the scan's first image); the trusted
range is the lower of the detector's count limit and 2^bits - 3, a broken link
into a `_meta.h5` said so; and two sessions once ran at once on this
repository -- check `git log` and `git status` before starting work.

**Done since the last handover (308):** `mxi_max` and `mxi_find --gpu-force`,
32-bit frames thresholded on a GPU of 16 bits; both of the detector's markers --
0xfffe a bad pixel, 0xffff a tile join -- now invalid in integration, where a
bad pixel had been a count (item 50); `--gpu` saying why no device is usable, in
CUDA's words (a GH200 failed for a CUDA 13.4 toolkit on a 12.9 driver); each
reflection's reason for not being profile fitted, `profile.failure`, and
`mxeq failures` (item 51); `mxi_find --subtract-background`, the profile model
subtracting what a shoebox records; `mxeq unique`, two scaled data sets one
reflection at a time, and `mxeq observations`, every observation of one;
`--d-max` in integration and scaling; and the background's dispersion, under
`mxi_integrate --save-background-parameters`.

**Done since 349 (to 414), 7 to 12 October:**

* **`.rflx`** (`docs/rflx.md`): one HDF5 file a step in dxtbx-h5's layout, every
  program's own; the pair as well; `mxi_convert` between them.
* **Several crystals** (`docs/multi-lattice.md`, Graeme's two crystals of
  insulin, ins2_27): `mxi_index --max-lattices`, every spot assigned again
  globally -- the first lattice had taken 1700 of the second's by chance --
  ambiguous spots left out of refinement; models shared as the list numbers
  them, refined once; a sweep's lattices integrated in one pass, a profile
  model and reference profiles each, `--gpu` a batch each; cosym in
  `mxi_symmetry`, the same coset as dials.cosym; a model each in scaling;
  batches in blocks of the next power of ten (1 to 600, 1001 to 1600).
  Against DIALS: the same lattices and indices for 20762 of 20763 spots,
  intensities CC 0.998, the scaled ones 0.9997.
* **Sweeps carried through** (dials/dials#1029), the default: each sweep's
  crystal-less experiment first, its unindexed spots its own, through every
  program (`src/views.cc`, the two views in one place); `mxi_convert` leaves
  them out for DIALS by default (`--keep-sweeps`). A lattice missed is
  recovered beside refined ones.
* **Neighbours masked** in integration, as DIALS's 3d_threaded and not its
  default; the search keyed on images too after it took minutes on ten turns.
* **The variances** (`docs/integration.md`, "The fitted variance";
  `docs/scaling.md`, "Variances from the model"): the fitted variance's
  background term through the fit's weights; a pixel's least variance 0.001
  for the fit's variance, 0.05 for its weights (item 61);
  summation's from expected counts; each variance with no signal written, and
  scaling's built on the model's intensity. Graeme's ten-turn insulin: outer
  shell merged I/sigma -5.7 to +3.8. Confirmed on Metal (Graeme).
* **`--grid-points 6`** by default (item 2); import takes a panel's size from
  the data, which some masters' metadata has swapped.

**Done since 414, 12 and 13 October:** what the programs say
(`docs/output.md`). Each .rflx carries the processing's summary, `/summary`,
each program's figures and tables after its input's (`src/processing.cc`);
`mxi_show` tells the state of a file; developer detail under `--verbose`;
indexing's candidates, the cell's e.s.d.s (now scaling's, from two-theta), the
scan-varying cell's range, symmetry's elements by order and axis, the scale by
image, the MTZ's summary, spots by resolution and ice rings. Mind the names:
`src/summary.*` is integration's summary tables -- the processing's summary
was once written over them, unchecked, and restored.

## Two lessons the whole repository rests on (from the README, October 2026)

**A round-trip test cannot detect a wrong convention.** A file written and read
by the same code agrees with itself whatever it means. Five separate format
bugs here were invisible until DIALS read the output: the `.refl` container
shape, the geometry conventions, `0` written where `0.0` was required, blocks
dropped on write, and flags never set. Everything that matters is checked
against a file somebody else wrote.

**When what you are measuring is at the sampling scale, render the model into a
measurement rather than comparing a measurement to a model.** Spots here are
two pixels across. Comparing their moments with a Gaussian in closed form
suggested the model was wrong by factors and that the spots were anisotropic;
integrating the same model over the same pixels inside the same mask put the
positions within a hundredth of a pixel, the widths within four per cent, and
predicted the anisotropy. Most of the disagreement had been the comparison.

**Done 16 to 22 October, 1.0.0 to 1.0.5:**
- The spot finder made one stage of mxi, in three commits (1.0.1): the flat
  tree and one build, MXI_* options; the library's table writers, experiment
  lists, arguments, HDF5 handles and namespaces; one Metal and one CUDA layer
  shared with profile fitting; its tests in mxi_tests; its documents mxi's.
  Strong spots byte for byte throughout.
- mxi_combine and mxi_split, shoeboxes kept (src/combine.*); every program's
  --help works (seven had been errors); the documents' commands checked
  against the programs' options by a test.
- mxi_scale: each sweep's scale tied weakly to its own level (W 1), a scale at
  or below zero a warning (item 71).
- Item 72: --no-joint had merged identical beams, and mxi_refine fitted one
  beam for 14 crystals; now each sweep keeps its own and independent
  experiments are refined apart. And a list's images come from its
  experiments' own image set, not the document's first.
- 18 October: mxi_refine refines each sweep on its own by default, even one
  crystal's sweeps sharing a beam and detector (Graeme: better results
  throughout); --joint for one fit until a test case needs it.
- 18 October: refinement's scan-varying pass starts from the static pass's
  reflections (it had re-admitted its outliers -- on a small molecule, a third
  of them, RMSDs 2-5x DIALS's and a wandering cell); the crystal tied weakly
  to the static one (--crystal-restraint, 1e6); --interval DEG; a WARNING when
  a crystal changes more than a crystal can. Insulin identical to the digit.
  Lesson: a tie that only bites at weights nobody would call weak, with a
  response that is not monotonic, is a symptom -- look for what is dragging
  the fit before tuning the tie.
- 19 October: mxi_index --known-cell (item 74) and --real-space-grid (Gildea
  et al. 2014): the known cell's lengths leaving only directions, the grid
  search's candidates through the same known-cell choice; identical to the
  FFT's basis on insulin and the small molecule (item 75: its advantage still
  to be measured). The cows-pigs-people tutorial drafted, its output to come
  from Graeme's run.
- 20 October: mxi_correlation normalises each data set by its Wilson
  fall-off before correlating (as DIALS): a weak cow, CIX10, had fallen out as
  an outlier. On the way, found: mxi's sigmas for weak reflections calibrate
  at about 1, DIALS's at 0.16 (2.5 times too large). The tutorial's output
  from Graeme's run; each sweep's integration no longer names its scratch
  files; docs/outstanding.md reviewed against the code (four closed, nine
  rewritten).
- 22 October: mxi_symmetry's settings in order of length, a <= b <= c, for
  any lattice but monoclinic -- an orthorhombic cell in order had been
  reindexed c,a,b -- with a regression test through the whole analysis;
  RELEASE_NOTES.md.
- **A release**, each point release the same: the version in five places
  (CMakeLists.txt's project(), vcpkg.json, python/pyproject.toml,
  python/src/mxeq/__init__.py, docs/spots.md's example output -- the configure
  checks the first two agree); an entry at the top of RELEASE_NOTES.md, what
  changes results or commands first; an annotated tag, `git tag -a vX.Y.Z`,
  with the same summary.
- 22 and 23 October, from Graeme's small molecule refined to R1 0.026: `mxeq
  truth` (python/src/mxeq/truth.py) matches each scaled observation to Fc^2 by
  its unique index and position, and found the outliers to be shadows (the
  backstop's arm on l-cyst_04, its edge elsewhere), reflections near the axis,
  and in l-cyst_03 images 600-800 every spot displaced together (item 78).
  Answers: `--min-zeta` (0.2 for several sweeps; the profile model now follows
  it); untrusted regions (src/untrusted.cc, the image set's "untrusted" list,
  merged into the pixel mask by Series::add_masked; drawn at first in a
  matplotlib `mxi_mask`, since replaced by `mxi_view` -- napari, optional,
  tools/mxi_view.py and the shared tools/mxi_tools.py -- writing masked.rflx
  through `mxi_untrusted`); the sample
  offset (Goniometer::sample_offset, d(phi) = S R(phi) c, two parameters a
  lattice, exact derivatives -- every other derivative carries it). The old
  diagnostic `mxi_mask` is now `mxi_foreground`.
- Lessons of 22 and 23 October: check a file's conventions before reading its
  numbers -- intensity.scale.value is before the inverse scale, and the flags'
  bits are mxi's, both guessed wrong at first and both made artefacts look like
  findings; a tool run from a scratch directory must make its paths absolute;
  a name already taken in the build is overwritten silently -- check the
  build's targets before naming a program; a new block of parameters must
  reach every place the old ones are gathered (the normal equations' spans),
  and a recovery test is what shows it does.
- Lessons of 20 to 22 October: run the reference program on mxi's own output
  before theorising -- dials.correlation_matrix on mxi's scaled cows settled
  in one run what weighting and attenuation theories had not; and the
  bundle is the repository: the working copy, the build and everything in
  /tmp can be gone at the start of a session (restore with `git clone
  /mnt/user-data/outputs/mxi.bundle mxi`, then the submodules, then build).
- Lessons: a tie must be to the right thing -- tied to one, the scale pulled a
  weak crystal's level, tied to its own level it does not; test with what
  differs in the user's case (different files, different geometry), not only
  what is to hand (one file split) -- three of today's bugs were invisible to
  one-file tests; and a print of a list's models that leaves one out (the
  beam) can hide exactly what is wrong.

**Done 15 October, to 494:**
- mxi_correlation, whole (`docs/correlation.md`): the sigma-weighted matrix as
  DIALS's weighted_cchalf, cosym coordinates and their dimension as
  dials.correlation_matrix, OPTICS ported from scikit-learn, the clusters as
  files and an HTML report with Ward dendrograms; the same clusters, dimension
  and min_samples as DIALS on Graeme's 36 insulin crystals.
- mxi_scale: partiality divides the summed intensity only; split reflections
  combined first (sum_partial_reflections); `--error-model-fit individual`,
  dials.scale's real default (regression had been called it, wrongly); the
  split's refusal said (`docs/scaling.md`).
- mxi_symmetry: a centred input taken to its primitive cell first (item 70);
  cosym's correlations by hash lookup across threads, its eigen-decomposition
  by Householder and QL (83 s to 1.4).
- mxi_integrate: a frame reader a worker for the run's images, not a
  thread_local -- the calling thread had kept sweep 1's reader for every
  sweep; workers' exceptions reported; the background fitted on a histogram
  of its counts; centroids over the foreground alone. Ferritin 60.7 s to 44.
- mxi_max reads by the series' keys and takes an .rflx; `--sigma-m-estimator`.
- Lessons: read DIALS's source for its defaults, never recall them -- two
  wrong today (sigma_m's scan limit, the error model's fit); measure on the
  machine that matters -- this one has one core and showed no contention; and
  a lesson written here (no thread_local state in parallel calls) is worth
  nothing if the code is not checked against it.

**Done 14 October, to 473:**
- The cell's e.s.d.s from two-theta in mxi_scale, as dials.two_theta_refine;
  refinement's covariance e.s.d.s removed (`docs/scaling.md`).
- The error model fitted by regression across intensity bins, then split, one
  a for the variance with no signal and one for the rest; the sigmas'
  calibration printed and in the report; the anomalous slope's control,
  random halves (`docs/scaling.md`). 6S0R parked (item 62).
- `mxeq regress`: runs over real data sets, history, compare, show, clean,
  tiers; Graeme runs it on his machines (`docs/regression.md`).
- Speed and memory on 36000 images: refinement's Jacobian sparse (82.7 GB to
  5.1); one table held, not three, in refinement, symmetry and scaling; the
  scaling normal equations a group at a time, the covariance from one
  factorisation (scaling 188 s to 22); profiles learned into sparse partial
  sets; `read_chunk`, H5Dread_chunk by pread outside HDF5's lock, one
  descriptor a file; a table of each panel's pixels in mm (`docs/integration.md`).
  `--timing` gives each phase's peak memory.
- mxi_symmetry: cosym only for crystals whose orientations differ; monoclinic
  cells in the convention's order. mxi_index `--no-joint`; mxi_import
  `--split N`. mxi_refine refines the beam, analytically, by default;
  mxi_convert leaves crystal-less sweeps out of DIALS files by default.
- The scale model written into scaled .rflx (mxi_physical, item 26); the
  report a point a degree, the model drawn.
- Lesson: on 6S0R and the 36000 images, measuring each phase (`--timing`,
  calibration tables) found the cause each time; guessing before measuring
  -- the foregrounds reused in opening -- made ferritin slower, and was undone.

**Done 13 October, to 434:** mxi_report (`docs/output.md`), a plot a sweep,
refinement's own cell over the scan; the cell's e.s.d.s from one
factorisation, scan-static only -- the first version inverted the normal
matrix a column at a time and made Graeme's four small-molecule sweeps very
slow; the fit's weights apart from its variance (item 61, closed: the 0.001
floor had wrecked strong reflections on a near-empty background --
R merge 0.036 to 0.069 -- found only on Graeme's data); `mxi_scale
--no-model-variances`; .rflx numbers bitshuffle-LZ4 (`docs/rflx.md`), the spot
finder's file 14 MB to 2.4; item 4 measured again -- the background's bias a
constant fraction, not curvature. Graeme plans a regression system over full
data sets for every change set.

**The cell's e.s.d.s are mxi_scale's, from two-theta** (Graeme, 13 October;
`docs/scaling.md`): refinement's covariance e.s.d.s are removed. Each centroid's
variance is floored at a pixel's, size^2 / 12 -- without it strong low-resolution
spots, their counting variances absurdly small, pulled insulin's a 0.02 A off.

**Paused: the backstop** (`docs/backstop.md`, "Where it stands"; items 52, 53).
Against DIALS on ferritin the innermost shell agreed worst; the cause, the
reflections beside the backstop -- its shadow attenuating signal and background
together, its flare raising the background -- and scaling's outlier rejection
keeping the attenuated over the unharmed. Tools built to see, model and filter
it (`mxeq unique`, `observations`, `background`, `background-model`,
`attenuation`; `--d-max`), none a default. On ferritin the filters raised Rmerge
slightly, and DIALS's and mxi's data refined about the same before any of it.
Paused for want of an absolute measure of processing quality -- simulated data
with the truth known, anomalous peak heights at known scatterers, or agreement
with an independent reference model are the candidates. Terms: the backstop
shadow is the low-count region, the flare the high-background region around it.

**Waiting on Graeme's data or hardware:**

* **The CUDA side** of the October fit changes (the background term, the
  floor, the variance with no signal, `FitF`'s new field): Metal built and
  ran them (Graeme); CUDA not built here.
* **The anomalous slope** on the ten-turn insulin rose from 1.04 to 1.52 with
  the model's variances; on the data here Friedel and random-half splits
  agree, so likely real signal. A random-half line beside the slope in
  `mxi_scale` would settle it on his data -- offered, not built.

* **8BXT, a third not profile fitted** (item 51): `mxeq failures` on its table,
  integrated without `--gpu`, says which reason and where.
* **`mxi_import` across data sets** (item 49): dials.import and mxi_import on
  each master, `mxeq compare-expt`.
* **The fused CUDA threshold** (`SPOTFINDER_GPU_FUSED=1`): compiled, never run.
* **The CUDA and Metal sides** of the device-reason messages: not compiled here.

**Measured and settled:** CUDA profile fitting on an RTX 4060, 94.7 s against
103.8 on dense ferritin; HDF5's lock not what limits reading (item 48); on
ferritin sigma_m agrees with its partials; and insulin's is not too wide
either -- at 0.086 rather than 0.128 degrees it merged worse (items 1, 3).

**Next (Graeme):** item 4, the background half to one per cent high at low
resolution on bare background -- a plane or a radial term to test the
curvature near the beam.

**Paused:** background and summation on a device (item 47).

**Removed, and why:** reading frames around HDF5's lock with pread, and not
zeroing shoeboxes on opening: each byte-identical, neither faster where it
mattered. A model of the images' background, for the backstop: a heuristic,
Graeme's objection, and never built.

### How work is done here

* **What a program prints describes the data it was given.** Results from the
  data this was developed on -- insulin's sigma_m, a MacBook's timings --
  belong in the documents and the commit messages, never in output, where they
  read as though they were about whatever is being processed: `mxi_profile`
  printed insulin's DIALS values on every run. Comments and docstrings may cite
  what was measured. `python/tests/test_output_names_no_data.py` holds the
  programs to it.

* **Commits** are authored and committed as Graeme Winter
  <graeme.winter@gmail.com>, with `Co-Authored-By: Claude <noreply@anthropic.com>`,
  and a message that says what was found and how it was checked.
* **Before a commit**: `clang-format --style=file -i` on every changed `.cc`,
  `.hh` and `.cu`, and none left out of format; GCC and Clang builds with no
  warnings; `mxi_tests` under both, `ctest`, and the Python suite; and every
  output a change should not touch compared byte for byte, before and after.
* **The Python suite** wants the programs and data named:
  `MXI_FIND=<build>/mxi_find MXI_TEMPLATE_EXPT=/mnt/user-data/uploads/refined.expt
  MXI_INTEGRATE=<build>/mxi_integrate MXI_TEST_EXPT=/mnt/user-data/uploads/integrated.expt
  MXI_TEST_REFL=/mnt/user-data/uploads/refined.refl
  MXI_TEST_IMAGES=/mnt/user-data/uploads/ins10_1.nxs PYTHONPATH=src:tests
  python3 -m pytest -q`, run in `python/`.
  `python/tests/test_documents.py` fails on any path or program a document names that
  does not exist. `python/tests/test_lattices.py` wants the two crystals of
  ins2_27 at each step: `MXI_TWO_LATTICES` (strong spots), `..._REFINED`,
  `..._INTEGRATED`, `..._SYMMETRIZED`, `..._SCALED`, with the programs named.
* **Delivery** is a bundle: `git bundle create <file> --all`, `git bundle
  verify`, copied to `/mnt/user-data/outputs/mxi.bundle`; then a fresh
  clone, `git submodule update --init`, `cmake -DMXI_FFTW=ON`, a build and every
  suite from it, before `present_files`.
* **The reference chain**, 300 images of insulin, from the uploads, one
  `.rflx` a step (`docs/rflx.md`; `--output-expt` and `--output-refl` for the
  DIALS pair instead):
  `mxi_find /mnt/user-data/uploads/imported.expt -j 4` (strong.rflx: the
  models with the spots);
  `mxi_index strong.rflx`;
  `mxi_refine indexed.rflx --analytic --scan-varying 5` (on 30
  degrees the default, scan-varying at one point per 10 degrees and at least
  five, gives the same bytes);
  `mxi_integrate refined.rflx --threads 4`;
  `mxi_symmetry integrated.rflx`;
  `mxi_scale symmetrized.rflx`. It gives I 2 3 by
  b+c,a+c,a+b, Rmeas 0.041, error model a 1.009 and b 0.0244 (1.017 and 0.0236
  with five scan blocks, the default before 29 September 2026).
* **A CUDA compile check** is possible here, not a run: `apt-get install
  nvidia-cuda-toolkit` (12.0, some 5 GB), then `cmake -DSPOTFINDER_CUDA=ON
  -DCMAKE_CUDA_ARCHITECTURES=89`. A GPU branch of `find_spots.cc` alone can be
  compiled syntax only with `-fsyntax-only -DSPOTFINDER_GPU -DSPOTFINDER_CUDA
  '-DSPOTFINDER_VERSION="check"'` and the include paths from the build's
  `flags.make`.

### Traps in this environment

* **An edit that matches exact text fails once clang-format has reflowed it.**
  Scripted edits here assert that their text occurs once before writing
  anything, so a failed one writes nothing -- and a build that then passes has
  built the unchanged file. Match by pattern, or by the braces of a block, and
  check the edit landed.
* **Commands run under `/bin/sh`**, which has no `<(...)` and no `${var/...}`.
* **A tool call stops at about five minutes.** Build in the background --
  `(build; echo EXIT $? >> log) & sleep 240` -- and run each suite in a call of
  its own.
* **One core**: determinism across thread counts can be tested here, speed-ups
  cannot.
* **The disk**: some 8 GB free at the start, the CUDA toolkit takes 5, and when
  the disk filled every command failed, even `df`, until space was freed.
  Remove build trees and sanitizer builds that are done with.
* `mxeq.refl.load(path)` reads a reflection table in Python.

## Hard constraints

**No third-party dependencies.** Same rule as the spot finder's hand-rolled
msgpack writer: this builds on a laptop and on a beamline machine with nothing
installed. Vec3, Mat3, the Cholesky solver and the test harness are all here
for that reason, and each is under a hundred lines.

**Thresholds come from measurement, not from comfort.** The real-data tests
assert against numbers that were printed first: parallax 0.96 px over the
sample, residual median 1.9e-4 and max 7.1e-4, agreement with DIALS' `rlp` at
4.3e-5. Median and maximum are asserted separately, because the residual
distribution has a tail and bounding only the maximum bounds only the tail.

**A fixture that cannot diffract is a silent test, not a failing one.** The
synthetic panel was at `+z`, consistent only with the wrong sign of s0.
Correcting the sign made it predict nothing, and the prediction tests passed
over an empty list. They now assert counts.

**Double, not float, throughout the geometry.** Positions are wanted to a
millipixel over a 4000-pixel detector, one part in 4e6, and float has seven
digits. Apple GPUs have no doubles at all, so any Metal port of prediction has
to demonstrate the precision it achieves against this reference rather than
assume it — the summation-order differences already seen in the spot finder's
centroids are the floor, not the ceiling.

**Build with Clang as well as GCC before calling a change warning-free.** The
maintainer builds on macOS with Clang, and GCC's -Wall does not enable
-Wunused-const-variable or catch the misleading indentation left when a
preprocessor block removes an `if ... else`. Two constants unused since the
spot finder moved into the tree, and a dangling `else` in the pinned-buffer
release, were invisible here and visible there. Clang is installable in this
container: `apt-get update && apt-get install -y --no-install-recommends clang`,
then `CC=clang CXX=clang++ cmake -S . -B /tmp/bclang`.

**Format before committing: `clang-format` for C++, `black` for Python.** Neither
was installed where most of this was written, and neither was run, so 77 of 105
C++ files and 16 of 46 Python files were out of format before anyone noticed.
The tree was formatted with clang-format 18.1.3 (`.clang-format`: LLVM) and
black 26.5.1; clang-format's output varies between versions, so a different
version reformats lines that were already right, and that churn is worth
avoiding. Both install here: `apt-get install -y --no-install-recommends
clang-format` and `pip install black --break-system-packages`. The CUDA and
Metal sources are not formatted, because nothing here can compile them to show
a format change left them working.

**A test run after a failed build tests the old binary.** `cmake --build` then
`mxi_tests` reports success from whatever was last linked if the build fails,
and the failure scrolls past above it. Check the build's exit status, or check
that the test count moved when a test was added.

**Prediction is the oracle.** Build order is prediction, then indexing, then
refinement, because prediction generates ground-truth reflection lists to test
the other two against.

## The convention trap, and what it caught

The closed-loop test **cannot** detect a wrong convention: predicting and
mapping back applies the error twice, once in each direction, and it cancels.
A green suite proves self-consistency only.

It hid two real errors until a real `.expt` arrived. `docs/conventions.md` has
the detail; the short version:

- **`s0 = -direction / wavelength`.** The sign was backwards. dxtbx points the
  beam direction back towards the source.
- **The parallax correction was missing entirely.** Worth 1.58 pixels on a
  0.45 mm sensor. Not a constant offset — a smooth function of scattering
  angle, which refinement partly absorbs and leaves as a radial residual.

Plus two structural surprises: the goniometer is multi-axis
(`axes`/`angles`/`scan_axis`, not `rotation_axis`/`fixed_rotation`) and
`scan.properties.oscillation` is a per-image array, not `[start, width]`.

`tests/test_real_geometry.cc` embeds forty real reflections and the geometry
that produced them, regenerable by `python/src/mxeq/fixtures/real_data.py`
(it was tests/make_real_data.py, until the Python moved into one package). Both parallax
directions reproduce DIALS to the last bit; `entering` agrees on 100% of 13072
reflections; the full chain lands on `A h` with median residual 1.9e-4.

The goniometer decomposition is now closed too, against the l-cysteine
four-sweep data in `tests/real_cysteine.hh`. Insulin could never have done it:
every setting angle there is zero, so both rotations are the identity and any
arrangement passes.

**Test the wrong versions, not just the right one.** `test_multi_axis.cc`
asserts that swapping fixed and setting is caught by exactly two sweeps, that
dropping the fixed rotation is caught by exactly two, and that folding the
scan-axis angle into the fixed rotation is caught by three -- including sweep
0, which is blind to the other two failures and is the only thing that pins
that rule. Without those, a green suite would not distinguish a correct
composition from a lucky one.

**Goniometer axes run from the sample outwards to the laboratory.** `axes[0]`
carries the sample; each subsequent axis carries the one before it; the last is
bolted to the floor. So an axis further out applies later and multiplies on the
**left**. The first implementation accumulated the other way, composing the
stack inside out -- and all 48 tests passed, because with one axis below the
scan axis the two orders are the same matrix. Pinned now against the physical
arrangement on a synthetic three-axis goniometer, including a check that the
two orders really do differ so the test is not vacuous.

**"Three circles" is not the condition for testing the composition order.** A
four-sweep phi/chi/omega dataset does put two axes below the scan axis, but phi
is zero in every sweep, so its rotation is the identity and both orders give
the same matrix to the last digit. What is needed is two axes below the scan
axis at *simultaneously non-zero* angles. Asserted, so that substituting better
data makes the test fail rather than silently start passing for a new reason.

**Still open, and a green suite says nothing about any of them:** the
composition order *against data*, `first_image != 1`, multi-panel detectors,
and the `hierarchy` block, which is ignored.

## Indexing: three things that produced plausible wrong answers

**Peak centroids are not a basis.** The FFT grid step is around 0.8 Angstrom,
so vectors straight off the transform are good to a few parts in a thousand,
which is a tenth of an index at the detector edge. The least squares fit of A
to the indexed reflections is what makes the cell numerically right rather than
merely recognisable. It tightens its tolerance over four rounds; fitting
throughout at the acceptance tolerance lets reflections a quarter of an index
out pull the basis as hard as ones that are exact.

**The fraction indexed is basis-dependent and cannot compare two bases of the
same lattice.** Acceptance asks whether every fractional index is within
tolerance, and a unimodular change of basis mixes the components: (0.2, 0.2,
0.2) in one basis is (0.4, 0, 0.2) in another describing the same lattice. A
guard that only accepted the reduction if it indexed at least as many
reflections duly rejected the correct cell. The invariant to test is the
volume.

**Coincident reciprocal lattice points are input, not noise.** Pooling several
sweeps puts the same reflection, measured at different goniometer settings, at
the same place in the crystal frame — that is the point of measuring more than
one sweep. Counting those pairs in the nearest-neighbour cell estimate made it
diverge to 9e15 Angstrom and find no candidates at all.

## Resolved: the 0.45 per cent cell edge

Refinement fixed it, and the explanation is the detector. Before refinement two
edges matched DIALS to 0.02 per cent and one was 0.45 per cent short; after
refining crystal and detector together all three agree with each other and the
asymmetry is gone. Refining with the detector *held fixed* leaves the asymmetry
in place, which is the evidence: the imported detector model was wrong, and the
cell was absorbing the error anisotropically because indexing had no other
parameter to put it in.

## Not the parallax correction, and how that was settled

A natural suspicion about the cell and distance discrepancy is that the
millimetre-to-pixel mapping is missing from the prediction path: omitting it
would make predictions fall short radially, and refinement would push the
detector further out to compensate, enlarging the cell. The magnitude fits --
the correction is 0.048 mm median, which at a typical 50 mm radius is about a
tenth of a per cent in distance, against the 0.16 observed.

It is not that. Predicting from DIALS' own refined model reproduces DIALS' own
`xyzcal.px` **exactly**, 0.0000 px median and maximum. Turning the correction
off displaces predictions by 0.68 px median, 1.0 px worst.

The test that should have settled this immediately had a tolerance of 0.5 px,
which is less than the effect it was guarding against -- so it would have
passed for more than half the reflections with the correction missing
altogether. It now asserts exactness, and a companion test measures what
switching the correction off does, so the bound is known to be sensitive to the
thing it guards. **A tolerance looser than the error it exists to catch is not
a weak test, it is not a test.**

## Residual structure: what is measured, and two claims retracted

There is real structure in the post-refinement residuals at the 0.2 px level.
Binned eight ways with 1600 reflections each, the medians are far above the
noise on the mean. `mxi_residuals` prints it.

**Retracted: "x and y disagree about the distance by 0.6 mm."** That came from
scanning the detector distance while holding the *whole panel rigid*, so the
crystal had to absorb every lateral and angular error too. Measured properly --
fitting the radial residual against radius separately along each axis, with the
panel refined -- the implied distance differs between fast and slow by 0.04 mm,
not 0.6. The 0.6 was an artefact of the scan's own constraint.

**Partly retracted: "a clean radial gradient" on this refinement's own
residuals.** There the tangential component mirrors the radial almost exactly,
which is a directional pattern pushed through a radial decomposition. On DIALS'
own output, with the affine part removed, the radial gradient is clean and the
claim stands -- but it stands on DIALS' residuals, not on this code's, and the
two should not have been conflated.

**Refuted: the absorption depth.** Eqn (6) is the unconditional first moment
and counts photons that pass through the sensor as contributing depth zero. The
conditional mean, dividing by 1 - exp(-mu t), is 24 per cent larger at normal
incidence and varies with angle, so it cannot be absorbed into the distance.
Implemented behind `parallax_conditional` and tested: rmsd goes from
0.3514/0.3388 to 0.3516/0.3409 and the cell moves by 0.015 per cent. It is not
the cause.

## The radial residual, measured properly

Measured on DIALS' own refined insulin model, 12907 indexed reflections, so it
is DIALS' model inadequacy with none of this code's mixed in. `mxi_residuals`
prints all of it.

**It is real and it is radial.** Median radial residual runs monotonically from
-0.155 px at the beam centre to +0.067 px at the edge, with a fitted slope of
+3.0e-4 px per px, about 0.48 px across the radius range.

**No change to the detector can remove it.** The affine part of the residual
field -- two scales, a rotation about the beam, a shear, a translation, which
is everything a flat panel can express -- comes to almost nothing: scales of
0.01 per cent, 0.07 mrad of rotation, 0.035 mrad of shear, worth 0.054 px at
the panel corner. Removing it leaves the radial gradient essentially untouched,
-0.155 to +0.067. Refinement was free to move the distance and did; the
gradient is not a distance error.

**It is geometric, not a centroid artefact.** Fitting the slope within each
intensity quartile separately gives 3.4, 2.9, 3.3 and 4.7e-4 across a factor of
fifty in intensity. A centroid-estimation bias -- thresholding, background,
shoebox truncation -- would scale strongly with signal to noise. It does not.
The mild rise in the strongest quartile tracks bounding-box depth (3.2e-4 at
one image, 5.0e-4 at seven) and is plausibly second order.

**Not the absorption depth.** See below; tested and refuted.

**Not fluorescence escape either.** Si K-alpha is 1.74 keV. A 13 keV photon
that loses it still deposits 11.26 keV, far above a threshold set at half the
photon energy, so it is counted at its own pixel -- and the escaped 1.74 keV
photon falls below threshold and is never counted at all. A photon-counting
detector with a half-energy threshold is immune to this by construction.

**Partly the module tiling.** This is an Eiger2 16M: 4 x 1028 + 3 x 12 = 4148
fast, 8 x 512 + 7 x 38 = 4362 slow, confirmed against the data by zero spots in
every predicted gap out of 13766. DIALS imports it as ONE flat panel, so
nothing about the individual modules' positions or tilts can be represented,
and whatever is wrong with the tiling has nowhere to go but the residuals.

`mxi_residuals --modules 1028,12,512,38` shows it as a step at each boundary:

    fast at 1028   step -0.127 px
    fast at 2068   step -0.203
    fast at 3108   step -0.121

Three boundaries out of three, same sign, 0.12 to 0.20 px, each many times the
standard error on medians of hundreds. The slow boundaries are four of five
negative and noisier.

How much it accounts for, with controls, on median |residual| of 0.2795 px:

    per-module constant            0.2412  (13.7 per cent)
    per-module affine              0.1856  (33.6)
    control, diagonal bands        0.2508  (10.3)
    control, random groups         0.2800  (-0.2)

The control matters. A per-module constant barely beats carving the detector
into the same number of arbitrary diagonal bands, so rigid module offsets are
not the right description -- it is mostly absorbing smooth spatial structure.
The per-module *affine* is much better, which points at each module having its
own tilt rather than its own position.

Removing per-module constants also drops the radial slope from 3.0e-4 to
1.9e-4, so tiling is about a third of the radial gradient and something else is
the rest.

## Where the cell and distance difference from DIALS actually came from

Not the refinement weights, which was the obvious suspect. On this detector the
centroid variances are nearly constant -- median 0.086, 0.087, 0.085, with a
floor at 1/12, the variance of a uniform distribution over one pixel -- so
statistical weighting is very nearly unit weighting. Measured: unit weights
give a distance of 170.080 against 170.088 for statistical. Eight microns.

**Given DIALS' reflection set, this refinement reproduces DIALS.**

    refinement input          distance     cell
    DIALS' indexed.refl       170.0714     67.424 67.468 67.475
    DIALS' own answer         170.0730     67.425 67.474 67.468
    our indexed.refl          170.3452     67.526 67.560 67.605

Two microns in distance and six thousandths of an Angstrom in cell. The entire
discrepancy is in the indexing, not the refinement.

**What the indexer admits that DIALS does not:** 380 reflections, all of them
weak. Median intensity 24 against 241 for the rest, median n_signal 4 pixels
against 27. The 13064 both index agree on every single Miller index.

**Outlier rejection cannot remove them.** At 2, 2.5, 3 and 4 sigma the refined
distance is 170.348, 170.347, 170.345, 170.344 -- it does not move. That is the
part worth remembering: *a misindexed population that is internally consistent
does not produce outliers, it moves the model.* Refinement adjusts until those
reflections fit, and then nothing about them looks anomalous. Rejection only
catches reflections that disagree with their neighbours, and these agree with
each other.

Tightening the assignment tolerance does work -- 0.30 and 0.20 both give
170.35, 0.15 gives 170.076, 0.10 gives 170.038 -- but lowering the default is
the wrong conclusion, because DIALS' own `hkl_tolerance` is also 0.3. The
difference is that DIALS *iterates*: it assigns, refines, re-assigns against
the improved model, and discards what stops fitting. This indexer assigns once.

**The fix is a macrocycle**, not a tighter tolerance, and it is now in:
assign, refine on the stronger half, re-assign under the improved model,
repeat. Strength is measured against the dataset's own median so the criterion
travels between detectors and spot finders.

    stage                                distance    volume vs DIALS
    assign once, refine on everything    170.345        +0.55 %
    indexing macrocycles                 170.101        +0.09 %
    ... and mxi_refine --strong-only     170.092        +0.08 %
    DIALS                                170.073          --

It converges in one cycle on this data and `rmsd_index` falls from 0.0757 to
0.0318. The same option exists on `mxi_refine`, because refining on everything
afterwards puts the bias straight back -- 170.101 becomes 170.243. It is off by
default there: throwing away half the data should have to be asked for.

## Is the radial residual the same in DIALS and here? Partly

Compared directly, on the same 12896 reflections, both models' residuals binned
by radius:

     radius      DIALS     here        radius      DIALS     here
     50- 150    -0.270   -0.250       800-1000    -0.007   +0.011
    150- 250    -0.227   -0.214      1000-1200    +0.051   -0.003
    250- 350    -0.184   -0.156      1200-1500    +0.098   -0.082
    350- 450    -0.140   -0.096
    450- 600    -0.070   -0.031
    600- 800    -0.045   +0.006

**Shared and robust:** an inward radial displacement at low angle, reaching
-0.27 px nearest the beam and decaying monotonically to zero by about 800 px.
The two models agree to 0.02 px in the innermost bins. It is independent of
signal to noise -- -0.125, -0.108, -0.109, -0.117 px across a factor of 4.7 in
I/sigma -- so it is geometric, not a centroiding or background bias.

**Not shared:** everything beyond about 800 px. Excluding the innermost bin the
two profiles correlate at 0.12, and at 0.03 after the affine part is removed;
over the outer eight bins alone the correlation is -0.16, and at 1200-1500 px
the two models disagree in sign. Peak to peak, DIALS 0.235 px against 0.156
here.

**So the single "radial slope" figure describes neither model well.** It mixes
one real shared feature at low angle with refinement-specific structure at high
angle, and its value -- +3.0e-4 for DIALS, +1.09e-4 here -- depends mostly on
where each refinement put the detector. Earlier entries in this file quote that
slope as though it were one phenomenon. It is not.

What is genuinely unexplained is therefore narrower than previously stated: an
intensity-independent inward radial displacement of a quarter of a pixel within
about twenty degrees of the beam, decaying to nothing beyond it. Module tiling
is a separate matter and lives at all radii.


### It is a fixed hardware fingerprint

The steps repeat. Across four independently refined sweeps of a separate
four-sweep insulin experiment on the same detector, the slow-axis pattern comes
out the same shape every time:

    boundary   sweep 0   sweep 1   sweep 2   sweep 3
      1062      +0.045    +0.070    +0.022    +0.090
      1612      -0.187    -0.097    -0.128    -0.142
      2162      -0.321    -0.362    -0.366    -0.326
      2712      +0.354    +0.239    +0.330    +0.361
      3262      -0.319    -0.132    -0.317    -0.348

and comparing that experiment against the single-sweep one -- different data,
different crystal orientation, separately refined detector models -- gives a
**correlation of 0.950 over seven boundaries, with an rms difference of 0.063
px**. The steps themselves are 0.229 px rms, which at 75 micron pixels is
**17 microns**.

Seventeen microns is a module placement tolerance, not a modelling error. This
is a property of the hardware, measurable from diffraction data, stable between
experiments, and currently absorbed into the residuals because the detector is
described as one flat panel.

**The actionable part:** a multi-panel model, one panel per module, would give
refinement somewhere to put this. That is a change to how the NXmx file is
imported, not to any algorithm. The per-module affine result suggests the
modules want tilts and not just positions.

### The falsification test, passed

l-cysteine on a PILATUS 2M at I19: 1475 x 1679, 172 micron pixels,
3 x 487 + 2 x 7 fast and 8 x 195 + 7 x 17 slow, confirmed from the data by zero
spots in every predicted gap out of 15477 where arbitrary bands of the same
width hold 46 and 141.

**It shows no module steps.** Median absolute step at its nine real boundaries
is 0.030 px, against 0.038 at sixty arbitrary slow rows and 0.024 at forty
arbitrary fast columns. The per-sweep values flip sign between sweeps. Nothing
at the Eiger2 spacings either, as there should not be.

The test has the power to have found it. With around four hundred reflections
each side and a 0.56 px spread the standard error on a median step is about
0.035 px, so an Eiger2-sized step of 0.12 to 0.35 px would have shown at four
to ten sigma.

One thing had to be fixed first. As supplied, this dataset is indexed with one
UB across all four sweeps and perfect goniometry, which is false, and the
residuals are dominated by it: median 1.05 px, with only 5218 of 11686
reflections inside a 1.5 px clip. Refining with a crystal per sweep --
`mxi_refine --separate` -- takes it to 0.56 px with 10791 reflections, and only
then is the test sensitive enough to mean anything. **A null result from an
underpowered test is not evidence of absence**, and this one would have been
underpowered by a factor of two.

So the fingerprint is detector-specific: present and reproducible on the Eiger2
at 0.229 px rms (17 microns), absent on the PILATUS 2M below about 0.06 px
(10 microns). An artefact of the analysis would appear on both.

So: a radial, geometric, intensity-independent displacement of about a fifth of
a pixel peak to peak, in a direction the detector model cannot represent.
Candidates not yet tested: the parallax offset is evaluated at the front-face
intersection rather than solved self-consistently (matching DIALS exactly is no
defence -- both could be wrong the same way, though the one-shot to converged
difference measured 0.0005 px, which is too small); a non-planar sensor; and
the interaction depth distribution being something other than a single
exponential, which for silicon near an absorption edge it is not.

## The degeneracy that replaces it

Detector distance and cell scale are very nearly degenerate. A detector further
away with a proportionally larger cell predicts the same spots in the same
places, so the two refinements settle at different points on a flat valley:

    distance  170.34 mm (here)  vs  170.07 (DIALS)   0.16 per cent
    cell      +0.18 per cent per edge, +0.55 in volume

which is arithmetically the same statement twice. Both fit the observations
about equally well, and through the reindexing operator the two models predict
the same positions to 0.13 px median. So this is not a bug to find; it is a
direction the data barely constrains. Anything that claims to have fixed it
should be checked against whether the residual actually improved.

Refining the beam as well moves it by almost nothing (170.340 against 170.344),
so the beam is not what breaks the degeneracy.

## Scan-varying refinement

Control points in A, evenly spaced over the scan, interpolated by a uniform
cubic B-spline clamped at both ends. Not DIALS' Gaussian smoother, which
weights three sample points per observation; both are smooth and local, and
they are not the same model.

It was linear at first. Cubic because a crystal does not change direction
abruptly at arbitrary points in a scan, and linear interpolation says it does --
continuous, but with a derivative that jumps at every control point. A B-spline
is C2.

B-spline rather than an interpolating spline because the support is local: four
control points per evaluation, so the Jacobian is banded rather than dense.
That is worth more than it sounds -- see `docs/gpu.md`.

The cost is that in the interior a control point is a coefficient, not the
value of the model at its own position. Only the two ends interpolate, by the
clamping. Anything reading an interior control point as "the setting matrix
there" is wrong.

What changed and what did not: on real l-cysteine the fit is identical, 0.2101
px at nine control points against 0.2100 for linear and DIALS' 0.2096. On
synthetic data with a planted drift the spline recovers it exactly at five
control points where linear needed three and could only approximate -- because
the clamped ends are pinned by the data near them, where linear's outermost
points wandered.

It matches DIALS. On l-cysteine, median |d| against number of control points:

    1 (static)  0.590      5  0.211
    3           0.252      9  0.210
                          15  0.217     DIALS (Gaussian smoother)  0.2096

Five to nine points reproduces DIALS to the fourth decimal, and fifteen is
worse -- the drift is captured and what is left is noise being fitted.

**Static first, always.** Control points started from an unrefined model
absorb errors that belong to the detector, fit well, and mean nothing.

**Outlier rejection is required, not optional.** About four per cent of
reflections have forward and reverse maps that disagree under a scan-varying
model. This was attributed to near-tangential geometry; **that attribution is
withdrawn**. The paper's own criterion for near-tangential -- the volume
(e x r_phi) . s0 below 0.05 -- removes 5.4 per cent of reflections and only 12
per cent of the disagreements. Nor is it convergence of the forward iteration,
nor reflections whose two Ewald roots are close. The population is
unexplained; rejection removes it and refinement then works, which is a
workaround. Once the model varies over the scan, the forward
and reverse maps select different roots for those, and their residuals run to
tens of images. With them in, a planted 0.5 degree drift is not recovered at
all; with them rejected it comes back as 0.5000. The fraction is 3.5 per cent
at a drift of 0.02 degrees and 4.3 at 0.5, so it is a property of those
reflections and not of how much the crystal moved.

**The outermost control points were weakly constrained under linear
interpolation**, by whatever lay at the very ends of the scan. Clamping the
B-spline fixed that: the end control points are on the curve, so the data near
the ends pins them directly.

## A bug the round-trip test could not see

The forward map (`predict`) and the reverse map (`centroid_residual`) must be
one model seen from two directions. For a scan-varying crystal they were not:
the iteration that finds the setting matrix was seeded from the start of the
scan for *both* Ewald roots, so a reflection late in the sweep converged on the
wrong one. Forward and reverse disagreed by 0.62 px rms, and since refinement
minimises the reverse map against data made by the forward one, no amount of
refining could reach the truth.

`prediction_round_trips_through_reciprocal_space` could not catch it: it goes
through `reciprocal_lattice_point`, not through the target function refinement
actually minimises. **Test the function being minimised, not a cousin of it.**
`prediction_and_the_refinement_target_agree_exactly` now does.

It was found by asking the dullest possible question -- does the truth model
reproduce its own predictions -- which should have been the first test written
and was not.

## Refinement notes

**Derivatives are numerical, deliberately.** Fifteen parameters over thirteen
thousand reflections is sixteen predictions per iteration and takes a second. A
wrong analytical derivative does not crash — it converges smoothly to the wrong
answer and reports a small residual doing it. If profiling ever demands the
analytical version, this one stays as the thing it is checked against.

**Choose the Ewald root by proximity to the observation, not by the `entering`
flag.** The flag assumes a model already close enough to trust, which at the
start of refinement it is not. Choosing wrongly puts a reflection tens of
images away, so there is a test for it.

**Rotate the detector about the panel centre, not the laboratory origin.** A
rotation about a point 200 mm away is mostly a translation, and the two
parameter groups would then be so correlated that the normal matrix is nearly
singular.

**Reject outliers between macrocycles, never within one.** Rejecting while the
model is still moving throws away reflections for being far from a prediction
that was wrong.

## Planned: how indexing should handle multiple sweeps

Index across **all** sweeps at once, then immediately split and refine against
each sweep individually, keeping the bulk matrix common.

The reason is visible in the l-cysteine residuals. A joint index must assume
one UB and perfect goniometry, and the goniometer does not return to precisely
the same place between sweeps, so no single matrix fits all four -- residuals
come out a hundredfold worse than insulin. That is not a failure of indexing;
it is the constraint doing what it must. Joint indexing is still the right
first step, because it is what guarantees a consistent basis across sweeps and
avoids four independently-chosen and mutually reindexed lattices. The
constraint is then broken at the first opportunity rather than carried forward.

Consequence for the code: the indexer works on the pooled reciprocal lattice
points from every sweep, but each sweep keeps its own goniometer, scan and
detector throughout, and the output is a list of experiments sharing a crystal
rather than one experiment.

## Things got right for a reason

**The setting matrix is the plain inverse of the real-space matrix.** Real
space vectors are the rows, so their product with A is the identity. The
inverse transpose passes on any cubic cell because both are diagonal; the test
uses a triclinic cell for exactly this reason. Same bug, same fix, as in
`mxeq`.

**Panel intersection has a nanopixel edge tolerance.** A ray aimed at the exact
corner of a panel returns -4e-13 pixels, because the intersection is solved in
millimetres and divided back. A bare `< 0` test rejects it. That loses
reflections predicted on an edge, and on a tiled detector loses rays striking
the seam between two panels, which then belong to neither. A nanopixel is nine
orders below anything physical and three above the round-off.

**The goniometer inverse negates the angle rather than inverting the matrix.**
Exact for the rotation part, and avoids a determinant on a matrix that is
orthogonal by construction. It does assume `fixed` and `setting` are
orthogonal, which is true for goniometers and is not checked.

**`Mat3::inverse` sets a flag instead of returning infinities.** A caller that
ignores it gets an answer that looks wrong, rather than NaNs that propagate
silently through a refinement and come out as a plausible-looking cell.

## Command lines take options in any position

All three programs originally read their file names from `argv[1]` and
`argv[2]` and scanned for options from `argv[3]` onwards, so

    mxi_refine --beam indexed.expt indexed.refl

took `--beam` as the experiment list, `indexed.expt` as the reflections, and
reported "unknown option 'indexed.refl'". The error named an argument that was
not the problem, which is worse than no message.

`src/args.hh` splits a command line into positionals and options wherever they
appear, and is in `src/` rather than `apps/` so that it can be tested. Argument
parsing is exactly the kind of code that never gets tested because it looks too
simple to get wrong.

An unknown option is refused rather than ignored. Silently accepting a
misspelling is how a run comes to use settings nobody chose -- `--strong-ponly`
would otherwise have been dropped and the refinement would have used every
reflection while the operator believed otherwise.

## Scan points are boundaries, so there are N + 1 of them

DIALS samples A at the start of the scan and at the end of every image: 1801
samples for 1800 images. This wrote 1800, which is what a per-image reading of
the name suggests, and dials.export asked for the point after the last image:

    DXTBX_ASSERT(index < A_at_scan_points_.size()) failure

**And `A_points` was doing double duty.** A refinement produces a handful of
spline control points; a file carries the evaluated samples. Re-evaluating the
spline over samples smooths them again -- 4e-5 relative on a real 1801-point
model, 0.003 Angstrom on a 67 Angstrom cell, compounding on every read and
write. Samples now carry a flag and go back out untouched. Reading and writing
somebody else's model must not change it.

The first version of that branch used an early `continue`, which skipped the
code that interns the crystal and appends the experiment: the file came out
with no experiments at all. **Every test passed**, because not one of them wrote
a scan-varying model. There are now three that do.

**The flag was only half the rule.** Later found: `A_at` re-smoothed samples as
control points on every READ, whatever the flag said (predictions from DIALS'
model 0.0003 images from DIALS' own); refinement replaced samples with control
points without clearing the flag, so a model refined twice was written with
five scan points for three hundred images; and a static refinement of a
scan-varying crystal refined an A that predicted nothing. The rules now: samples
are interpolated linearly between image boundaries, a scan-varying refinement
always converts to control points and clears the flag, a static refinement of
the crystal makes it static, and the reader refuses any count but N + 1.

## What we cannot decode is not ours to drop

Shoeboxes were read as opaque bytes and thrown away on write. The argument was
that nothing here can subset a shoebox and one that silently stopped matching
its table would be worse than its absence. The argument is sound; the
conclusion was not.

Indexing and refinement do not remove rows. They add columns and set flags, so
the bytes stay correct, and dropping them turned an 84 MB table into a 20 MB
one that `dials.integrate` refuses outright:

    Error: shoebox data missing from reflection table

They are now kept and written back verbatim -- 13604086 bytes in, the same
sha256 out, through `mxi_index` and `mxi_refine` both.

The condition the original argument was right about is checked rather than
assumed: an opaque column whose row count no longer matches the table cannot be
written, and that **throws** rather than dropping it. Doing it silently a second
time, for a better reason, would be no better than the first time.

Third instance of the same pattern in this file, after the `.expt` blocks and
the flags: **a format is not only the parts of it you understand.**

## dials.* asks the flags, not the Miller indices

The `flags` column is a bitmask and every DIALS tool filters on it. Indexing
here set correct Miller indices and left the flags alone, so the output
processed perfectly and was then invisible to every selection downstream --
which is a failure mode worth naming, because nothing errored.

`mxi_index` now sets `indexed` (1 << 2) and `mxi_refine` sets
`used_in_refinement` (1 << 3), both cleared as well as set so that the flag and
the Miller index cannot disagree. An indexed insulin row now carries 36, strong
| indexed, which is what a real DIALS indexed.refl carries; after refinement
the fitted ones carry 44.

**`Table::int_column` REPLACES a column with a zeroed one.** That is right for
a derived column recomputed in full -- leaving stale predictions in the rows a
pass skips would be worse than clearing them -- and wrong for any column that
must be read before it is written. Setting the indexed bit that way threw away
the strong bit dials.find_spots had set, and nothing failed until something
filtered on it. `modify_int_column` returns the existing column;
`flags` is the one place that needs it.

`RefineResult::rows_used` reports which rows the fit actually used, after
outlier rejection and after the ill-conditioned ones were dropped, because
which reflections a residual was averaged over is not a detail and this is
where DIALS records it.

`centroid_outlier` is bit 17, settled from a real DIALS indexed.refl rather
than guessed: 8230 rows carry it, every one indexed, none also marked
used_in_refinement, and their median |xyzcal - xyzobs| is 0.843 px against
0.310 for the rest. The four values DIALS writes are 32, 36, 44 and 131108, and
this now writes the same four.

DIALS sets `used_in_refinement` on exactly 18000 of 75406 indexed reflections
on a 180 degree scan -- a hundred per degree, its sampling default. This sets
it on everything it fitted, which is a real difference in what the flag means
between the two and is not a defect in either.

## Hold the detector during the scan-varying pass

A scan-varying crystal and a refinable detector distance are degenerate in
position. Measured on a truth that fits exactly, both scaled by two tenths of a
per cent:

    distance and cell together     0.062  0.072  0.318
    distance alone                 0.778  1.100  0.000
    cell alone                     0.836  1.152  0.318

Twelve times smaller in position when they move together. The rotation angle is
*not* degenerate -- it sees the cell and not the distance -- which is why a
static refinement pins the pair, and why a scan-varying one does not: a hundred
and sixty crystal parameters can absorb the angular residual by drifting the
orientation, leaving the cell and the distance free to slide together.

On 1800 images of insulin, eighteen control points:

    --beam                       distance 169.977  V 236088   rmsd 0.309 0.295 0.320
    --beam --scan-varying 18     distance 169.705  V 234704   rmsd 0.233 0.218 0.193
    ... detector held            distance 170.010  V 236014   rmsd 0.231 0.234 0.192
    dials.refine                                   V 236092   rmsd 0.206 0.208 0.203

The distance drifted 0.27 mm and the volume fell 0.59 per cent, for five
thousandths of a pixel. Holding the detector brings the cell to within 0.013 per
cent of DIALS' edges and 0.03 per cent of its volume.

**The static pass has already placed the detector**, which is what makes
holding it safe rather than a constraint on a quantity nobody has determined. A
detector does not move during a sweep, so there is nothing scan-varying about
it; letting it move while the crystal is free only gives crystal drift
somewhere else to go. `--detector-in-scan-varying` restores the old behaviour.

A first attempt to test this refined a static truth with the detector free and
with it held, and asserted that holding it kept the cell closer. It does not,
and should not: with exact data and a detector starting in the wrong place,
refining it recovers the truth exactly and holding it forces the error into the
cell. The degeneracy is invisible there because there is a unique exact answer.
The test now measures the degenerate direction itself.

## The rotation angle has to be determined before it is worth fitting

Waterman eqn (40) divides by the volume of the parallelepiped formed by the
rotation axis, the reciprocal lattice vector and the beam. Near the rotation
axis it goes to zero: the angle at which such a reflection diffracts is
arbitrarily sensitive to the model, and the Lorentz factor has the same
asymptote, so its observed angular centroid is poorly determined too. Fitting
them puts noise into the rotation-angle residual that no model can remove.

DIALS discards below 0.05 by default. This did not, and the cost shows in the
rotation-angle residual first:

    insulin, 30 degrees, static     0.310 0.290 0.245  ->  0.274 0.254 0.220
    l-cysteine, 170 deg, sv 9       0.168 0.212 0.133  ->  0.161 0.194 0.124

The fraction removed depends on the geometry, not only on the cutoff: 2.4 per
cent of indexed insulin over thirty degrees, 0.1 per cent of l-cysteine, whose
small cell puts every reflection far from the axis.

**A residual averaged over a different set of reflections compares with
nothing**, so `mxi_refine` now says how many it dropped and how many it
averaged over. dials.refine applies the same cutoff and reports over what is
left, which is part of why its numbers looked better than they were.

## Which model each observed column was computed through

`dials.refine` stops without `xyzobs.mm.value`, and reasonably: DIALS measures
its residual in millimetres and radians -- Waterman eqn (25) is in X, Y and phi
-- so that column is the observation it minimises against. `dials.index`
produces it and this did not.

The columns split by **which model they were computed through**, and the split
was measured against DIALS' own output rather than assumed:

    xyzobs.mm.value       the model as IMPORTED    median difference 3e-17
    xyzobs.mm.variance    the model as imported    essentially exact
    s1, rlp, entering     the CURRENT model        3e-3 from imported,
                                                   9e-5 from DIALS' refined

`dials.find_spots` writes the millimetre centroids through the detector as
imported and nothing recomputes them; `dials.index` calls
map_centroids_to_reciprocal_space after refining, so s1 and rlp follow the
refined model. **That is why xyzobs.mm.value is stale after refinement and must
never be a join key** -- a fact recorded here long before the reason for it was
understood.

Recomputing the millimetre centroids from the refined model would look like a
correction and would be a silent change to the observations refinement had just
been fitted to.

`entering` is `s1 . (m2 x s0) > 0`, the opposite sign to the test as usually
written, because s0 here points source to sample and dxtbx's beam direction
points the other way. Determined against DIALS' flags: 13760 of 13766 agreed,
the rest having a triple product within rounding of zero. The wrong sign gave
six.

## Write back everything that was there, not only what we model

`dials.refine` then failed with

    IndexError: list index out of range        experiment_list.py:603

because the experiment said `"imageset": 0` and the imageset list had been
written out empty. The imageset block is the only link from an experiment list
to the images; `profile`, `scaling_model` and `history` are equally not this
package's to discard.

`ExperimentList` now keeps the document it was read from, and writing replaces
the models it understands and leaves everything else alone. A list built in
memory has no source, and then the absent models are **null** -- never `-1`,
which is a valid index into a Python list and reaches for the last element of
one that may be empty.

Carrying a block through then exposed the other half of the number problem
below: the reader did not record whether a value had been written with a
decimal point, so an `"imageset": 0` came back out as `"imageset": 0.0` and
dxtbx indexed a list with it. Both halves are needed, and each is useless
alone.

## `0` is not `0.0`, and dxtbx knows the difference

`dials.refine` refused an `indexed.expt` written here:

    DXTBX_ASSERT(obj_type == "float") failure      dxtbx scan.cc:80

A scan starting at zero degrees has an oscillation array beginning 0.0, and the
JSON writer collapsed any whole-valued double to an integer, so the file said
`[0, 0.1, 0.2, ...]`. dxtbx reads the type of an array from its first element
and requires float.

The fix is not to write every number with a decimal point -- `image_range` of
`[1.0, 300.0]` would be just as wrong. Whether a number is a count or a
measurement is known at the point it is constructed and nowhere else, so
`json::Value` now carries it: built from `int`, `long`, `long long` or
`std::size_t` it writes without a point, built from `double` it writes with one.

**Nothing here could have caught this, and that is the recurring lesson.** The
document round-tripped through this reader perfectly, because this reader parses
`0` and `0.0` into the same double. A round-trip test cannot detect a wrong
convention; only a file written or read by something else can. That is now five
times: the `.refl` container shape, the geometry conventions, the scan
oscillation form, the number types, and the blocks dropped on write.

`tests/test_expt_format.cc` therefore asserts on characters rather than values,
and was checked by reverting the fix -- three of its tests fail without it.

## An instrument that changes what it measures

The phase timers in indexing were added with a regular expression that inserted
the closing assignment before every `return` in the function. One of those
returns was inside a lambda -- the accessor that reads a grid point -- which
the peak search calls about four hundred and fifty million times. Every call
read the clock.

So the peak search was reported at three seconds when it takes seven hundred
milliseconds, and the extra was the instrument. The number was then used to
argue about what to optimise, and the answer would have been to parallelise the
clock.

This is the GPU benchmark bug again, which passed ordinary `std::vector` memory
where the device wanted `gpu::host_alloc` and inflated every timing equally.
Both times the measuring apparatus was wrong and the conclusion drawn from it
was confident. **Check what an instrument costs before believing what it
reports**, and be suspicious of any edit that places code by pattern rather
than by position: a `return` inside a lambda is not the function's return.

## Render the model into a measurement; do not compare a measurement to a model

Spots on this detector are two pixels across. Every comparison of the profile
model with the data took a number from the data -- through a pixel grid that
truncates and quantises it, inside a mask that cuts it off -- and a number from
the model in closed form. At that scale the difference between the two ways of
computing is most of what is being measured.

It produced three wrong conclusions in a row, each argued carefully from what
had been measured: that one sigma held 0.74 of the density where a Gaussian
holds 0.47, so the model was far too wide; that the spots were anisotropic by a
factor of 1.24; and that the anisotropy was not the sensor, because the excess
was flat in obliquity where absorption predicts growth.

Integrating the same model over the same pixels, over the same images, inside
the same mask, and reducing both the same way:

    positions      agree to a hundredth of a pixel
    widths         agree to within four per cent
    anisotropy     observed 1.074, model 1.092 -- predicted, not missing

Whatever the grid does to the data it does to the model, and what is left is
the model being wrong. It was barely wrong. **When the thing being measured is
at the sampling scale, the comparison has to be made in the data's own terms.**

## A silent replace leaves two documents that disagree

`README.md` sat for a dozen commits as a seventy per cent duplicate of this
file, still titled `mxi-index`, describing none of the six tools
written after it. Edits meant for it had been written against text that was not
in it, matched nothing, and changed nothing -- the same failure that let a
reformatted line silently escape a fix in `mxeq` twice.

An edit that matches nothing is not a no-op, it is a change that did not
happen while the commit says it did. **Anything that rewrites a file must
assert that its anchor was found**, and fail loudly when it was not.

Written down, and then done again three commits later: a regex meant to add a
test target to spotfinder/CMakeLists.txt matched nothing, the commit said the
test was wired in, and it was not. Knowing the failure mode is not the same as
guarding against it. The guard is an assertion in the edit, every time.

## Relative precision belongs to the arithmetic, not to the answer

This file once claimed that float32 would leave about four digits in a
numerical derivative of the target. It leaves none. The reasoning was that a
1e-6 relative parameter step changes the residual by 1.5e-3 of its own size and
float epsilon is 1.2e-7, so four digits survive.

The residual is a difference of detector positions of order two thousand
pixels. Its absolute error in float32 is epsilon times the *position*, about
2.4e-4 px, not epsilon times the residual. The change being measured is about
5e-4 px. Signal and noise are the same size.

Measured rather than estimated, in `tests/test_precision.cc`: median relative
error 1.00 at the step the refinement uses, and 1.4e-2 at the best step float
can manage. The analytical derivative in the same precision gives 1.3e-7.
Analytical derivatives are therefore a precondition for a device port and not
an optimisation.

Both `target.h` and `derivatives_t.h` are templated on the scalar type for this
reason, and both are checked against the double-precision originals first --
the derivatives bit for bit. Writing a second implementation to measure the
first is worse than useless: a disagreement could be either thing.

## Integration: what has to stay true

The full account is `docs/integration_history.md`; `docs/integration.md` is the
reference. These are the rules that came out of it, each one broken once.

**Validity and region are separate mask bits.** A pixel was measured if and
only if `mask & kValid`. Never test `mask == 0` for "not measured": a bad pixel
keeps its region bit so the fit knows part of a reflection is missing, and
when that changed five places broke together, including the transform that
learns profiles, which began putting bad pixels on the grid as real zeroes.

**Flags are DIALS', read from `dials/array_family/reflection_table.h`.** Not
from memory, and not from a table here: bit 19 was once taken for an exclusion
flag, and `mxeq`'s own table had three wrong entries for as long as nothing
read it. A foreground reaching a masked pixel gets
`foreground_includes_bad_pixels` and `failed_during_summation` and NOT
`integrated_sum`, because the sum is missing counts; setting it anyway sent
truncated sums into `dials.scale`, which rejected them along every module edge.

**Files DIALS reads follow DIALS' conventions,** even where the internal ones
differ. Saved shoeboxes carry no region bit on an invalid voxel; DIALS rejects
a table that does.

**An angle goes into the scan before the crystal is looked up at it.** The
Ewald solver answers in any 2 pi interval. At -260 degrees on a 0 to 180 scan
the crystal came from image zero, and leaving reflections were up to a frame
late from the first scan-varying version onward. Wrap into a window centred on the scan: one that
starts at the scan start sends a root just before it to the far end.

**A scan-varying test needs a crystal that moves.** The same A at every scan
point is the static geometry with a scan-varying label; two separate bugs
passed tests written that way. Test physical invariants where they exist --
every prediction on the Ewald sphere under its own crystal; a pull of one on
spots planted at known positions -- since those cannot be satisfied by an
answer that is merely consistent.

**A variance is judged by the deviations it predicts, not by a floor.** It was
written here that a fitted variance below about half the summed one is a
missing term. A term was missing -- equation 34's background -- but put in as
a sum's it was several times too large, and the weak reflections' fitted
deviations between equivalents then had variance 0.58 where the summed ones
had 1.04. Carried through the fit's own weights, the weakest fitted variance is
a quarter of the summed, and its deviations 1.04: what the data say the fit
achieves. Check a variance against equivalents' exact deviations, or bare
background (a crystal turned off its spots), never against a remembered bound.

**A weight must not move with the observation's own noise.** A variance from
an observation's own counts -- the counted foreground, the fit's own intensity
-- gives a low fluctuation a small variance and a large weight, and an
inverse-variance weighted mean leans low: at a multiplicity of 210, a merged
I/sigma of -5.7 where CC half said there was signal. Variances from expected
counts in integration, from the model's intensity in scaling.

**DIALS is a reference, not the truth.** Where they disagree, find something
that can judge both: a reflection's symmetry equivalents judged DIALS' gap
recoveries worse than ours. The comparison says two programs differ, not which
is right.

**Before comparing, check the model and the images are one dataset,** and read
the matched count before anything else. One sweep's images integrated with
another's model gives reasonable-looking numbers, and a whole stretch of the
real-data measurements made here was of that kind. A comparison whose
matched count is far below the smaller table's rows is measuring the matcher.

**A median over a mixture hides a split.** Entering and leaving reflections,
pooled, gave a prediction difference of 0.0007 images; split, it was 0.0003
and 0.227. Split by the obvious categories before concluding two things agree.

**No per-thread state through `thread_local` in a parallel call.** The calling
thread takes part in every call and outlives them all, so a lane it chose
once was shared with a later call's new threads -- a data race in profile
learning that ThreadSanitizer, watching one run, did not see. Index per-thread
state by the worker number `in_parallel_by_worker` passes; and make a reduction
add in an order fixed by the data, so that the thread count changes no byte of
the output.

**When a fast version replaces a slow one, keep the slow one as the
specification** and test the fast one against it. `profile_on_pixels_direct`
is that for `profile_on_pixels`; the fast one's first version differed in 147
of 3000 boxes, which a test with a flat profile could not have seen.

**A number measured is a number attributed.** Performance and agreement figures
in the documents name their dataset. The documents that said "integration not
started" were accurate when written; nothing marked them as dated.

## One pass over the images: what keeps it the two passes' answer

`mxi_integrate` fits in the same pass it learns in, and the table is the two
passes' byte for byte (`python/tests/test_single_pass.py`). What that rests on:

* **A block is final only when every box learning into it has closed** -- the
  latest closing frame of the boxes whose learning region is in that block,
  computed from the planned boxes before a frame is read. A reflection learning
  into a final block stops the program: that would be a wrong rule.
* **A reflection waits for the last block its interpolation weights reach**,
  taken from `neighbours_of` exactly as the fit uses it.
* **Blocks are finalised in order by `finalise_reference`'s own rule**: a sparse
  cell takes its block's average, an empty block the latest earlier one's, else
  the first later one's. A whole-scan fallback would hold boxes to the end.
* **The first pass changes only a box's background**, to the GLM's mean, which
  the fit sets anyway; so a held box drops its background and is the box the
  second pass would rebuild. If that ever changes -- `close()` touching the mask
  or pixels -- the held box must be copied first.
* **A held box's pixel cells are learning's own**: the geometry of where its
  subdivisions fall in the grid, from its extent, s1, phi and the grid, all the
  same for its fit. If a fit ever used a different s1 or phi or grid -- a
  refined model between learning and fitting, say -- it must compute them.
* **A chunk's save skips the boxes one pass holds**: their release saves them.
  Saving them in the chunk overwrote boxes fitted in the same chunk with the
  emptied ones the move left, and the table failed to write.

## Symmetry and scaling: what has to stay true

Each of these was got wrong once, or had to be found. `docs/symmetry.md` and
`docs/scaling.md` have the measurements.

* **Score symmetry to a resolution limit, with Friedel mates apart.** Without
  the limit, noise beyond the diffraction pulled every element's CC toward zero
  and a cubic sweep came out P -1; with mates merged, the identity compared each
  reflection with itself, gave CC 1, and set the expected CC of every true
  element too high.
* **Candidate space groups are the chiral ones.** A centrosymmetric group is its
  own Patterson group; without the rule I m -3 is offered for I m -3.
* **gemmi names a group only in a tabulated setting.** The change of basis to a
  reference setting is searched for (`reference_setting`); for a group with one
  axis the conventional a and b lie perpendicular to it, where no rotation of
  the group fixes them.
* **The error model's deviations take the exact variance** of a difference from
  a mean the observation is part of, v_i (1 - 2 c_i) + sum c_j^2 v_j -- not
  sqrt((n-1)/n), which biases a low by (n-1)/n. The tests fail with that form.
* **The scale's uncertainty goes into the variance before the error model**,
  I^2 var(g) / g^2, so that a and b correct what remains after it.
* **Differentiate with <I_h> following the parameters** (variable projection):
  held fixed, the gradient is right and the curvature too large, and fits crept
  for fifty steps. With the curvature right, a nearly free constant in B
  wandered, so B is centred on its mean as the scale is on its.
* **One d for completeness**: observed and possible reflections both from the
  crystal at the unique index, or a shell comes out above 100 per cent.
* **Change a column by copying it and setting it back.** `int_column` makes a
  column of zeroes; used on the flags, it wiped every flag integration had set.
* **Parallel work is cut into fixed blocks and combined in block order**
  (`src/parallel.hh`), never into one block a thread: then the answer does not
  depend on the thread count, as integration's does not. And a flag written
  from several threads is a byte, not a `std::vector<bool>` bit, whose
  neighbours share one word and one write.
* **A model reference is an integer.** `json::Value(0.0)` writes a float that
  Python's json reads as one, and a float does not index a list. See also
  "`0` is not `0.0`" above.

## Style

Row-major Mat3, because that is how DIALS serialises a matrix into a flat nine,
and a transpose hidden in the I/O layer is the easiest possible day to lose.

Tests register themselves by static constructor, so there is no list to forget
to add to. Comments explain why, not what.
