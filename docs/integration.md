# Integration

`mxi_integrate` integrates rotation data by summation and by profile fitting,
and writes a reflection table that DIALS reads -- `dials.scale`,
`dials.export` and `dials.image_viewer` take it as they take their own, as the
DIALS pair: written by `-o integrated.refl --output-expt integrated.expt`, or
from `integrated.rflx`, what it writes otherwise (`docs/rflx.md`), by
`mxi_convert integrated.rflx`.

This is the reference: what it does, how to run it, what it writes, the
conventions it follows, and what is still open. **Why** it is the way it is --
the measurements behind each choice and the mistakes made on the way -- is in
`docs/integration_history.md`, a notebook kept as the work happened.

**Several sweeps.** Given an experiment list of several, `mxi_integrate`
integrates each sweep alone -- a whole run on a list of that one, its
reflections chosen by `id` with their shoeboxes -- and joins the tables, each row
its sweep's `id` and its image set's `imageset_id` in the joined list, as
dials.integrate writes them and dials.image_viewer reads them, and the
experiment lists, each sweep with its own profile model; a crystal they shared is
written once. Every option applies to each sweep, `--postrefine` included:
`--save-shoeboxes` keeps every sweep's shoeboxes, decoded and encoded again as
one column in the joined table, and `--save-profiles F` writes each sweep's
reference profiles to a file of its own, its index before F's extension --
profiles.txt as profiles_0.txt, profiles_1.txt (`docs/multi-sweep.md`).

## How well it works

Measured, each against the thing that can actually judge it:

* **Against DIALS, through scaling.** On 1800 images of insulin, after the
  fitted-variance fix, `dials.scale` gave Rmerge 0.038 and Rpim 0.009 on this
  output and 0.038 and 0.009 on DIALS'. With the module-gap flags corrected,
  the merging statistics on a 3600 image Eiger 16M run are close to DIALS'.
* **Against the data's own symmetry.** A reflection crossing a module gap is
  judged by its clean symmetry equivalents, which involves no other program.
  Where a tenth to a third of a reflection is lost, this output comes to 0.967
  and DIALS' to 0.777, each relative to clean reflections judged the same way
  (1800 insulin, matched inputs, unscaled).
* **Against physics.** Every prediction lies on the Ewald sphere under the
  crystal model at its own predicted position; position uncertainties give a
  pull of one on spots planted at known positions. Both are tests.

Speed, on a 3600 image Eiger 16M sweep, 1.08 million reflections, 16 threads
(the maintainer's machine):

    phase                    before    now
    reading frames           19.6 s   11.0 s
    profile fitting          35.8 s   12.2 s
    total                    67.7 s   36.9 s

## Running it

    mxi_integrate refined.rflx                                    # integrated.rflx
    mxi_integrate refined.expt refined.refl -o integrated.refl    # the DIALS pair

`refined.expt` is a refined experiment, scan-static or scan-varying.
`refined.refl` is optional: given, the profile model is estimated from its
indexed spots; otherwise it is read from the `.expt`'s profile block. It has to
be INDEXED reflections -- the estimate needs each spot's predicted diffracted
beam -- so the spot finder's `strong.refl` will not do, and says so. The images
are found from the `.expt`'s imageset template, or named with `--images`.

Every run writes `integrated.rflx` -- or the pair, `integrated.refl` and
`integrated.expt`, when `-o` names a `.refl` or `--output-expt` is given -- the
reflections and the models integrated with, and the profile model used as DIALS' `gaussian_rs`
block -- what `mxi_symmetry` and anything after it read. Integrating again from
that `.expt` alone takes the profile model from its block, and reproduces the
first integration byte for byte.

`mxi_integrate --help` is the authority on options and defaults. By purpose:

| purpose | options |
| --- | --- |
| what to integrate | `--d-min`, `--d-max` (a low resolution limit), `--first-image`, `--last-image`, `--min-zeta` (0.05) |
| the profile model | `--sigma-b`, `--sigma-m`, `--n-sigma` (3), `--box-scale` (1.9), `--gain` (1) |
| reference profiles | `--scan-blocks` (one per 10 degrees), `--regions` (3), `--reference-signal` (10), `--grid-points` (6; 4 before, `docs/outstanding.md` item 2), `--subdivisions` (5) |
| profile fitting | `--least-measured` (0.6), `--summation-only`, `--two-pass`, `-g`/`--gpu` (single precision on the device, `docs/gpu.md`), `--gpu-emulate` |
| speed and memory | `--threads` (0, one per core), `--window` (64), `--max-boxes` (6000 a thread, at least 20000) |
| output | `-o`, `--save-shoeboxes`, `--save-profiles`, `--save-background-parameters`, `--timing` |

Choosing a few of them:

* **`--scan-blocks`** divides the integrated range, not the whole scan, into
  blocks of their own reference profiles. One block per 100 images or so
  follows a crystal that changes. Too many and a block holds too few strong
  spots and borrows the detector average; the run reports how many did.
* **`--min-zeta`** leaves out reflections close to the rotation axis. Their
  extent in phi goes as 1/|zeta| and the smallest are forty images deep; below
  |zeta| of 0.3 summation agrees with DIALS far worse than above it.
* **`--window`** is frames read per chunk. It bounds memory, not re-reading:
  a frame is read once a pass whatever it is. Reading speed is flat from 32 to
  256; the default keeps the open boxes few.
* **`--postrefine`** integrates, refines against the centres integration
  measured -- scan-static, then scan-varying with the detector held -- and
  integrates again with the refined models, which are then what `--output-expt`
  holds. It uses only reflections that were summed and have a
  centre of mass: an undetected one's observed position is its prediction, and
  a residual of exactly zero pulls the refinement back to where it began.
  `--postrefine-points` sets the scan-varying control points; the default is
  one per 10 degrees, at least five, as `mxi_refine --scan-varying` with no
  number (36 degrees and two more before 30 September 2026). It removes the
  z offset described under open questions, at the cost of integrating twice.
* **`--save-shoeboxes`** keeps pixels and masks, about 31 kB a reflection. A
  msgpack binary cannot exceed 4 GB and the writer refuses rather than wrap; a
  slice with `--first-image` and `--last-image` keeps it under.

## How it works

1. **The profile model.** sigma_b and sigma_m, from the command line, else the
   strong spots, else the `.expt`, in that order. From the strong spots, each
   pixel's count less the background its shoebox records, as Kabsch section 3.1
   step (v) says and DIALS does not: on a table from `dials.find_spots`, or
   `mxi_find` by default, that background is zero and nothing changes; under
   `mxi_find --subtract-background` it is the threshold's local mean, and the
   widths are the spots', not the spots' and their background's. On the 300
   image insulin sweep it took sigma_b from 0.0273 to 0.0270 degrees and left
   sigma_m at 0.1286 -- so background in the spots is not why sigma_m is wider
   than DIALS's.
2. **Prediction.** Every reflection's crossings of the Ewald sphere, once per
   turn on a sweep of more than one, with a scan-varying crystal looked up at
   each reflection's own position in the scan.
3. **Boxes.** A shoebox around each prediction: foreground within `--n-sigma`
   in Kabsch's (eps1, eps2, eps3) frame, the box `--box-scale` times wider for
   the background.
4. **The first pass.** The scan is read in chunks. A box opens with the chunk
   its first frame is in and closes when its last frame is filled, so each
   frame is read once. On closing: a robust Poisson GLM background (Parkhurst
   et al. 2016), summation with Leslie's variance (1999, equation 11), the
   centroid, and -- for strong, nearly whole reflections -- a contribution to
   the reference profile of its detector region and scan block.
   A cell with too few spots borrows its own scan block's average across the
   detector: the profile drifts along the scan, and a block's average is known
   once its reflections have been seen. (It borrowed the whole scan's until
   29 September 2026. On the 300 image sweep with 5 x 5 regions and ten blocks,
   28 of 250 cells borrow; 21 of 21031 fitted intensities moved by more than 0.1
   sigma, none by more than 0.31, and the merging statistics not at all.)
   The scan is divided into one block per 10 degrees of the range integrated
   unless `--scan-blocks` says otherwise: a block should be a fixed angle, since
   that is what the profile drifts over, and in a single pass its length is how
   long a shoebox waits to be fitted. (Five blocks whatever the scan until 29
   September 2026. On the 300 image sweep, 30 degrees, three blocks against five:
   the same Rmeas 0.041 and CC half 0.987, the lowest shell's CC half 0.977
   against 0.979, and the error model 1.009 and 0.0244 against 1.017 and 0.0236.)
5. **Fitting, in the same pass.** A box is held until the scan blocks its
   profile is interpolated from are final, then fitted by weighted least
   squares against its own pixels, with the reference profile interpolated
   between the neighbouring cells and carried onto the pixels, and let go. The variance is
   the fit's plus the background term of Leslie's equation 34, the background
   estimate's error carried through the fit's own weights (below, "The fitted
   variance").
6. **Writing.** One row per prediction.

**Both transforms share one geometry.** Learning carries counts onto the
grid (`transform_shoebox`) and fitting carries the profile back onto the
pixels (`profile_on_pixels`), and both get it from `pixel_cells` and
`plane_weights`: eps1 and eps2 at the pixel corners, interpolated, with any
subdivision near a cell boundary computed exactly, so the cells are the ones
a direct computation gives. The direct versions are kept as
`transform_shoebox_direct` and `profile_on_pixels_direct`, and the tests hold
the fast ones to them.

**The output does not depend on the thread count.** Every parallel reduction
adds in an order fixed by the data, not by the scheduler: profile learning sums
sixteen fixed blocks of reflections, each in index order, and adds them in block
order. One thread and eight give the same bytes, and a test holds it.

Profile fitting is on the pixels, not the grid, because the grid has more
points than a shoebox has pixels and one pixel reaches several of them, so
treating grid points as independent overstates the information. That and a
missing background term together had the fitted variance at 0.43 of DIALS' at
high resolution, and scaling saw it before any comparison did.

A reflection crossing a module gap is fitted to the pixels it has, and its
intensity is the fitted scale times the WHOLE profile, not the visible part --
which is what lets a fitted intensity recover what the gap took.

## Why a reflection is not profile fitted

Scaling takes only profile-fitted observations, so a reflection the fit loses is
lost to the merged data. After its summary the integrator says how many were
lost and why, and writes each reflection's reason to the table as
`profile.failure`, 0 if it was fitted (DIALS ignores the column):

| code | reason |
| --- | --- |
| 1 | the box rejected before fitting |
| 2 | no reference profile for its place and scan block |
| 3 | the profile not carried onto its pixels |
| 4 | the profile zero over the whole foreground: off its box |
| 5 | none of the foreground measured: wholly in a gap, say |
| 6 | under `--least-measured` of the profile measured |
| 7 | the least squares had no solution |
| 8 | the device's fit failed; a run without `--gpu` says which |

`mxeq failures integrated.refl` counts them and bins the failures by position in
the scan, |zeta|, the box's depth and width, place on the detector, resolution
and summed I/sigma, so that a cause shows as the variable whose bins differ. On
the 300 image insulin sweep, 1011 of 21031: 728 wholly in a gap and 283 under
0.60 measured -- the same with `--gpu-emulate`, on the device's batching path,
and with sigma_m or sigma_b forced small.

## The detector's markers

A detector marks two kinds of pixel with the largest values of its width, and
neither is a count: max() - 1 a bad pixel, max() a tile join, the gap between
modules -- 0xfffe and 0xffff in 16 bits, 0xfffffffe and 0xffffffff in 32. The
spot finder has always masked both, on the CPU and the GPU. The integrator took
only max() for a marker until 2 October 2026, so a bad pixel was filled into a
shoebox as a count of 65534: in a reflection's foreground, that much added to
its intensity; in its background, left to the robust fit. Now both are markers
(`src/fill_row.hh`): the voxel left at 0 and its validity cleared, as a tile
join's always was. Data with no bad pixels -- the insulin sweeps here -- integrate
byte for byte as before; data with them integrate differently, and rightly.

## Chunks on dense data, and filling shoeboxes

A chunk's frames are read in parallel, a frame a thread, and a chunk ends where
opening its boxes stops: so `--max-boxes` sets how many threads have a frame to
read. A fixed 20000 cut ferritin's 64 frame chunks to some 11 frames -- 6.2
million reflections over 3600 frames, 1725 opening a frame -- and on 32 threads
reading was a third busy. The default is now 6000 a thread, at least 20000:
192000 on 32 threads, measured there at 200000 to take the integration from
112.0 s to 102.2, threads 90 per cent busy, 2.2 GB more held.

With every thread busy each did its work at half the speed -- 16 cores of two
threads each -- and filling the shoeboxes, two thirds of reading's work, was a
pixel at a time with a branch each. A row of a box is now converted in one loop
with no branch, which the compiler makes SIMD, checked for the detector's markers
by an OR across it, and gone over again only where it has one: on one thread
here 1.41 s to 0.355, the same bytes.

**Not zeroing a box when it opens was tried and reverted.** Every voxel is set
to 0 at opening and then overwritten by the fill, some 140 GB of zeros a run on
dense data, so leaving them unset -- zeroing only slices no frame fills -- looked
free. It was byte-identical and changed nothing on the 32 thread machine
(opening 15.05 s against 15.09), and was slower on a MacBook: opening's cost is
building each voxel's mask, not writing its zero.

**Chunking moves the last bits of the profile fit.** The reference profiles are
sums over boxes in the order they close, which follows where chunks end: with
`--max-boxes 1000` against the default on the 300 image sweep, `intensity.prf`
differs by at most 2.5e-13 relative, the summation not at all. So the new
default changes dense data's profile-fitted intensities at that level, where
the cap binds; on the 300 image sweep it does not bind, and nothing changes.

## One pass over the images

Each frame is read once. A reflection's profile is interpolated from the cells
around it, trilinearly between cell centres, so it needs the scan block whose
centre is at or before it and the next; and a block is complete once every box
learning into it has closed, which is known before a frame is read, from the
predicted boxes. So after each chunk the blocks whose boxes have all closed are
finalised, in order, against the rule `finalise_reference` applies to the whole
scan -- a sparse cell takes its block's own average -- and every held box whose
blocks are final is fitted and dropped. A box is held without its background
array, which the fit restores from the GLM's mean as the second pass did: 5
bytes a voxel. A box learned from keeps the cells its transform computed -- where each
pixel's subdivisions fall in the grid -- and its fit uses them rather than
computing them again.

The answer is the two passes' own, byte for byte. The profiles are learned in
the same chunks and the same order, a finished cell only has zeroes added to it
afterwards, and the first pass never changes a box's pixels or mask -- only its
background, to the same constant -- so a held box is the one the second pass
would rebuild. On the 300 image sweep: identical with the defaults, with 5 x 5
regions and ten blocks and 28 cells borrowing, with one block holding every box
to the end, with shoeboxes saved, and on one thread and four; 300 frames read
where two passes read 600. On a smaller data set on a MacBook, the same SHA-1
both ways and 7.0 s against 8.1. `--two-pass` reads twice, to compare against, and
`python/tests/test_single_pass.py` holds the two to the same bytes. If a reflection ever
learned into a block already final the program stops rather than give a
different answer.

`--timing` reports the most shoeboxes held and their memory: 14060 and 0.32 GB on
the 300 image sweep, most of it, since three 100 image blocks and boxes up to
148 images deep keep nearly everything waiting. On the 360 degree 16M sweep, on a
MacBook with 16 threads: at most 79584 boxes and 2.66 GB, where an estimate from
the 300 image sweep had said 49000 and 1.1 -- boxes wait longer than the estimate
allowed. Integration there took 19.4 s against two passes' 28.4, and decompressing
41 thread-seconds against 80.

## Several lattices, and neighbours masked

**Neighbours.** Each reflection's shoebox has the voxels inside any other
reflection's foreground marked Overlapped -- under that reflection's own
experiment and profile model, so another lattice's spots are masked as they
are -- left out of its background, kept in its foreground and flagged
`overlapped_bg` or `overlapped_fg` (bits 11 and 12, as DIALS's); a reflection
whose foreground is overlapped is not learned from. DIALS's 3d_threaded
integrator does this; its default, 3d, finds the overlaps and masks nothing.
On a single well-separated lattice nothing overlaps and nothing changes; on
two crystals of insulin, a sixth of each lattice's boxes meet another.
`--no-mask-neighbours` leaves them in. Which boxes meet: a hash of tiles of 32 pixels and 16
images, looked up on every thread. Keyed on the panel's pixels alone, a tile
held every reflection that fell there through the whole scan: on a ten-turn
scan of 2.5 million boxes the search took 178 s on one thread, a couple of
minutes before Graeme's 7.4 million began integrating; now 4 s on two.

**Lattices.** The experiments of one sweep -- on one image set and one scan --
are integrated in one pass over its images, each lattice with its own profile
model, estimated from its own spots, its own profile grid and reference
profiles, finalised apart, so that one lattice never borrows another's
average. Each experiment's profile model is written as its own profile block.
The same, row for row, as a run a lattice (`--each-alone`), but for the row
numbers and rounding in the fitted values. Several sweeps are integrated a run
a sweep, each its lattices in one pass and its own crystal-less experiment
carried, the runs joined in the list's order; postrefinement and
`--each-alone` a run an experiment. `docs/multi-lattice.md`.

## What it writes

Every column DIALS' own `integrated.refl` has that this computes:
`miller_index`, `entering`, `panel`, `id`, `imageset_id`, `partial_id`, `bbox`,
`flags`, `d`, `s1`, `zeta`, `lp`, `qe`, `partiality`,
`intensity.sum.value/variance`, `intensity.prf.value/variance`,
`profile.correlation`, `background.mean`, `background.sum.value/variance`,
`num_pixels.foreground/background/background_used/valid`, `xyzcal.px/mm`,
`xyzobs.px.value/variance` and `xyzobs.mm.value/variance`. `lp`, `qe`, `d` and
`partiality` were each checked against DIALS' values for the same reflections.

Columns of this package's own:

* **`intensity.sum.variance_background`** and
  **`intensity.prf.variance_background`** -- each variance with no signal: the
  background's counts and its estimate's term; the fit's with weights b
  alone. Scaling builds each observation's variance on them and the model's
  intensity ("The fitted variance"; `docs/scaling.md`, "Variances from the
  model").

* **`xyzres.px.value`** and **`.variance`** -- the observed centre of mass less
  the prediction, in fast pixels, slow pixels and images. Both ends are in the
  frame the images are in, parallax included, so it is a plain difference. The
  centre is the UNCLIPPED one, which is not the centre `xyzobs.px` holds; its
  uncertainty is the one checked against spots at known positions. NaN where
  the reflection was not detected at three sigma, or its centre fell outside
  its own box. See the open questions for what it cannot tell you in z.
* **`profile.measured`** -- the fraction of the reflection's profile the
  detector recorded; one unless it crosses a masked pixel.
* **`profile.failure`** -- why it was not profile fitted, 0 if it was; the codes
  are under "Why a reflection is not profile fitted".
* **`background.dispersion`** and **`background.dispersion_trimmed`**, with
  **`num_pixels.background_trimmed`**, under `--save-background-parameters`
  only, so that by default the table is what it was -- the background pixels'
  sample variance over their sample mean: about 1 on a flat background, Poisson
  on a photon counter, far more on a ramp such as a backstop shadow's edge. Over
  every valid background pixel, and trimmed of those beyond five Poisson
  standard deviations of the robust mean (a neighbour's spot, a hot pixel), with
  how many that leaves. NaN where it cannot be said. Measured; nothing is
  decided by it yet (`docs/backstop.md`).

### Flags

Named and numbered as in DIALS' Flags enum,
`dials/array_family/reflection_table.h`:

| flag | bit | set when |
| --- | --- | --- |
| `predicted` | 0 | always |
| `integrated_sum` | 8 | summation succeeded AND the foreground reached no masked pixel |
| `integrated_prf` | 9 | profile fitting succeeded and `profile.measured` is at least `--least-measured` |
| `foreground_includes_bad_pixels` | 14 | the foreground reached a masked pixel |
| `background_includes_bad_pixels` | 15 | the background did |
| `failed_during_summation` | 19 | with bit 14: the sum is missing part of the reflection |

A reflection whose foreground crosses a module gap is not flagged as summed,
because its sum is not its intensity. That is DIALS' convention, and it
matters: `dials.scale`'s default intensity needs `integrated_sum` and
`integrated_prf` together, so these reflections are left out of scaling. To
use their recovered profile-fitted intensities, run
`dials.scale intensity_choice=profile`.

### Saved shoeboxes

In DIALS' mask convention: `valid` 1, `background` 2, `foreground` 4, and a
voxel with no measurement is zero. Internally a bad pixel keeps its region bit
and loses only `valid`, which is what the fit needs to know a part is missing;
it is translated on writing, because DIALS never sets a region bit on an
invalid voxel and rejects a table that does.

## What it prints

Standard output is for results; warnings and errors go to standard error, and
both go to `mxi_integrate.log` in the working directory, in order. A run says
what profile model it used and where the images came from, how many
reflections it predicted and why any were not integrated, and how many
reference profiles it learned; then two tables, modelled on the summaries
`dials.integrate` prints and limited to what is measured here:

* **against resolution**, in shells of equal volume in reciprocal space:
  counts fully and partially recorded, summed and fitted; mean background;
  mean I/sigma by summation and by profile fitting; mean profile correlation;
  and RMSD XY;
* **overall, lowest shell and highest**: the same, with the counts whose
  foreground or background reached a masked pixel, and the Pearson and
  Spearman correlations between the summed and fitted intensities.

There are no overload or ice-ring columns, because nothing here detects either.
Thread counts, chunking and frames read are under `--timing`, with the phases.

Two numbers are defined so as to be read beside DIALS', and one is not:

* **RMSD XY** is the observed centre, `xyzobs.px`, against the prediction, over
  reflections integrated by summation and detected. On a 300 image insulin
  sweep it is 0.267 pixels where DIALS gives about 0.25. It is not built from
  `xyzres.px`, the unclipped centre, which is honest about its noise and so twice
  as scattered for weak spots -- 0.544 on the same data. That is the one to
  judge predicted positions by, and the wrong one to set beside DIALS'.
* **Partially recorded** counts reflections whose partiality is below 0.99. It is
  not DIALS' count, which is several times larger: DIALS splits a reflection
  crossing one of its job boundaries into separate partial rows, and counts
  reflections running off the ends of the scan.

## Conventions that caused bugs

Each of these was once wrong here. They are the ones to check first when
something disagrees.

* **Validity and region are separate bits.** A pixel is measured if and only
  if `mask & valid`. Testing `mask == 0` stopped meaning "not measured" the
  moment a bad pixel kept its region bit, and five places broke at once.
* **z is the array index,** image n at n - 1. A frame arrives numbered by its
  index in the file; `single_file_indices` in the `.expt` says which file index
  each scan image is.
* **A crystal's scan points are samples in a file and control points in a
  refinement.** A file holds A at the image boundaries, N + 1 for N images, and
  they are interpolated linearly, as DIALS does; a scan-varying refinement moves
  the control points of a spline and writes it back as N + 1 samples. Mixing
  the two once had predictions from any file re-smoothed, and a model refined
  twice written with five scan points for three hundred images.
* **An angle must be in the scan before the crystal is looked up at it.** The
  Ewald solver answers in whatever 2 pi interval its arithmetic lands in; a
  crystal looked up at -260 degrees on a scan from 0 to 180 is the crystal from
  the start of the scan.
* **Flags come from DIALS' source,** not from memory or from a table here:
  a guess once took bit 19 for an exclusion flag.
* **Fit on the pixels, not the grid,** and give the fitted variance its
  background term -- carried through the fit's weights, not as a sum carries
  it. That the fitted variance could not fall below about half the summed was
  wrong: equivalents' deviations say a weak reflection's is a quarter
  ("The fitted variance").

## The fitted variance

Judged by the deviations between equivalents -- one for sigma that are right
-- the weak reflections' fitted variance was too large, and so were their
sigma after scaling, in mxi's pipeline and DIALS's own (`docs/scaling.md`).
Integrating insulin's sweep with a crystal turned 3 degrees, every prediction
on bare background, the summed I/sigma had variance 1.12 to 1.25: summation's
sigma right, or a little small. Two things in the fit were not:

* **The background term.** The background's error is one number under every
  pixel; through the fit it moves the intensity by profile_sum sum(p/v) /
  sum(p^2/v), not by the m foreground pixels a sum carries it by. Taking
  summation's (m/n) I_bg overstated it by the square of their ratio. Carried
  through the weights, insulin's weak fitted deviations go from 0.58 and 0.71
  to 1.04 and 1.06, summation's being 1.04 and 1.02; the weakest fitted
  variance is a quarter of the summed, which is what the deviations say it
  achieves; scaled, the weakest sextile goes from 0.66 to 0.95.
* **The least variance of a pixel, 0.05 counts.** Two crystals of insulin
  (`docs/multi-lattice.md`) have a background of 0.016 to 0.034 counts a
  pixel, below it everywhere, so the floor, not the counts, set a weak
  reflection's weights and variance. Now 0.001, below any background met:
  insulin's, about 0.3, changes 17 variances of 20018. Scaled, the two
  crystals' weakest sextiles go from 0.26 and 0.58 to 0.42 and 0.80, the rest
  near one.

**The summed variance from counts expected, not counted.** Graeme's ten-turn
insulin, multiplicity 210: merged I/sigma -5.7 in the outer shell. Leslie's
equation 11 was taken with the counted foreground for I + I_bg; at a hundredth
of a count a pixel an empty foreground has I = -I_bg and next to no variance,
one count I > 0 and a variance of one, so each sigma moved with its own noise.
On the two crystals' bare background (a crystal turned 3 degrees): mean I
+0.011 counts, unbiased, but mean I/sigma -0.35, -0.53 at the edge, and the
inverse-variance weighted mean of 210 at -12.9 sigma -- merging's negative
intensities. With the expected counts, max(I, 0) + I_bg: mean I/sigma -0.07,
the weighted mean of 210 -1.9 sigma. The rest is the max(I, 0), a positive
fluctuation still raising its own variance; only a model's intensity, the
same for every observation of a reflection, removes it -- in scaling. DIALS
(summation.h) takes |I| + I_bg (1 + m/n): noise-dependent too, but evenly either
side of zero, so it does not lean a mean low as the counted form did. Insulin's
sweep: 877 variances of 21030 change, merged I/sigma 18.9 to 18.7.

**Variances from the model, in scaling.** Graeme's rerun was unmoved, digit for
digit: its weak intensities are profile fitted, and the fit's weights b +
max(I, 0) p make its variance move with its own noise the same way --
simulated at 0.005 to 0.5 counts a pixel, a weighted mean of 210 at -3 to
-18.6 sigma, worst for weak but real signal. No variance from an observation's
own pixels is free of it. So integration writes each observation's variance
with no signal -- `intensity.sum.variance_background`, gain I_bg and the
background estimate's term; `intensity.prf.variance_background`, the fit's
with weights b alone, gain b P^2 / sum(p^2) and the background's term carried
by P sum(p) / sum(p^2) -- and scaling (`docs/scaling.md`) builds each
variance on it and the model's intensity.

What remains is the weakest sextile of the weakest crystal: shoeboxes of a
handful of counts, where a variance taken from the fit's own noisy intensity
-- b + max(I, 0) p -- is least to be trusted. Taking it linear in I, which is
right on average, moved the weakest little and left the next too small, the
weights unsettled by negative intensities: not done. The device's fit
(fit_device.hh, fit_cuda.cu, fit_metal.metal) carries both changes; the
emulation matches the CPU, the CUDA and Metal kernels are not built here.

## Checking a run

`mxeq` compares two integrations and explains the differences:

| command | what for |
| --- | --- |
| `mxeq trend A B --value col` | where two runs disagree, binned against resolution, position, frame, intensity |
| `mxeq html A B -o report.html` | the same, drawn |
| `mxeq explain A B` | why observations went unpartnered; run it first when the match rate is low |
| `mxeq disagree A B -o out.refl` | only the reflections that disagree, for `dials.image_viewer` |
| `mxeq residuals A -o r.html` | how well positions were predicted, from `xyzres.px` |
| `mxeq profiles profiles.txt -o p.png` | the learned reference profiles, from `--save-profiles` |

Matching pairs observations, not reflections: Miller index and entering flag
group them and the frame separates them, since a sweep of several turns sees
each reflection many times. The `matched` count in every report should be close
to the smaller table's rows; if it is not, nothing downstream of it can be
trusted, and `mxeq explain` says why.

Before comparing, check that the model and the images are the same dataset.
Integrating one sweep's images with another's refined model gives numbers that
look entirely reasonable.

## Open questions

The whole list of open work, integration's among it, is `docs/outstanding.md`.

What is known not to be right yet, with what is known about each.

* **Gap-crossing reflections are 3 per cent low** where a tenth to a third is
  lost, judged by their clean equivalents. The likely cause is profile
  learning: a reflection that is itself truncated can still contribute to a
  reference profile.
* **A z residual oscillates along the scan and sits below zero** by about 0.1
  images, with a period of 180 degrees of rotation, when the models were
  refined against the spot finder's centres. Those depend on strength -- a weak
  spot shows the finder only its peak, a strong one its whole rocking curve,
  which has a tail toward earlier images -- so the fitted model sits between.
  `--postrefine` refits against the integrator's centres: on a 300 image sweep
  the offset over I/sigma of ten or more goes from -0.115 images to -0.006, flat
  along the scan. What remains is a dependence on strength, +0.018 images for
  the weakest to -0.036 for the strongest, which is the integrator's own centre
  on an asymmetric rocking curve. Whether post-refinement should be the default,
  and what it does to scaled intensities, is open.
* **For a thin spot, `xyzres.px` in z reports the foreground window,** which is
  placed on the prediction -- so it understates the prediction error in z, and
  most for the spots that ought to measure it best. On the detector it is
  sound.
* **The background is about one per cent high at low resolution.** The signal
  reaches beyond the foreground -- measured falling from 0.74 counts a pixel at
  the foreground's edge to 0.58 at two and a half radii -- and a guard ring at
  1.5 radii would exclude it for about fourteen per cent of the background
  pixels. Not built.
* **Overlapping reflections are masked, not resolved** (Leslie sections 6.3,
  6.7.1): each reflection's neighbours' foregrounds are kept out of its
  background and flagged ("Several lattices, and neighbours masked"); a
  foreground two reflections share is summed and fitted as one.
