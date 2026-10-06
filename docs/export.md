# Export: an unmerged MTZ

STATUS: built, 6 October 2026, and its output understood by AIMLESS (Graeme).
Modelled on
dials.export's MTZ writer (`dials/util/export_mtz.py`, with
`filter_reflections.py`, `batch_handling.py` and the C++ `ub_to_mosflm_u`),
read for this; each rule below is that code's.

`mxi_export EXPT REFL -o scaled.mtz` writes the unmerged MTZ that AIMLESS,
POINTLESS and the rest of CCP4 read: one row an observation, a batch header an
image. One sweep or several. Written through gemmi's C++ `Mtz`, as dials.export
writes through gemmi's Python one -- `third_party/gemmi`'s `mtz.cpp`, `gz.cpp`
and `sprintf.cpp`, with zlib.

## Which intensities, corrected how

As dials.export's `intensity = auto`:

* **A scaled table** -- `intensity.scale.value` present -- exports the scaled
  intensity: `intensity.scale.value` and its variance divided by the inverse
  scale and its square, as I and SIGI, with SCALEUSED and SIGSCALEUSED the
  inverse scale and its sigma. Rows bad for scaling, or with an inverse scale
  not above zero, are left out.
* **An integrated table** exports profile-fitted and summation intensities
  together, only rows integrated both ways: IPR, SIGIPR from profile fitting,
  corrected by LP / QE; I, SIGI from summation, corrected by LP / QE /
  partiality -- a profile-fitted intensity is already the whole reflection, a
  summed one is the part recorded. Rows with QE or partiality not above zero
  are left out.

Then, as dials.export does: rows whose variance is not above zero left out;
partiality below 0.4 left out (`--partiality-threshold`); I/sigma below -5 left
out after the corrections (`--min-isigi`); `--d-min` if asked. Partials are not
combined: mxi integrates whole reflections, each row its own `partial_id`.

## The columns

H, K, L, M/ISYM, BATCH; then the intensities above; BG and SIGBG where the
table has `background.sum.value`; FRACTIONCALC, the partiality (1 where there is
none); XDET, YDET, the predicted position in pixels; ROT, the rotation angle at
the predicted frame, from the reflection's own sweep's scan; LP; QE (1 where
there is none). Indices then put in the asymmetric unit, M/ISYM recording the
operation, and the rows sorted on H, K, L, M/ISYM, BATCH.

## Batches

A reflection's batch is its observed frame, floor(z) + 1, plus its sweep's
batch offset. The offsets as dials' `_calculate_batch_offsets`: each sweep keeps
its own image numbers as batch numbers unless they would overlap another's --
those moved to the next number ending in 01 past the highest so far; sweeps
that only meet, 1 to 150 and 151 to 300, keep theirs -- and batch zero never
used.

## Batch headers

One a image, as dials.export writes them. First every sweep's models are turned
into the Cambridge frame -- the rotation axis along +Z, the beam along +X --
by rstbx's `align_reference_frame`: the rotation taking the axis onto Z, then
the one about Z bringing the beam as near X as it goes. Then for each image:

* **the cell and U**: U the Mosflm U, UB B^-1 with B Busing and Levy's from the
  reciprocal cell, of UB = S F A in the Cambridge frame. For a scan-varying
  crystal, A at the frame's centre: B the mean of the frame's two scan points',
  U the first's turned half way to the second's. (dials.export composes this
  with R(axis datum, phi) and its inverse about the turned axis, which cancel
  to S F whatever phi is -- just as well, as it passes phi in degrees where
  radians are wanted. mxi writes S F U B.)
* **phi start, end and range** of the image, from its sweep's scan;
* **the scan axis**, S times the axis datum, as scanax and e1;
* **the beam**: the sample-to-source direction, and its idealised form, -1 on
  its most negative component;
* **the detector**: the first panel's directed distance and its size in pixels;
* the fixed flags dials.export sets -- one crystal, a 3D batch, one goniometer
  axis named AXIS, the scan axis the first, one detector, batch scale 1 -- the
  wavelength, and the dataset.

The file's cell is the median over the sweeps of each cell parameter, as
dials' `determine_best_unit_cell`; it gives the d-spacings `--d-min` cuts on.
One dataset a wavelength, sweeps within 1e-4 A one wavelength, the dataset's
the mean; crystal XTAL, project mxi (`--crystal-name`, `--project-name`). Every
sweep must have the same space group.

## A pitfall, met

gemmi puts indices in the asymmetric unit, and sets M/ISYM, only for data it
has been told hold their original indices: `switch_to_original_hkl()` called
before the data are set -- which, with no data yet, only marks it, as
dials.export does -- or `switch_to_asu_hkl()` returns false and does nothing,
leaving every M/ISYM zero. mxi_export now marks the data and refuses if the
switch reports failure.

## How it is judged

`python/tests/test_export.py`, read back with gemmi in Python -- scaled: every
row's I and SIGI the scaled intensity and sigma over the inverse scale, its
SCALEUSED, XDET, YDET and FRACTIONCALC the table's, every row's batch a header,
each header's U orthonormal and the scan axis on +Z; integrated: IPR and SIGIPR
by LP / QE, I and SIGI by LP / QE / partiality; and two sweeps overlapping in
their image numbers, the second's batches moved to 201. And the test that
counts: AIMLESS reads mxi_export's output and understands it (Graeme). A
comparison with dials.export column by column and header by header has not
been made, and is there if a difference ever needs explaining. As designed: the columns and their values against the table,
corrections and all; the batch numbering across sweeps against dials' rules;
each header's U and cell against the experiment's models, by the same
construction independently; the frame -- axis on Z, beam on X. And with DIALS
at hand, against dials.export on the same data, column by column and header by
header: the test that it is what CCP4 is used to.
