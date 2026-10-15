# What the programs say, and for whom

STATUS: built, 12 and 13 October 2026, as Graeme decided (below); the GUI
deferred.

Graeme: the programs' output is technical -- where the data are -- and light on
crystallography, refinement's progress and the like. DIALS's logs are mostly of
no use either, but some of each is what the tutorials send people to; and
`dials_ui` (github.com/graeme-winter/dials_ui), a GUI over the DIALS programs,
could be ported to mxi.

## What was read

**DIALS's logs** for a 1800 image thaumatin sweep (`dials.import` to
`dials.scale`), 6000 lines between them. The useful part is a few tables a
program, among configuration echoed back and per-image counts: indexing's
candidate solutions, its crystal model and the indexed fraction, each
macrocycle's residuals, outliers and refinement steps; refinement's RMSDs and
the cell with e.s.d.s; integration's profile model (sigma_b, sigma_m) and its
summaries by resolution and overall; symmetry's elements, scored with their
order and axis, the subgroups and the choice; scaling's model, its error model,
the merging statistics and the resolution limits.

**The tutorials** (dials_tutorials, slac-2026/WORKFLOW.md, the GUI workflow)
point people at: from import, whether the wavelength, distance, beam centre
and image count look right; from spot finding, its progress and spots per
image; from indexing, the cell and the fraction indexed, "ideally close to
100%"; from refinement, the RMSDs against indexing's, and the cell and
orientation varying a little over the scan; from symmetry, the individual
elements and the subgroups' likelihoods, and the conclusion; from scaling, the
merging statistics, the error model and the resolution; from merging, the
MTZ's summary and its columns.

**`dials_ui`** builds each step's command line, runs the real program, streams
its output, and reads `dials.<program>.log` back with regular expressions --
table borders, headings like "Refinement steps" and "Merging statistics by
resolution bin" -- for a Summary tab and live Plots (`dialsgui/parsers.py`,
860 lines; the steps in `dialsgui/steps.py`).

**mxi's own output**, the insulin sweep from import to export (260 lines in
all). Better than feared where it matters -- indexing's and refinement's
cycles with their RMSDs and rejections, the cell, integration's summaries by
resolution and overall, symmetry's subgroups, scaling's merging statistics
and resolution -- but uneven.

## Graeme's decisions

* The tables a program prints go into the .rflx it writes, under a top-level
  `summary` group -- less clutter than a JSON file beside it, and the file says
  how it was made. Each program carries its input's steps and adds its own.
* Developer detail under `--verbose`, to the terminal and the log alike;
  without it the crystallography alone, to both (`cat mxi_find.log`).
* No hint of the next step.
* `mxi_show state.rflx`: the state of processing, summary statistics with it.
* The GUI later.

## What was built

**The summary** (`src/processing.cc`): `{"steps": [...]}` at `/summary`, in the
encoding the experiment list has; one step a program run, `program` and `time`
first, then its figures and tables, a table column by column. A program reads
its input's steps and writes them before its own; `mxi_convert` carries them
between .rflx files and adds none; the DIALS pair has no place for them. On
insulin, from import to scaling: seven steps. `mxeq.rflx.summary_dict` reads
it in Python.

**mxi_show FILE [REFL]**: each sweep's images, rotation, wavelength, distance
and beam centre once, each experiment on it -- its cell and space group, or the
sweep's own crystal-less one -- and its profile model; the reflections by
experiment, strong, indexed, summed, fitted, scaled, outliers; every step of
the summary in order.

**Each program:**

* **mxi_import** -- the total rotation; the first step, each sweep as imported.
* **mxi_find** -- "Finding spots on N images" where the images' path, the
  experiments, the series, threads, mask and grouping were (now `--verbose`);
  spots by resolution, ten shells of equal reciprocal volume, from the
  experiment's geometry; the ice rings the spots crowd, three times a uniform
  share and twenty spots or more; its step, with the spots per image.
* **mxi_index** -- the candidate solutions the basis search weighed, its ten
  best distinct reduced cells, the chosen one marked; the beam centre and
  distance as given and after indexing; the FFT grid under `--verbose`.
* **mxi_refine** -- the cell (its e.s.d.s are scaling's, from two-theta, since
  14 October: `docs/scaling.md`); scan-varying,
  the cell's range over the scan and the orientation's turn; the geometry as
  given and refined.
* **mxi_integrate** -- progress each tenth of the way on runs of 100000
  reflections or more (`--progress-from`); the image path, neighbours and
  reference profiles' bookkeeping under `--verbose`; several sweeps, their
  summaries as one after each sweep's own.
* **mxi_symmetry** -- each element by its order and axis as cctbx describes it,
  2 |(1,1,0), beside its letter, the likely ones marked.
* **mxi_scale** -- the mean scale over each tenth of each sweep, on a line.
* **mxi_export** -- the MTZ's summary: space group, cell, resolution, batches,
  every column with its type.

**mxi_report FILE.rflx [-o REPORT.html]** (13 October): the processing a file
records as a page, after dials.scale's report (Graeme's example), the file
named in its title. The steps run and when, the experiments; then each step's
section from its summary, a plot a sweep wherever the data are a sweep's:

* spot finding -- a row a sweep (images, signal pixels, spots), spots per
  image a plot each, spots by resolution in bins equal in 1/d^2, ticked in d;
* indexing -- the candidates, the RMSDs by cycle;
* refinement -- the cell with e.s.d.s, the RMSDs cycle by cycle, and each
  crystal's cell over the scan, lengths beside angles, as refinement left it
  (its step records every scan point: the file's own crystals are by then
  constrained by symmetry, a cubic cell's angles 90.000);
* integration -- the profile model; a sweep at a time, <I/sigma> fitted and
  summed beside profile correlation and background by resolution, and a row
  each overall; then every sweep together;
* symmetry -- the elements and subgroups, the choice;
* scaling -- the summary laid out as DIALS lays it; CC half, <I/sigma>,
  completeness (0 to 100 per cent, multiplicity on its own axis from zero)
  and the R factors (log scale) by resolution; and from the reflections the
  inverse scale and <I/sigma> by image, each sweep from its first image.

A point a degree along every scan (14 October; Graeme: fifty blocks a sweep
drew nothing on 3600 degrees): spots per image as each degree's mean, the cell
over the scan sampled a degree apart, scaling's inverse scale and <I/sigma>
binned a degree. And the scale model the scaled file holds (mxi_physical):
the scale C and the relative B over each sweep, a point a degree, and each
lattice's absorption surface, S over the directions of the crystal frame, five
degrees a step.

Plotly from its CDN (one file, to mail; it draws with a network), plot text in
characters, not entities, which plotly does not know; an axis spans at least
what makes a constant flat (0.01 A, 0.1 degree). A step it has no section for
gets its tables; long or nested lists are left to the plots. Not yet: the error
model's normal probability plot, anomalous plots, Wilson plots -- dials.scale's analysis plots past the merging statistics.

## The GUI, later

Port `dials_ui` to mxi as a second back end in the same program, not a copy:
its steps, runner and frame are about running programs in order and showing
their output, and only its step definitions and its parsers are DIALS's. An
mxi back end is a step list (`mxi_import` to `mxi_export`, their options) and a
reader of `mxi_<program>.json` in place of the log parsers. The viewers
(image viewer, reciprocal lattice viewer) stay DIALS's, given the DIALS pair
through `mxi_convert` where they need it.

## --timing

Every program's `--timing` table gives each phase's wall time, its share, and
the process's peak memory as the phase ended (`src/timing.hh`); the peak only
rises, so the phase after which it jumps is the one that asked for it.
