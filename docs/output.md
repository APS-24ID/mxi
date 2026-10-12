# What the programs say, and for whom

STATUS: plan, for Graeme's agreement, 12 October 2026. Nothing here is built.

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

## What changes

Each program's printed output for a crystallographer first. What only a
developer needs moves behind `--verbose` (it is still in `mxi_<program>.log`).

* **mxi_import** -- as it is: wavelength, distance, beam centre, panel, scan.
  Add the total rotation, and the exposure where the file gives it.
* **mxi_find** -- the preamble (images' path, identifier, thread count, mask
  count, grouping parameters) to `--verbose`. Keep the counts and the
  histogram of spots per image. Add spots by resolution, and ice rings
  flagged if the strong spots crowd them.
* **mxi_index** -- add the candidate solutions considered, as dials.index's
  table; the cell with e.s.d.s; the refined beam centre and distance beside
  import's, so a wrong beam centre is seen; and per sweep, or per lattice, the
  fraction indexed. With several lattices, a table a lattice.
* **mxi_refine** -- add the cell with e.s.d.s; the refined geometry (distance,
  beam centre, wavelength) and what refinement moved it by; and, scan-varying,
  how far the cell and the orientation wander over the scan, the check the
  tutorial makes in dials.report.
* **mxi_integrate** -- the image path to `--verbose`; progress on a long run
  (images done, the rate, the time left); the summaries kept.
* **mxi_symmetry** -- each symmetry element described as DIALS does, its order
  and axis ("2 |(1,1,0)", "3 |(1,1,1)"), not a letter; the rest kept.
* **mxi_scale** -- kept: the model, the error model, the merging statistics, the
  resolution. Add the scale's range by image (decay at a glance).
* **mxi_export** -- the MTZ's summary: space group, cell, resolution range,
  batches, the columns.
* **Every program** ends by naming what it wrote and the step that takes it
  next ("Next: mxi_refine indexed.rflx"), as a tutorial would.

## For a GUI: each program's results as JSON

Beside `mxi_<program>.log`, each program writes `mxi_<program>.json`: the
numbers its printed tables hold, keyed and typed -- indexing's candidates and
cycles, refinement's steps and models, integration's summaries, symmetry's
scores and choice, scaling's merging statistics by shell. A GUI then reads
data, not text, and the printed layout can change without breaking it; scripts
(a pipeline, a beamline's automation) read the same file. One writer in the
library, each program filling a document as it goes and writing it on exit,
failure included (with the reason).

## The GUI

Port `dials_ui` to mxi as a second back end in the same program, not a copy:
its steps, runner and frame are about running programs in order and showing
their output, and only its step definitions and its parsers are DIALS's. An
mxi back end is a step list (`mxi_import` to `mxi_export`, their options) and a
reader of `mxi_<program>.json` in place of the log parsers. The viewers
(image viewer, reciprocal lattice viewer) stay DIALS's, given the DIALS pair
through `mxi_convert` where they need it.

## Questions for Graeme

1. **JSON beside the log** -- `mxi_<program>.json` -- for the GUI and scripts,
   rather than DIALS-shaped tables for its parsers?
2. **Where the GUI lives** -- a second back end in `dials_ui` itself (one GUI,
   DIALS or mxi), or a port inside this repository?
3. **`--verbose`** for the developer's lines, or kept in the log file only
   (printed nowhere, written to `mxi_<program>.log` always)?
4. **The "Next:" line** -- helpful, or noise for those who know the chain?
