# The documents

What each is for. A **reference** says what a part does, how to run it and how
well it works, and is kept true to the code; **notes** record invariants and
lessons for whoever changes the code; a **history** is kept as it was written,
its conclusions sometimes overturned later, with an index of those.

| document | kind | about |
|---|---|---|
| [`README.md`](../README.md) | reference | building, the chain, the programs, where it stands |
| [`docs/review.md`](review.md) | guide | for reviewers: what this is, how it differs from DIALS, where things are |
| [`docs/algorithms.md`](algorithms.md) | reference | each stage's methods, their sources, where they live and how they differ from DIALS |
| [`docs/outstanding.md`](outstanding.md) | list | every open task, with its evidence |
| [`docs/backstop.md`](backstop.md) | plan | the backstop shadow, its flare, and scaling's outlier rejection: what was found and what is to be done |
| [`docs/multi-sweep.md`](multi-sweep.md) | plan | more than one sweep, from import to scaling: where each program stands, the steps, how each is tested |
| [`docs/multi-lattice.md`](multi-lattice.md) | plan, built | several lattices in one data set, with shared models and their own crystals and profiles; sweeps carried through indexing, as dials/dials#1029 proposes: the model, each program, the stages, the questions |
| [`docs/rflx.md`](rflx.md) | design, built | `.rflx`, one HDF5 file a step in dxtbx-h5's layout, as mxi's primary format: the decisions, the format, the programs, the converter, and how it will be judged |
| [`docs/output.md`](output.md) | built; the GUI later | what each program prints for a crystallographer, the processing's summary in each .rflx, mxi_show, mxi_report's HTML report, and porting dials_ui to mxi |
| [`docs/regression.md`](regression.md) | built; bisecting next | regression runs over real data sets: mxeq regress, its directory, the figures kept and how runs are compared |
| [`docs/correlation.md`](correlation.md) | built | mxi_correlation: clustering data sets by their intensities, as dials.correlation_matrix -- the matrix, the coordinates, OPTICS |
| [`docs/export.md`](export.md) | design | mxi_export, an unmerged MTZ as dials.export writes it: the intensities and their corrections, the columns, the batches and their headers |
| [`docs/import.md`](import.md) | reference | `mxi_import`, and checking it against dials.import |
| [`docs/integration.md`](integration.md) | reference | `mxi_integrate` |
| [`docs/symmetry.md`](symmetry.md) | reference | `mxi_symmetry` |
| [`docs/scaling.md`](scaling.md) | reference | `mxi_scale` |
| [`docs/spots.md`](spots.md) | reference | the spot finder, `mxi_find` |
| [`docs/gpu.md`](gpu.md) | reference | where integration's, indexing's and refinement's time goes, and what a device would take; the refinement target's port, designed and not written |
| [`docs/conventions.md`](conventions.md) | reference | conventions checked against DIALS' files, and what is not |
| [`CLAUDE.md`](../CLAUDE.md) | notes | invariants and lessons across the pipeline |
| [`docs/spots_notes.md`](spots_notes.md) | notes | the spot finder's, from when it was a repository of its own |
| [`docs/spotfinder.md`](spotfinder.md) | history | how the spot finder came into this tree |
| [`docs/integration_history.md`](integration_history.md) | history | integration's development, with its overturned conclusions |
| [`python/README.md`](../python/README.md) | reference | `mxeq`, the comparison tools |
| [`python/CLAUDE.md`](../python/CLAUDE.md) | notes | `mxeq`'s |

Every path and program a document names is checked by
[`python/tests/test_documents.py`](../python/tests/test_documents.py), which fails on a stale one.
