# Integrating the spot finder into mxi

STATUS: commit 1 of 3 done, 16 October 2026.

**Remove this document when the series is complete:** it is a plan, and once
done it belongs to the project's history (git log), not to its documents.

The spot finder was a project of its own -- DIALS's spot finding on a Metal GPU,
where it began -- and lived in `src/spots/` with its own build, options, tests,
namespaces and conventions. This series makes it one stage of mxi like the
others: one build, one source tree, one set of conventions, its common
functions shared with the rest.

Decisions (Graeme): the options renamed MXI_*, the old names taken with a
warning through 1.0.x; a flat `src/`; `mxi_find` run as every program is run;
the standalone spot-finder build dropped.

## Commit 1: structure (done)

Moves and renames, the spots unchanged to the bit:

* `src/spots/` is gone; its files are in `src/`: the image layer
  (`nxmx_reader.cc`, `image_series.hh`, `decompress.*`), the threshold
  (`dext*`, `signal_order.*`, `signal_pixel.hh`, `histogram.hh`, `queue.hh`),
  the spots' shapes and centroids (`spots.*`, was `dials_spots.*`), and for
  now `spot_expt.*` and `spot_refl.*` (were `expt.*` and `refl.*`, whose names
  the library's own files have).
* `mxi_find` is `apps/mxi_find.cc`, as every program.
* The build: `spotfinder_*` targets named `mxi_*`; every `SPOTFINDER_*` option
  and macro `MXI_*`, the old option names taken with a deprecation warning,
  once, through 1.0.x; the two sanitizer settings one, `MXI_SANITIZE` (ON, or a
  list: address, undefined, thread).
* Checked: strong spots on insulin and the two-sweep set, every column the same
  as before; the DIALS `.refl` byte for byte; integration on insulin the same;
  a fresh configure and build with no warnings; every test.

## Commit 2: the shared functions (next)

* `spot_expt.*` goes: `mxi_find` reads the library's experiment list.
* `spot_refl.*` goes: `mxi_find` writes through the library's writers, `.rflx`
  by default and DIALS's `.refl` on request; the golden test kept against them.
* `mxi_find`'s own argument parsing goes: `parse_arguments`, as every program;
  `mxi_find imported.rflx`, a master given directly kept as a convenience, `-e`
  dropped.
* The work queue goes: `parallel.hh`, a reader a worker.
* One Metal layer and one CUDA layer -- the device, its library, its errors --
  shared by the threshold's kernels and profile fitting's; `--gpu` and
  `--gpu-force` alike in `mxi_find` and `mxi_integrate`.
* One HDF5 handle wrapper; the detector's pixel mask shared with integration's
  masking where they duplicate.
* The namespaces (`series`, `decompress`, `dials_spots`, `gpu`, ...) into `mxi`.
* `tests/spots/` into `mxi_tests`, the `TEST` framework.

## Commit 3: the documents and the framing

* `docs/spots.md` rewritten as mxi's spot finding -- the method, its source, its
  implementations -- not DIALS's on a GPU; `docs/spotfinder.md` and
  `docs/spots_notes.md` folded into it or retired.
* No `dials.import` anywhere: `mxi_find` takes `mxi_import`'s output.
* mxeq's `spotfinder_*` checks renamed.
* This document removed.
