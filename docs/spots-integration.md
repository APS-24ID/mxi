# Integrating the spot finder into mxi

STATUS: commit 1, steps 2a and 2b done; 2b-GPU written, to be built on the GPU machines; 16 October 2026.

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

## Commit 2: the shared functions (in steps)

**2a (done): mxi_find as every program.** `spot_refl.*` gone: the spots are a
library table (`spot_table.*`), written by the library's writers -- which now
write every column in name order, shoeboxes among them, as DIALS does, so that
the .refl is byte for byte what it was. `spot_expt.*` is the library's
`experiment_images.*`, still tolerant of stills (the full reader refuses a
scan with no oscillation). The shared argument parser: an experiment list or a
master positionally, `-e` and `-gpu` gone. Checked: the .refl byte for byte,
both .rflx every column, the golden test through the library, ctest, the
Python suite.

**2b (done): the HDF5 handles; the queue and the mask kept.** One handle class
for the library (`hdf5_handle.*`, with the lock around HDF5's API moved there
from the reader), where the reader, the NXmx import and the .rflx code had one
each -- Handle, H5 and H, the reader's alone closing under the lock. Not
changed, on reading them: the spot finder's work queue is not a copy of
`parallel.hh` but a bounded producer and consumer, frames taken as they are
written -- finding spots during collection -- which a known count of work up
front cannot do; it stays, into `mxi` with 2c. The detector's pixel mask is
read once by the reader and applied by `apply_pixel_mask`, which integration
already uses: nothing duplicated. Checked: the .refl byte for byte, the two
spot tables and insulin's integration every column, ctest, the Python suite.

**2b-GPU (written, for Graeme to build): one Metal layer and one CUDA layer.**
`metal_device.*` (library `mxi_metal`): the device, created once for the
process; an embedded shader library loaded; a pipeline with its threads per
threadgroup checked; and metal-cpp's implementation, compiled there alone.
`cuda_device.*`: whether there is a device, why not (the reason had been the
same forty lines twice), its name. The threshold and profile fitting each
keep their own error policy -- the threshold throws, fitting falls back to
the CPU. Neither toolchain is here: the CUDA layer syntax-checked against a
stand-in for the runtime's header, the Metal changes read line by line.

**2c, after:**



* `spot_expt.*` goes: `mxi_find` reads the library's experiment list.
* `spot_refl.*` goes: `mxi_find` writes through the library's writers, `.rflx`
  by default and DIALS's `.refl` on request; the golden test kept against them.
* `mxi_find`'s own argument parsing goes: `parse_arguments`, as every program;
  `mxi_find imported.rflx`, a master given directly kept as a convenience, `-e`
  dropped.
* One Metal layer and one CUDA layer -- the device, its library, its errors --
  shared by the threshold's kernels and profile fitting's; `--gpu` and
  `--gpu-force` alike in `mxi_find` and `mxi_integrate`.
* The namespaces (`series`, `decompress`, `dials_spots`, `gpu`, ...) into `mxi`.
* `tests/spots/` into `mxi_tests`, the `TEST` framework.

## Commit 3: the documents and the framing

* `docs/spots.md` rewritten as mxi's spot finding -- the method, its source, its
  implementations -- not DIALS's on a GPU; `docs/spotfinder.md` and
  `docs/spots_notes.md` folded into it or retired.
* No `dials.import` anywhere: `mxi_find` takes `mxi_import`'s output.
* mxeq's `spotfinder_*` checks renamed.
* This document removed.
