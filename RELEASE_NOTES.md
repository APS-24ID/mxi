# Release notes

Newest first. Each release lists first what changes results or commands -- a
default changed, an option renamed or gone -- then what is new, then what is
fixed. The tag of each (`git tag -n99 v1.0.5`) carries the same summary.

## Unreleased -- since 1.1.4

**Changed**
- `mxi view` is mxi ui's image viewer, in a window of its own, no longer
  napari's: the same command and options (`-o`, `--sum`, `--sweep`), and Save
  writing `masked.rflx` and `mxi_mask.json` where it is run, as before. It
  needs Python with tkinter, numpy, h5py and hdf5plugin, as mxi ui does, and
  no napari; the wheel's `[napari]` extra is gone. Regions are rectangles,
  circles and polygons, drawn, chosen and deleted, undone; what napari's
  shapes layer had besides -- ellipses, freehand polygons, moving a region or
  its vertices -- is gone: a region is deleted and drawn again.

**New**
- mxi ui's viewer draws the detector's untrusted regions: rectangles, circles
  and polygons on a sweep's images, chosen and deleted, undone a step at a
  time, the same as the previous sweep's; saved as `mxi view` saves them, into
  `masked.rflx` through `mxi untrusted` and `mxi_mask.json`, and Find spots
  then reads `masked.rflx`. Closing with regions not saved asks first.
- mxi ui's viewer sums images: from the one shown, as many as are asked for,
  so that a shadow or faint spots show; the spots on any of them marked, the
  contrast set afresh.
- `mxi_tools` finds mxi by the `MXI` environment variable first, as mxi ui
  does.

## 1.1.4 -- 4 November 2026

**Changed**
- A cluster `mxi correlation --write-clusters` writes carries its own
  crystals' history, not the whole run's: each step's record cut to its sweeps
  and experiments, the whole run's figures left out and named, scaling's
  record left out (a cluster is scaled afresh), and correlation's own record
  -- the cluster, its data sets, the file it was cut from. Its report, and mxi
  ui's summaries, are its crystals'.
- mxi's Python tools are one wheel, `mxi-python`, where mxeq's was its own:
  mxeq, the referee, under the same import name and command; `mxi_tools`,
  which reads mxi's files and images; the viewer `mxi-view` (napari through
  the `[napari]` extra); and `mxi-ui`. matplotlib is now a dependency. The
  release's wheel is `mxi_python-1.1.4-py3-none-any.whl`; where the `mxeq`
  wheel was installed, `pip uninstall mxeq` first.

**New**
- `mxi process MASTER.nxs...`, for routine data: the whole sequence, import to
  export, a step at a time, each the command a script of them would run, with
  its own output and log. `-j` goes to every step, `--gpu` to spot finding and
  integration, `--cell` and `--space-group` to indexing (the cell as
  `--known-cell`) and the space group to scaling, `--no-joint` to indexing for
  several crystals, `--d-min-auto` to scaling. Given one `.rflx` instead, `mxi
  process FILE.rflx` picks up after the step that wrote it and finishes: a
  cluster from `mxi correlation --write-clusters` is scaled afresh and
  exported. A step that fails stops it, saying which, its exit code, its log,
  and how to go on -- its own command, or `mxi process` on its input;
  `--dry-run` prints the commands and runs none.
- `mxi ui`, an interface to mxi, for processing data and for teaching
  (docs/ui.md). The workflow's steps and tools down the side, each marked done
  from the working directory, so that a directory processed before opens as it
  stands. Each step a form -- its inputs, main options, output, threads and
  further options, each field's option a tooltip -- building the mxi command
  as it is filled in, shown to copy, and run with its output streamed into the
  Log, cancellable. Each step's Summary and Plots read from its file's
  records, not its log: its figures and tables, the merging statistics as
  DIALS lays them out; spots per image and by resolution, RMSDs by cycle, the
  cell over the scan, <I/sigma>, CC half, completeness and the R factors by
  resolution, with matplotlib. And View images: a step's images in a window of
  their own, each spot drawn as its box, indexed apart from unindexed,
  predictions; zoom by buttons, keys, a right-button drag or a scroll;
  contrast and colour maps; the pixel, its counts and its resolution under the
  cursor.
- `mxi ui` runs with tkinter, matplotlib, numpy, h5py and hdf5plugin -- no Qt
  -- and where tkinter is missing it says what to install for the system (`mxi
  ui --check`). It is in the release archives beside mxi, and in the wheel as
  `mxi-ui`.

**Not yet**
- Untrusted regions are still drawn in `mxi view` (napari) or set with `mxi
  untrusted`; the interface's viewer will draw them, and then replace `mxi
  view`.

