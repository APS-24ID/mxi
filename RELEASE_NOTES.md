# Release notes

Newest first. Each release lists first what changes results or commands -- a
default changed, an option renamed or gone -- then what is new, then what is
fixed. The tag of each (`git tag -n99 v1.0.5`) carries the same summary.

## Unreleased -- since 1.1.3

**New**
- `mxi process MASTER.nxs...`, for routine data: the whole sequence, import to
  export, a step at a time, each the command a script of them would run, with
  its own output and log. `-j` goes to every step, `--gpu` to spot finding and
  integration, `--cell` and `--space-group` to indexing (the cell as
  `--known-cell`) and the space group to scaling, `--no-joint` to indexing for
  several crystals, `--d-min-auto` to scaling. A step that fails stops it,
  saying which, its exit code, its log, and the command to run it again;
  `--dry-run` prints the eight commands and runs none.
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

**Changed**
- mxi's Python tools are one wheel, `mxi-python`, where mxeq's was its own:
  mxeq, the referee, under the same import name and command; `mxi_tools`,
  which reads mxi's files and images; the viewer `mxi-view` (napari through
  the `[napari]` extra); and `mxi-ui`. matplotlib is now a dependency.

**Not yet**
- Untrusted regions are still drawn in `mxi view` (napari) or set with `mxi
  untrusted`; the interface's viewer will draw them, and then replace `mxi
  view`.

