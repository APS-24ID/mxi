# mxi ui: the interface

An interface to mxi, for processing data and for teaching: the workflow's
steps down the side, each running the real mxi command and showing what it
wrote; the images, with their spots, predictions and untrusted regions; the
results as plots. It follows the DIALS workflow interface written before it --
a step a form that builds the command, the command run as it would be by
hand, its output streamed -- and grows to replace `mxi view`.

```sh
mxi ui                 # in the current directory
mxi ui ~/data/run1     # in another
mxi ui --check         # whether it can run here, and what is missing
```

From the wheel, `mxi-python`, the same is `mxi-ui`.

## What it needs

Python 3.10 or later with tkinter, Tk 8.6 or later; numpy, h5py, hdf5plugin
and matplotlib, which the wheel brings; a display. tkinter comes with Python,
but some systems package it apart, and `mxi ui` says what to install where it
is missing: `sudo apt install python3-tk` on Debian and Ubuntu, `sudo dnf
install python3-tkinter` on Fedora, RHEL and its clones, `brew install
python-tk@3.12` (the version of Python's) with Homebrew; the python.org
installers include it. There is no Qt anywhere, which is what made `mxi view`
heavy to install.

## How it works

The interface runs mxi, it does not reimplement any of it. Each step builds
the command one would type -- shown, so that it can be copied and learnt -- and
runs it in the working directory, in the background, its output streamed into
the Log tab, cancellable. It finds mxi beside itself, as a release archive
lays it out, else on the PATH, or where the `MXI` environment variable says.
What a step produced is read from its files directly, the images with h5py and
the `.rflx` files' tables and summaries likewise; nothing is scraped from a
log. A directory processed before, here or by hand, opens as it stands: a step
is done where its output is.

The package, `python/src/mxi_ui`, keeps tkinter to the window's modules, so
that the steps, the runner, the directory's state and the dependency check
import without it, and are tested where there is no display.

## A step's form

The window's right side is the chosen step's tabs: Setup, its form, first;
then its Log, Summary, Images and Plots. The command the form makes, with
Copy and Run, is in a bar beneath the tabs, in view whichever is shown, so
that a step can be run while its log or summary is read; running shows the
Log. Choosing another step keeps the tab, so that going through the steps'
summaries stays on the summaries.


Each step's form has its inputs (the previous step's output by default, a
Browse button to choose others; for Import, the images' master files), its
main options -- a tick box for a flag, a field for a value, a list where the
choices are few -- its output where the command takes `-o`, its threads (`-j`;
empty or 0, every core the job may use), and a line for any further options,
as `mxi COMMAND --help` lists them. Each field's option and what it does are
a tooltip, shown when the pointer rests on it -- `--d-min: nothing beyond
this d` -- not lines beside the fields. The command they make is shown beneath as
it is built, to copy or to learn from, and Run -- named for the step, "Run
scale", or Cmd-Return on a Mac and Ctrl-Return elsewhere -- runs exactly that.
A value that
cannot be right -- a resolution that is not a number -- is said in red, and
Run waits for it. Each step keeps its entries while the window is open.

The forms show each command's main options, not all of them; every option on
a form is checked against the command's own `--help` by the tests, so that
they cannot offer what mxi does not take.

## A step's results

What a step produced is read from its output file -- `strong.rflx` for spot
finding, `scaled.rflx` for scaling -- whose summary records every step so
far, each program's figures and tables (`docs/rflx.md`, `docs/output.md`);
the step's own record is its program's. Nothing is read from the log.

**Summary**: the step, its command, its file and when it ran, its figures --
those it knows by name first (spots, indexed fraction, cell, RMSDs, Laue
group, observations, the resolution at CC½ 0.3, the error model, ...), then any
other the record holds, so that a figure a command adds later is shown
without a change here -- and its tables: scaling's merging statistics as DIALS
lays them out, a row a figure and a column overall, low and high; the
subgroups symmetry weighed; integration's overall figures.

**Plots**, with matplotlib, its toolbar beneath to zoom, pan and save -- those
mxi report draws from the same records:

| Step | Plots |
| --- | --- |
| Find spots | spots per image, a line a sweep; spots by resolution |
| Index | the RMSDs, cycle by cycle |
| Refine | the RMSDs, cycle by cycle; each crystal's cell lengths and angles over the scan |
| Integrate | <I/σ> summed and profile-fitted, profile correlation, background, by resolution |
| Scale | CC½, <I/σ>, completeness with multiplicity on its own axis, the R factors on a log axis, by resolution |

A resolution axis is spaced evenly in 1/d², ticked in d. Not yet: each
sweep's own integration figures, where there are several; the clusters of
`mxi correlation`.

## The images

**View images**, beneath the steps and tools, opens a window of its own on
the chosen step's file -- its output where it has run, else what it reads, so
that after Import there are the images alone -- as the DIALS interface opened
its viewer for a step. On each image the spots and predictions the file
holds, on exactly the images each spans: spot finding's strong spots;
indexing's spots, those indexed (green) apart from those not (red); the
predictions (blue) of refinement and integration. A marker centred on the
image shown is solid, one spanning it from a neighbour dashed. Each spot is
drawn as its box, the pixels it covers, as DIALS's viewer draws them -- never
too small to see zoomed out -- and a prediction without a box (refinement's)
as a small circle; each layer can be hidden.

Step through the images with the arrows or the arrow keys, or type one; the
zoom is kept from image to image. Zoom with + and −, or the keys + and -
about the pointer; or drag up and down with the right button (or Control and
the left), a level each 30 pixels, about where the drag began; or a wheel or
a trackpad's scroll (Tk 9's own trackpad event too), a level a tenth of a
second at most. Drag with the left button to pan; Fit, 0 or a double-click
for the whole image. The keys are the image's, not the image box's: typing
an image number does not zoom. Contrast is set from the image, its darkest counts on the slider, with
three colour maps. Under the cursor: the pixel, its counts and its
resolution.

The frame is reduced to the window by the maximum over each block of pixels,
so that a spot one pixel wide shows at any zoom; the scale is a whole number
of pixels either way, so that the picture and its markers agree pixel for
pixel. On a 16M frame a view takes some 40 ms at the fit and under 10 zoomed
in; the contrast and colours re-colour what was reduced.

## Phases

1. **Foundations** (done): the package; the dependency check, its advice for
   each system; the main window -- the steps and tools down the side, done or
   not from the directory, tabs, a status line; running a step, streamed,
   cancellable; `mxi ui`, `mxi-ui` and `python -m mxi_ui`.
2. **The forms** (done): each step's inputs, main options, output, threads and
   further options, the command built as they are filled in.
3. **Results** (done, but for each sweep's integration and the clusters): each
   step's summary and plots from its file (above).
4. **The viewer** (begun): a window of its own for a step, its images, spots
   and predictions (above). Still to come: untrusted regions drawn, undone
   and saved through `mxi untrusted`; then `mxi view` opens it, and the napari
   viewer is retired.
5. **Teaching**: a short explanation of each step, and the guides linked.
6. **Packaging and documents**, released as 1.2.0.
