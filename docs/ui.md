# mxi ui: the interface

An interface to mxi, for processing data and for teaching. The workflow's
steps are down the side; each is a form that builds the real `mxi` command,
shows it, runs it and shows what it wrote -- its log, its summary, its plots
-- and each step's images open in a viewer of their own, with its spots and
predictions on them. It runs mxi rather than reimplementing it, as the DIALS
workflow interface before it ran DIALS: what it does by a click can be done,
and learnt, by hand.

## Starting it

```sh
mxi ui                 # in the current directory
mxi ui ~/data/run1     # in another
mxi ui --check         # whether it can run here, and what is missing
```

From a release archive `mxi ui` runs the interface beside `mxi`; from the
wheel, `mxi-python`, the same is `mxi-ui`, or `python -m mxi_ui`.

## What it needs

Python 3.10 or later with tkinter, Tk 8.6 or later; numpy, h5py, hdf5plugin
and matplotlib, which the wheel brings; a display. tkinter comes with Python,
but some systems package it apart, and `mxi ui` says what to install where it
is missing: `sudo apt install python3-tk` on Debian and Ubuntu, `sudo dnf
install python3-tkinter` on Fedora, RHEL and its clones, `brew install
python-tk@3.12` (the version of Python's) with Homebrew; the python.org
installers include it. There is no Qt.

## The window

Down the left, the steps of processing -- Import, Find spots, Index, Refine,
Integrate, Symmetry, Scale, Export -- and the tools beside them -- Verify a
lattice, Correlation, Report -- each marked done where its output is in the
working directory. A directory processed before, here or by hand, opens as
it stands, on the first step not yet done. Beneath them, **View images**
(below).

On the right, the chosen step's tabs: **Setup**, its form; **Log**, what its
command printed; **Summary** and **Plots**, what its file records. Choosing
another step keeps the tab, so that going through the steps' summaries stays
on the summaries. Beneath the tabs, in view whichever is shown, the command
the form makes, with **Copy** and **Run**; at the foot, the status and
**Cancel**.

## Running a step

A step's form has its inputs -- the previous step's output by default, Browse
to choose others; for Import, the images' master files -- its main options,
a tick box for a flag, a field for a value, a list where the choices are few;
its output where the command takes `-o`; its threads (`-j`; empty or 0, every
core the job may use); and a line for any further options, as `mxi COMMAND
--help` lists them. Resting the pointer on a field shows its option and what
it does -- `--d-min: nothing beyond this d`.

The command they make is shown beneath as it is built, to copy or to learn
from, and Run -- named for the step, "Run scale"; Cmd-Return on a Mac,
Ctrl-Return elsewhere -- runs exactly that, in the working directory, in the
background, its output streamed into the Log, cancellable. A value that
cannot be right -- a resolution that is not a number -- is said in red, and
Run waits for it. Each step keeps its entries while the window is open.

The forms show each command's main options, not all of them; the rest go in
the line for further options. Every option on a form is checked against the
command's own `--help` by the tests, so that a form cannot offer what mxi
does not take.

## A step's results

What a step produced is read from its output file -- `strong.rflx` for spot
finding, `scaled.rflx` for scaling -- whose summary records every step so far,
each command's figures and tables (`docs/rflx.md`, `docs/output.md`). Nothing
is read from the log.

**Summary**: the step, its command, its file and when it ran, its figures --
those it knows by name first (spots, indexed fraction, cell, RMSDs, Laue
group, observations, the resolution at CC½ 0.3, the error model, ...), then any
other the record holds, so that a figure a command adds later is shown -- and
its tables: scaling's merging statistics as DIALS lays them out, a row a
figure and a column overall, low and high; the candidate cells indexing
weighed; the subgroups symmetry weighed; integration's overall figures. A file
that cannot be read says why; a file without the step's record says which
steps it does record.

**Plots**, with matplotlib, its toolbar beneath to zoom, pan and save -- those
`mxi report` draws from the same records:

| Step | Plots |
| --- | --- |
| Find spots | spots per image, a line a sweep; spots by resolution |
| Index | the RMSDs, cycle by cycle |
| Refine | the RMSDs, cycle by cycle; each crystal's cell lengths and angles over the scan |
| Integrate | <I/σ> summed and profile-fitted, profile correlation, background, by resolution; with several sweeps, <I/σ> and profile correlation a line a sweep |
| Scale | CC½, <I/σ>, completeness with multiplicity on its own axis, the R factors, by resolution |
| Correlation | the correlation matrix in OPTICS's order; the reachability, coloured by cluster; the data sets on their first two principal components; the functional by dimension -- from its `mxi_correlation.json`, the tool always asking for it |

A resolution axis is spaced evenly in 1/d², ticked in d.

## The images

**View images** opens a window of its own on the chosen step's file -- its
output where it has run, else what it reads, so that after Import there are
the images alone -- as the DIALS interface opened its viewer for a step;
several can be open at once. On each image, the spots and predictions the
file holds, on exactly the images each spans: spot finding's strong spots;
indexing's spots, those indexed (green) apart from those not (red); the
predictions (blue) of refinement and integration. Each spot is drawn as its
box, the pixels it covers, as DIALS's viewer draws them -- never too small to
see zoomed out -- solid where centred on the image shown, dashed where it
spans it from a neighbour; a prediction without a box (refinement's) is a
small circle; each layer can be hidden.

- **Images**: ◀ and ▶, the arrow keys, or a number typed; where there are
  several sweeps, a list of them. The zoom is kept from image to image.
- **Sum**: the images summed from the one shown -- 10 or 50 shows a shadow
  against the summed background, and faint spots, where one image shows
  neither; the spots on any of them marked, the contrast set afresh, the sum
  stopping at the sweep's last image.
- **Zoom**: + and −, or the keys + and - (about the pointer); or drag up and
  down with the right button, or with Control and the left, a level each 30
  pixels, about where the drag began; or a wheel or a trackpad's scroll.
  Fit, 0 or a double-click shows the whole image, as the window opens.
- **Pan**: drag with the left button.
- **Contrast**: set from the image, its darkest counts on the slider; three
  colour maps.
- **Under the cursor**: the pixel, its counts and its resolution.

### Untrusted regions

Beneath the layers, the detector's untrusted regions -- a backstop's shadow,
its arm, a cryo nozzle -- which spot finding and integration then leave out
(`docs/mask.md`). The bar leads with how many the sweep has, and whether they
are saved. The left button draws what is chosen:

- **Look**: drag to pan; click a region to choose it (yellow), Delete or
  Backspace removes it.
- **Rectangle**: drag corner to corner. **Circle**: drag from its centre to
  its edge. **Polygon**: click each vertex, a double-click or Return closing
  it, Escape abandoning it. Shift-drag pans while drawing.

**Undo** (Ctrl-Z, Cmd-Z on a Mac) puts the sweep's regions back as they were
before the last change, a step at a time, each sweep its own; **Clear**
removes the sweep's; **Same as previous** copies those of the sweep shown
before, for sweeps sharing their shadows. **Save** writes `masked.rflx` beside
the file shown, with every sweep's regions, through `mxi untrusted`, and keeps
them in `mxi_mask.json`; the file shown is left as it was. Find spots then
reads `masked.rflx`. Closing the viewer with regions not saved asks first. The
regions a file holds already are shown, to change or add to.

Each pixel shown is the maximum of the detector pixels it stands for, so that
a spot one pixel wide shows at any zoom; the scale is a whole number of
pixels either way, so that the picture and its boxes agree pixel for pixel.
On a frame the size of an Eiger 16M's a view takes some 40 ms at the fit and
under 10 zoomed in (measured on a synthetic frame); a change of contrast or
colours re-colours what was reduced.

## How it works

The package, `python/src/mxi_ui`, in the wheel and beside `mxi` in the
release archives: `steps.py` and `model.py` the steps and their forms'
options, and the command a form builds; `state.py` the steps a directory has
done; `runner.py` a command in the background; `results.py` a step's facts,
tables and plots from its file's summary; `render.py` the viewer's
reduction, scale and geometry; `deps.py` the check for what is needed;
`window.py` and `viewer.py` the windows. Only the windows import tkinter, so
that everything else imports, and is tested, where there is no display. It
finds mxi beside itself, as a release archive lays it out, else on the PATH,
or where the `MXI` environment variable says. It reads images and files
through `mxi_tools`, as `mxi view` does.

## Not yet

- The resolution under the cursor from a multi-panel detector's first panel
  only.
- A short explanation of each step, and the guides linked, for teaching.
