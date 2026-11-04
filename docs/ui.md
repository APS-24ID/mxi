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

Each step's form has its inputs (the previous step's output by default, a
Browse button to choose others; for Import, the images' master files), its
main options -- a tick box for a flag, a field for a value, a list where the
choices are few -- its output where the command takes `-o`, its threads (`-j`;
empty or 0, every core the job may use), and a line for any further options,
as `mxi COMMAND --help` lists them. The command they make is shown beneath as
it is built, to copy or to learn from, and Run runs exactly that. A value that
cannot be right -- a resolution that is not a number -- is said in red, and
Run waits for it. Each step keeps its entries while the window is open.

The forms show each command's main options, not all of them; every option on
a form is checked against the command's own `--help` by the tests, so that
they cannot offer what mxi does not take.

## Phases

1. **Foundations** (done): the package; the dependency check, its advice for
   each system; the main window -- the steps and tools down the side, done or
   not from the directory, tabs for images, log, summary and plots, a status
   line; running a step, streamed, cancellable; `mxi ui`, `mxi-ui` and
   `python -m mxi_ui`.
2. **The forms** (done): each step's inputs, main options, output, threads and
   further options, the command built as they are filled in (above).
3. **The viewer**: on a Tk canvas, each frame reduced to the window by the
   maximum over each block of pixels, so that a spot one pixel wide shows at
   any zoom; sweeps and frames, zoom and pan, contrast and colour maps, the
   pixel, value and resolution under the cursor; spots and predictions on
   exactly the images they occupy; untrusted regions drawn, undone and saved
   through `mxi untrusted`. Then `mxi view` opens it, and the napari viewer
   is retired.
4. **Results**: each step's summary from its file; several sweeps, several
   lattices, and clusters scaled one by one; plots with matplotlib as each
   step's results allow -- spots per image, refinement's RMSDs, I/σ and CC½
   against resolution, scale factors, the clusters.
5. **Teaching**: a short explanation of each step, and the guides linked.
6. **Packaging and documents**, released as 1.2.0.
