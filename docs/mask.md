# Viewing the images, and masking the detector: mxi_view and mxi_untrusted

A backstop, its arm, a cryo nozzle, the goniometer: whatever shadows the
detector makes observations there wrong -- low where the spot is shadowed, high
where its background was -- and nothing in the detector's own pixel mask knows
of it. On a small molecule refined to R1 0.026, most of the observations beyond
5 sigma against Fc^2 lay in the backstop's shadow and along its arm
(`docs/outstanding.md`, item 78; `mxeq truth`). These two programs mark such
regions **untrusted**, so that spot finding and integration leave them out.

## Viewing and drawing them: mxi_view

```sh
mxi_import ../data/*.nxs
mxi_view imported.rflx        # Save writes masked.rflx
mxi_find masked.rflx
mxi_view strong.rflx          # the strong spots on the images
```

`mxi_view` shows an experiment list's images in napari, a sweep at a time,
chosen in the panel on the right. The images are a stack with a slider; each
is shown as the sum of the number set in the panel from it (one to begin
with), so that a shadow shows against the background; contrast and colour
are napari's own. Whatever the file holds is drawn on them: the strong spots
(`xyzobs.px.value`) as red circles, and the predictions (`xyzcal.px`), hidden
to begin with, as cyan crosses -- each on its own images.

The untrusted regions are a shapes layer, `untrusted`, over every image of the
sweep: draw **rectangles**, **ellipses** -- a circle where round, a polygon of 48
sides where not -- and **polygons** with napari's tools, and edit, select,
delete and undo them as napari does. **Same as last** copies the regions of
the sweep shown before, for sweeps that share their shadows. **Save** writes
the experiment list with every sweep's regions to a new file, `masked.rflx`
(`-o` to choose), and keeps them in `mxi_mask.json`. The input is left as it
was, to compare with: masked and not, processed side by side.

napari is an optional dependency: `pip install "napari[all]" h5py
hdf5plugin`. Nothing else in mxi needs it, and `mxi_view` says what to install
where it is missing. It is built and installed with the other programs
(`cmake --install`, the README), and finds `mxi_convert`, `mxi_untrusted` and
its module, `mxi_tools.py`, beside itself. Over a remote desktop it needs one
that carries OpenGL, as NoMachine does.

## Where they are kept, and applying them again

Each image set of the experiment list keeps its own:

```json
"untrusted": [{"panel": 0, "circle": [729, 866, 40]},
              {"panel": 0, "rectangle": [720, 0, 760, 866]},
              {"panel": 0, "polygon": [[0, 0], [9, 0], [0, 9]]}]
```

in pixels as DIALS counts them, pixel i spanning i to i + 1: a pixel is
untrusted when its centre is inside a region. Every program carries the image
sets forward whole, so the regions go with the data from the imported file to
the scaled one; `mxi_find` and `mxi_integrate` each say how many pixels they
left out ("Untrusted regions: 1 circle, 1 rectangle, 12345 pixels").

`mxi_untrusted` lists, sets and clears them with no window -- the mask kept
from one run applied to the next:

```sh
mxi_untrusted imported.rflx                       # list each image set's
mxi_untrusted imported.rflx --set mxi_mask.json   # set them: masked.rflx
mxi_untrusted masked.rflx --clear -o clear.rflx   # remove them
```

Setting or clearing writes a new file, `masked.rflx` unless `-o` names another;
the input is changed only when named. A mask file is `{"imagesets": {"0": [...],
"1": [...]}}` for the image sets it
names, or a list of shapes alone for every image set.

## How they are applied

Rasterised once a sweep (`src/untrusted.cc`), the regions join the detector's
own pixel_mask in the image series (`Series::add_masked`), so every frame read
afterwards carries them: each such pixel becomes the bad-pixel marker, which
spot finding and integration already take for no measurement, on the CPU and
on a device alike. `mxi_max` leaves them in, to show what was masked. Panel 0
only, as frames are read here; a region on another panel is refused.
