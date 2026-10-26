# Masking the detector: mxi_mask and mxi_untrusted

A backstop, its arm, a cryo nozzle, the goniometer: whatever shadows the
detector makes observations there wrong -- low where the spot is shadowed, high
where its background was -- and nothing in the detector's own pixel mask knows
of it. On a small molecule refined to R1 0.026, most of the observations beyond
5 sigma against Fc^2 lay in the backstop's shadow and along its arm
(`docs/outstanding.md`, item 78; `mxeq truth`). These two programs mark such
regions **untrusted**, so that spot finding and integration leave them out.

## Drawing them

```sh
mxi_import ../data/*.nxs
mxi_mask imported.rflx        # Save writes masked.rflx
mxi_find masked.rflx
```

`mxi_mask` shows each sweep in turn as a sum of images -- the first image and
how many set at the bottom, 10 to begin with, and the contrast -- so that a
shadow shows against the background. Draw:

* **rectangles**, dragged from a corner;
* **circles**, dragged from the centre;
* **polygons**, a click at each vertex, closed with a right click or Enter.

**Undo** takes back the last change (u or Ctrl-Z too); **Clear** removes this
sweep's regions; **Same as last** copies the previous sweep's -- or, on the
first, the mask last saved -- for sweeps that share their shadows; **Previous**
and **Next** move between sweeps; **Save** writes the experiment list with every
sweep's regions to a new file, `masked.rflx` (`-o` to choose), and keeps them in
`mxi_mask.json`. The input is left as it was, to compare with: masked and not,
processed side by side. To carry on editing, open `masked.rflx` and name it with
`-o` to write over it. While the toolbar's zoom or
pan is on, clicks are the toolbar's: zoom in to draw a fine edge.

It needs numpy, h5py, hdf5plugin and matplotlib. `--snapshot PNG` draws the
first sweep and its regions to a file, with no window.

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
