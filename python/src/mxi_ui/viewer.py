"""The image viewer: a window of its own, opened for a step from the main
window's View images, showing that step's file -- its images, and on them the
spots and predictions it holds: strong spots after spot finding, indexed and
not after indexing, predictions after refinement and integration -- each on
exactly the images it spans (mxi_tools.Markers).

The frame is reduced to the window by the maximum over each block of pixels
(render.py), so that a spot one pixel wide shows at any zoom; zoom with the
wheel about the cursor, drag to pan, double-click to fit. Under the cursor:
the pixel, its counts and its resolution. Tk only, with numpy and h5py.
"""

import os
import sys
import time
import tkinter as tk
from tkinter import ttk

import numpy as np

import mxi_tools

from mxi_ui import gestures
from mxi_ui import regions as shapes
from mxi_ui import render
from mxi_ui.tooltip import Tooltip

MIN_RADIUS = 5  # canvas pixels: the least a marker is, however small its spot
REGION = "#ff6d00"  # an untrusted region's outline
CHOSEN = "#ffd600"  # the one chosen, to delete
MODES = (
    ("look", "Look"),
    ("rectangle", "Rectangle"),
    ("circle", "Circle"),
    ("polygon", "Polygon"),
)
DRAG_LEVEL = 30  # canvas pixels dragged a zoom level, with the right button
WHEEL_INTERVAL = 0.1  # seconds: at most a zoom level each, from a scroll
# Several lattices on a sweep: a colour each, its spots strong, its
# predictions a lighter tint of it -- matplotlib's tab20 pairs, red left out,
# red being the unindexed spots' (Graeme: boxes of different lattices told
# apart).
LATTICE_COLOURS = (
    ("#2ca02c", "#98df8a"),  # green
    ("#1f77b4", "#aec7e8"),  # blue
    ("#ff7f0e", "#ffbb78"),  # orange
    ("#9467bd", "#c5b0d5"),  # purple
    ("#8c564b", "#c49c94"),  # brown
    ("#e377c2", "#f7b6d2"),  # pink
    ("#bcbd22", "#dbdb8d"),  # olive
    ("#17becf", "#9edae5"),  # cyan
    ("#7f7f7f", "#c7c7c7"),  # grey
)
COLUMNS = ["id", "xyzobs.px.value", "xyzcal.px", "bbox", "miller_index"]


def _at_least(lo, hi, least):
    """lo to hi widened about its middle to `least` where it is narrower."""
    if hi - lo >= least:
        return lo, hi
    middle = (lo + hi) / 2
    return middle - least / 2, middle + least / 2


def typing(event):
    """Whether a key came while typing in a field -- the image box, a list --
    where it is the field's, not the viewer's."""
    return event.widget.winfo_class() in (
        "TSpinbox",
        "TEntry",
        "TCombobox",
        "Entry",
        "Spinbox",
    )


def masters(document, path):
    """Each sweep (mxi_tools.sweeps), its master made absolute: a relative one
    is the file's directory's."""
    here = os.path.dirname(os.path.abspath(path))
    out = mxi_tools.sweeps(document)
    for s in out:
        if s["master"] and not os.path.isabs(s["master"]):
            s["master"] = os.path.join(here, s["master"])
    return out


def layers_for(table, experiments, lattices=()):
    """The marker layers a sweep's rows make, and how each is drawn: ({key:
    mxi_tools.Markers or None}, [(key, label, colour)] in drawing order) --
    each key one word, as a canvas tag must be.
    With one lattice or none: spots indexed and not where the table has
    indices -- else all its spots as one -- and predictions where it has them.
    With several (`lattices`, the sweep's experiments with a crystal): each
    lattice's indexed spots and its predictions, a colour each, the spots
    unindexed apart."""
    if "id" not in table:
        return {}, []
    ids = np.asarray(table["id"])
    mine = np.isin(ids, list(experiments))
    hkl = (
        np.asarray(table["miller_index"]).reshape(len(ids), -1)
        if "miller_index" in table
        else None
    )
    indexed = mine & np.any(hkl != 0, axis=1) if hkl is not None else None
    out, styles = {}, []

    def layer(key, label, colour, rows, column):
        sub = {k: np.asarray(v)[rows] for k, v in table.items()}
        markers = mxi_tools.Markers(sub, column, experiments)
        out[key] = markers if len(markers) else None
        styles.append((key, label, colour))

    if len(lattices) > 1:
        for n, e in enumerate(lattices):
            spots, predictions = LATTICE_COLOURS[n % len(LATTICE_COLOURS)]
            if "xyzobs.px.value" in table and indexed is not None:
                layer(
                    f"lattice{n + 1}",
                    f"Lattice {n + 1}",
                    spots,
                    indexed & (ids == e),
                    "xyzobs.px.value",
                )
            if "xyzcal.px" in table:
                layer(
                    f"predicted{n + 1}",
                    f"Lattice {n + 1}, predicted",
                    predictions,
                    mine & (ids == e),
                    "xyzcal.px",
                )
        if "xyzobs.px.value" in table:
            layer(
                "unindexed",
                "Spots",
                "#d32f2f",
                mine & ~indexed if indexed is not None else mine,
                "xyzobs.px.value",
            )
        return out, styles
    if "xyzobs.px.value" in table:
        if indexed is not None:
            layer("indexed", "Indexed spots", "#2e7d32", indexed, "xyzobs.px.value")
            layer("unindexed", "Spots", "#d32f2f", mine & ~indexed, "xyzobs.px.value")
        else:
            layer("unindexed", "Spots", "#d32f2f", mine, "xyzobs.px.value")
    if "xyzcal.px" in table:
        layer("predicted", "Predictions", "#1565c0", mine, "xyzcal.px")
    return out, styles


class Viewer:
    def __init__(
        self,
        parent,
        path,
        title=None,
        on_saved=None,
        output="masked.rflx",
        where=None,
        summing=1,
        sweep=0,
    ):
        self.path = path
        self.on_saved = on_saved  # told masked.rflx's path when saved
        # Save's file, and the directory it and mxi_mask.json go in: beside
        # the file shown (mxi ui), or where mxi view is run.
        self.output = output
        self.where = where or os.path.dirname(os.path.abspath(path))
        self.document = mxi_tools.read_experiments(path)
        self.sweeps = masters(self.document, path)
        if not self.sweeps:
            raise ValueError(f"{os.path.basename(path)} names no images")
        self.table = (
            mxi_tools.read_reflections(path, COLUMNS) if path.endswith(".rflx") else {}
        )
        self.images = {}
        self.layers = {}
        # The untrusted regions each sweep's image set holds already; drawn,
        # undone and saved here (regions.py, docs/mask.md).
        self.regions = shapes.Regions(
            {s["index"]: s.get("shapes", []) for s in self.sweeps}
        )
        self.mode = None
        self.chosen = None
        self._shape = None  # a rectangle's or circle's first point, being drawn
        self._vertices = []  # a polygon's, being drawn
        self._previous_sweep = None
        self.saved = False
        self.top = tk.Toplevel(parent)
        self.top.title(title or f"mxi ui: {os.path.basename(path)}")
        self.top.geometry("1000x900")
        self.top.protocol("WM_DELETE_WINDOW", self.close)
        self.top.bind(
            "<Command-w>" if sys.platform == "darwin" else "<Control-w>",
            lambda _e: self.close(),
        )
        self.sweep = tk.IntVar(value=0)
        self.image = tk.IntVar(value=1)
        self.summing = tk.IntVar(value=max(1, int(summing)))
        self.ceiling = tk.DoubleVar(value=10.0)
        self.scheme = tk.StringVar(value="grey")
        self.showing = {}  # each layer's tick, made as layers are met
        self.styles = {}  # each sweep's layers: (key, label, colour)
        self.readout = tk.StringVar(value="")
        self.scan_at = tk.DoubleVar(value=1.0)  # the scan slider's place
        self._scan_wanted = None  # the image the slider last asked for
        self._scan_due = False
        self.view = None
        self.photo = None
        self.frame = None
        self._drag = None
        self._reduced = None
        self._last_wheel = -1e9
        self._zooming = None
        self._auto_fit = True
        self._controls()
        self.canvas = tk.Canvas(self.top, background="#808080", highlightthickness=0)
        self.canvas.pack(fill=tk.BOTH, expand=True)
        self._scan_bar()
        ttk.Label(
            self.top, textvariable=self.readout, padding=(8, 3), font="TkFixedFont"
        ).pack(fill=tk.X)
        self._bindings()
        self.choose_sweep(min(max(0, int(sweep)), len(self.sweeps) - 1))

    # --- building ---------------------------------------------------------

    def _scan_bar(self):
        """A slider the width of the window, under the image, through the
        sweep's images (Graeme: to scan through them). Dragged, it shows the
        image it is at as fast as images are read -- each drawn when the
        window is next idle, at the latest place, the events of a slow read
        come to one -- and it follows the arrows, the box and the keys."""
        row = ttk.Frame(self.top, padding=(8, 4, 8, 0))
        row.pack(fill=tk.X)
        self.scan = ttk.Scale(
            row,
            orient=tk.HORIZONTAL,
            from_=1,
            to=1,
            variable=self.scan_at,
            command=self._scanned,
        )
        self.scan.pack(fill=tk.X, expand=True)
        Tooltip(self.scan, "Drag through the sweep's images")

    def _scanned(self, value):
        try:
            wanted = int(round(float(value)))
        except (TypeError, ValueError):
            return
        self._scan_wanted = min(self.last, max(self.first, wanted))
        if not self._scan_due:
            self._scan_due = True
            self.top.after_idle(self._scan_show)

    def _scan_show(self):
        self._scan_due = False
        wanted = self._scan_wanted
        if wanted is None or wanted == self._image():
            return
        self.image.set(wanted)
        self.show()

    def _controls(self):
        bar = ttk.Frame(self.top, padding=(8, 6))
        bar.pack(fill=tk.X)
        if len(self.sweeps) > 1:
            ttk.Label(bar, text="Sweep").pack(side=tk.LEFT)
            names = [
                f"{s['index']}: {os.path.basename(s['master'])}" for s in self.sweeps
            ]
            box = ttk.Combobox(bar, values=names, state="readonly", width=24)
            box.current(0)
            box.bind(
                "<<ComboboxSelected>>", lambda _e: self.choose_sweep(box.current())
            )
            box.pack(side=tk.LEFT, padx=(4, 12))
        ttk.Button(bar, text="◀", width=3, command=lambda: self.step(-1)).pack(
            side=tk.LEFT
        )
        self.spin = ttk.Spinbox(
            bar, textvariable=self.image, width=7, command=self.show
        )
        self.spin.pack(side=tk.LEFT, padx=2)
        self.spin.bind("<Return>", lambda _e: self.show())
        ttk.Button(bar, text="▶", width=3, command=lambda: self.step(1)).pack(
            side=tk.LEFT
        )
        self.of = ttk.Label(bar, text="")
        self.of.pack(side=tk.LEFT, padx=(4, 8))
        # Images summed from the one shown: a shadow shows against the summed
        # background, and faint spots, where one image shows neither.
        ttk.Label(bar, text="Sum").pack(side=tk.LEFT)
        sum_box = ttk.Spinbox(
            bar,
            textvariable=self.summing,
            from_=1,
            to=1000,
            width=4,
            command=lambda: self.show(fit_contrast=True),
        )
        sum_box.pack(side=tk.LEFT, padx=(2, 12))
        sum_box.bind("<Return>", lambda _e: self.show(fit_contrast=True))
        Tooltip(
            sum_box,
            "Images summed from the one shown: 10 or 50 shows a shadow "
            "against the background; the spots on any of them marked",
        )
        ttk.Label(bar, text="Contrast").pack(side=tk.LEFT)
        self.slider = ttk.Scale(
            bar,
            from_=1,
            to=1000,
            variable=self.ceiling,
            length=140,
            command=lambda _v: self.draw(),
        )
        self.slider.pack(side=tk.LEFT, padx=4)
        Tooltip(self.slider, "The counts drawn darkest: lower to see faint spots")
        maps = ttk.Combobox(
            bar,
            textvariable=self.scheme,
            values=render.COLOUR_MAPS,
            state="readonly",
            width=9,
        )
        maps.bind("<<ComboboxSelected>>", lambda _e: self.draw())
        maps.pack(side=tk.LEFT, padx=(4, 12))
        # Buttons and keys to zoom, not a wheel alone: there need be no wheel
        # (Graeme, a trackpad).
        zoom_in = ttk.Button(bar, text="+", width=3, command=lambda: self.zoom(1))
        zoom_in.pack(side=tk.LEFT)
        zoom_out = ttk.Button(bar, text="−", width=3, command=lambda: self.zoom(-1))
        zoom_out.pack(side=tk.LEFT, padx=(2, 4))
        fit = ttk.Button(bar, text="Fit", width=4, command=self.fit)
        fit.pack(side=tk.LEFT)
        # Close, in sight -- the program's end under mxi view -- asking first
        # where regions are not saved (close()).
        done = ttk.Button(bar, text="Close", command=self.close)
        done.pack(side=tk.RIGHT)
        Tooltip(done, "Close the viewer: Cmd-W or Ctrl-W")
        Tooltip(
            zoom_in,
            "Zoom in: + or =; or drag up with the right button (or Control-drag)",
        )
        Tooltip(zoom_out, "Zoom out: -; or drag down with the right button")
        Tooltip(fit, "The whole image: 0, or double-click")
        self.toggles = ttk.Frame(self.top, padding=(8, 0, 8, 4))
        self.toggles.pack(fill=tk.X)
        self._regions_bar()

    def _regions_bar(self):
        """Untrusted regions: what the left button draws, and what to do
        with them (docs/mask.md)."""
        bar = ttk.Frame(self.top, padding=(8, 0, 8, 4))
        bar.pack(fill=tk.X)
        # Its count leading the bar, in view however narrow the window: at its
        # end, it had been cut off.
        self.count = ttk.Label(bar, text="Untrusted regions")
        self.count.pack(side=tk.LEFT, padx=(0, 6))
        self.mode = tk.StringVar(value="look")
        tips = {
            "look": "Drag to pan; click a region to choose it, Delete to remove it",
            "rectangle": "Drag corner to corner; Shift-drag pans",
            "circle": "Drag from its centre to its edge; Shift-drag pans",
            "polygon": "Click each vertex; double-click or Return closes it, Escape abandons it",
        }
        for key, label in MODES:
            b = ttk.Radiobutton(
                bar,
                text=label,
                value=key,
                variable=self.mode,
                command=self._mode_changed,
            )
            b.pack(side=tk.LEFT)
            Tooltip(b, tips[key])
        for label, action, tip in (
            (
                "Undo",
                self.undo,
                "The sweep's regions as before the last change: Ctrl-Z or Cmd-Z",
            ),
            ("Clear", self.clear, "Remove the sweep's regions"),
            (
                "Same as previous",
                self.same_as_previous,
                "The regions of the sweep shown before, for sweeps sharing their shadows",
            ),
            (
                "Save",
                self.save,
                "masked.rflx, beside the file, with every sweep's regions; "
                "and mxi_mask.json: the file shown is left as it was",
            ),
        ):
            b = ttk.Button(bar, text=label, command=action)
            b.pack(side=tk.LEFT, padx=(8 if label == "Undo" else 2, 0))
            Tooltip(b, tip)

    def _bindings(self):
        c = self.canvas
        c.bind("<Configure>", lambda _e: self._resized())
        c.bind("<ButtonPress-1>", self._press)
        c.bind("<B1-Motion>", self._move)
        c.bind("<ButtonRelease-1>", self._release)
        c.bind("<Shift-ButtonPress-1>", lambda e: self._press(e, pan=True))
        c.bind("<Shift-B1-Motion>", self._move)
        # Zoom by dragging up or down with the right button -- Button-2 on a
        # Mac, Button-3 elsewhere -- or Control and the left: no wheel needed,
        # and a trackpad's two fingers had done nothing (Graeme).
        for press, motion in (
            ("<ButtonPress-2>", "<B2-Motion>"),
            ("<ButtonPress-3>", "<B3-Motion>"),
            ("<Control-ButtonPress-1>", "<Control-B1-Motion>"),
        ):
            c.bind(press, self._zoom_press)
            c.bind(motion, self._zoom_drag)
        c.bind("<Double-Button-1>", self._double)
        c.bind("<Motion>", self._hover)
        c.bind("<MouseWheel>", lambda e: self._wheel(e, 1 if e.delta > 0 else -1))
        c.bind("<Button-4>", lambda e: self._wheel(e, 1))  # X11's wheel
        c.bind("<Button-5>", lambda e: self._wheel(e, -1))
        # Tk 9's trackpad scroll (Homebrew's Python has it), its own event;
        # Tk 8.6 knows no such event.
        try:
            c.bind("<TouchpadScroll>", self._touchpad)
        except tk.TclError:
            pass
        # The keys, but not while typing in a field: bound to the window, 0
        # typed into the image box -- image 10, 100 -- had fitted the view.
        keys = {
            "<Left>": lambda: self.step(-1),
            "<Right>": lambda: self.step(1),
            "<KeyPress-0>": self.fit,
        }
        for key in ("plus", "equal", "KP_Add"):
            keys[f"<KeyPress-{key}>"] = lambda: self.zoom(1, at_pointer=True)
        for key in ("minus", "KP_Subtract"):
            keys[f"<KeyPress-{key}>"] = lambda: self.zoom(-1, at_pointer=True)
        keys["<Control-z>"] = self.undo
        keys["<Command-z>" if sys.platform == "darwin" else "<Control-Z>"] = self.undo
        keys["<Delete>"] = self.delete_chosen
        keys["<BackSpace>"] = self.delete_chosen
        keys["<Return>"] = self._close_polygon
        keys["<Escape>"] = self._abandon
        for key, action in keys.items():
            self.top.bind(key, lambda e, act=action: None if typing(e) else act())

    # --- the sweep and image ----------------------------------------------

    def choose_sweep(self, k):
        if self.view is not None and k != self.sweep.get():
            self._previous_sweep = self.sweep.get()
        self.chosen = None
        self._vertices = []
        self.sweep.set(k)
        s = self.sweeps[k]
        if s["master"] not in self.images:
            self.images[s["master"]] = mxi_tools.Images(s["master"])
        images = self.images[s["master"]]
        last = s["last"] if s["last"] is not None else len(images)
        self.first, self.last = s["first"], last
        self.spin.configure(from_=self.first, to=self.last)
        self.of.configure(text=f"of {self.first} to {self.last}")
        self.scan.configure(from_=self.first, to=self.last)
        if k not in self.layers:
            lattices = [e for e in s["experiments"] if self._has_crystal(e)]
            self.layers[k], self.styles[k] = layers_for(
                self.table, s["experiments"], lattices
            )
        for child in self.toggles.winfo_children():
            child.destroy()
        for key, label, colour in self.styles[k]:
            self.showing.setdefault(key, tk.BooleanVar(value=True))
            markers = self.layers[k].get(key)
            if markers is None:
                continue
            box = tk.Checkbutton(
                self.toggles,
                text=f"{label} ({len(markers)})",
                variable=self.showing[key],
                command=self.draw,
                foreground=colour,
                activeforeground=colour,
            )
            box.pack(side=tk.LEFT, padx=(0, 12))
        h, w = images.frame_shape
        self.view = render.View((w, h), self._canvas_size())
        self._auto_fit = True
        self.image.set(self.first)
        self.show(fit_contrast=True)

    def step(self, by):
        self.image.set(min(self.last, max(self.first, self._image() + by)))
        self.show()

    def _which(self):
        """The image shown, or the images summed: "images 3 to 12"."""
        n, k = self._image(), self._summed()
        return f"image {n}" if k == 1 else f"images {n} to {n + k - 1}"

    def _summed(self):
        """How many images the one shown is summed with, from it: at least
        one, and none past the sweep's last."""
        try:
            k = int(self.summing.get())
        except (tk.TclError, ValueError):
            k = 1
        k = max(
            1, min(k, self.last - min(self.last, max(self.first, self._image())) + 1)
        )
        if k != self._summing_shown():
            self.summing.set(k)
        return k

    def _summing_shown(self):
        try:
            return int(self.summing.get())
        except (tk.TclError, ValueError):
            return 0

    def _has_crystal(self, experiment):
        try:
            return "crystal" in self.document["experiment"][experiment]
        except (KeyError, IndexError, TypeError):
            return False

    def _image(self):
        try:
            return int(self.image.get())
        except (tk.TclError, ValueError):
            return self.first

    def show(self, fit_contrast=False):
        """The image numbered in the spin box, read: frame image - 1 of the
        sweep's master, as mxi numbers them."""
        n = min(self.last, max(self.first, self._image()))
        self.image.set(n)
        self.scan_at.set(n)  # the slider follows, its command not called
        s = self.sweeps[self.sweep.get()]
        images = self.images[s["master"]]
        # The sum ends at the sweep's last image (mxi_tools.Images.summed;
        # its bad and overloaded pixels zero, adding nothing).
        images.window = self._summed()
        self.frame = images.summed(n - 1) if images.window > 1 else images.frame(n - 1)
        if fit_contrast:
            top = mxi_tools.display_top(self.frame)
            self.slider.configure(to=max(10.0, 20 * top))
            self.ceiling.set(top)
        self.draw()

    # --- drawing ------------------------------------------------------------

    def _canvas_size(self):
        return max(1, self.canvas.winfo_width()), max(1, self.canvas.winfo_height())

    def fit(self):
        if self.view is not None:
            self.view.fit(self._canvas_size())
            self._auto_fit = True
            self.draw()

    def _resized(self):
        """The canvas's size changed -- first of all when the window is first
        drawn, before which it measures a pixel: refitted, until the view is
        zoomed or panned, after which it is kept. Fitted once at a pixel's
        size and never again, as at first, the window had opened grey, Fit
        needed to see the image (Graeme)."""
        if self.view is not None and self._auto_fit:
            self.view.fit(self._canvas_size())
        self.draw()

    def draw(self):
        if self.frame is None or self.view is None:
            return
        size = self._canvas_size()
        if size[0] < 2 or size[1] < 2:
            return
        # The reduction kept while only the contrast or colours change: the
        # slider then re-colours alone, a few milliseconds.
        key = (id(self.frame), self.view.region(size), self.view.scale)
        if self._reduced is None or self._reduced[0] != key:
            self._reduced = (key, render.reduce(self.frame, key[1], key[2]))
        reduced, (x0, y0) = self._reduced[1]
        rgb = render.colour(reduced, self.ceiling.get(), self.scheme.get())
        # Pixel for pixel on the canvas (the scale is whole: render.snap),
        # its top left where its first detector pixel is.
        self.photo = tk.PhotoImage(data=render.ppm(rgb), format="PPM")
        cx, cy = self.view.to_canvas(x0, y0)
        self.canvas.delete("all")
        self.canvas.create_image(
            round(cx), round(cy), image=self.photo, anchor=tk.NW, tags="image"
        )
        self._markers()
        self._draw_regions()

    def _markers(self):
        """Each spot as its box -- the pixels it covers, its edges a pixel
        outside so that the outline frames the spot -- as DIALS's viewer draws
        them; never less than MIN_RADIUS about the box's middle, so that a
        spot shows zoomed out. A circle where there is no box (refinement's
        predictions). Circles their spot's size, as just before, had stood for
        boxes the file has (Graeme)."""
        layers = self.layers.get(self.sweep.get(), {})
        frame = self._image() - 1
        width, height = self._canvas_size()
        for key, _label, colour in self.styles.get(self.sweep.get(), []):
            markers = layers.get(key)
            if markers is None or not self.showing[key].get():
                continue
            yx, centred, boxes = markers.on_boxed(frame, self._summed())
            for i, ((y, x), own) in enumerate(zip(yx, centred)):
                style = dict(
                    outline=colour,
                    width=2 if own else 1,
                    dash=() if own else (2, 2),
                    tags=key,
                )
                if boxes is not None:
                    bx0, bx1, by0, by1 = boxes[i]
                    x0, y0 = self.view.to_canvas(bx0, by0)
                    x1, y1 = self.view.to_canvas(bx1, by1)
                    x0, x1 = _at_least(x0 - 1, x1 + 1, 2 * MIN_RADIUS)
                    y0, y1 = _at_least(y0 - 1, y1 + 1, 2 * MIN_RADIUS)
                    if x1 >= 0 and x0 <= width and y1 >= 0 and y0 <= height:
                        self.canvas.create_rectangle(x0, y0, x1, y1, **style)
                else:
                    # Markers' positions are pixel centres', half a pixel less than mxi's
                    cx, cy = self.view.to_canvas(x + 0.5, y + 0.5)
                    r = MIN_RADIUS
                    if -r <= cx <= width + r and -r <= cy <= height + r:
                        self.canvas.create_oval(cx - r, cy - r, cx + r, cy + r, **style)

    # --- the mouse ----------------------------------------------------------

    def _press(self, e, pan=False):
        self.canvas.focus_set()  # the keys the image's, not the image box's
        self._moved = False
        mode = self.mode.get() if self.mode is not None else "look"
        if pan or mode == "look":
            self._drag = (e.x, e.y)
            self._shape = None
            return
        self._drag = None
        x, y = self.view.to_detector(e.x, e.y)
        if mode == "polygon":
            self._vertices.append((x, y))
            self.draw()
        else:
            self._shape = (x, y)

    def _zoom_press(self, e):
        self._zooming = (e.x, e.y, e.y)  # where, and the last y a level taken
        self.canvas.focus_set()
        return "break"

    def _zoom_drag(self, e):
        """Up a level in, down a level out, each DRAG_LEVEL pixels, about
        where the drag began."""
        if self._zooming is None:
            return "break"
        x, y, last = self._zooming
        while e.y <= last - DRAG_LEVEL:
            last -= DRAG_LEVEL
            self.zoom(1, at=(x, y))
        while e.y >= last + DRAG_LEVEL:
            last += DRAG_LEVEL
            self.zoom(-1, at=(x, y))
        self._zooming = (x, y, last)
        return "break"

    def _touchpad(self, e):
        """Tk 9's trackpad scroll: its vertical part zooms, as a wheel's."""
        _dx, dy = gestures.touchpad_delta(self.top, e.delta)
        if dy:
            self._wheel(e, 1 if dy > 0 else -1)

    def _move(self, e):
        self._moved = True
        if self._shape is not None:
            self._rubber(e)
            return
        if self._drag is None or self.view is None:
            return
        self.view.pan(e.x - self._drag[0], e.y - self._drag[1])
        self._auto_fit = False
        self._drag = (e.x, e.y)
        self.draw()

    def _release(self, e):
        """A rectangle or circle drawn, from where the press was; in Look, a
        click that did not pan chooses the region under it."""
        if self._shape is not None:
            x0, y0 = self._shape
            x, y = self.view.to_detector(e.x, e.y)
            self._shape = None
            made = (
                shapes.rectangle(x0, y0, x, y)
                if self.mode.get() == "rectangle"
                else shapes.circle(x0, y0, x, y)
            )
            if self.regions.add(self._imageset(), made):
                self.saved = False
            self.draw()
            return
        if self._drag is not None and not self._moved:
            x, y = self.view.to_detector(e.x, e.y)
            self.chosen = self.regions.at(self._imageset(), x, y)
            self.draw()
        self._drag = None

    def _rubber(self, e):
        """The shape being dragged, drawn as it goes."""
        self.canvas.delete("rubber")
        x0, y0 = self.view.to_canvas(*self._shape)
        if self.mode.get() == "rectangle":
            self.canvas.create_rectangle(
                x0, y0, e.x, e.y, outline=REGION, dash=(4, 2), tags="rubber"
            )
        else:
            r = ((e.x - x0) ** 2 + (e.y - y0) ** 2) ** 0.5
            self.canvas.create_oval(
                x0 - r,
                y0 - r,
                x0 + r,
                y0 + r,
                outline=REGION,
                dash=(4, 2),
                tags="rubber",
            )

    def _double(self, e):
        if self.mode is not None and self.mode.get() == "polygon":
            self._close_polygon()
        else:
            self.fit()

    def _close_polygon(self):
        made = shapes.polygon(self._vertices)
        self._vertices = []
        if self.regions.add(self._imageset(), made):
            self.saved = False
        self.draw()

    def _abandon(self):
        self._vertices = []
        self._shape = None
        self.chosen = None
        self.draw()

    def _mode_changed(self):
        self._vertices = []
        self._shape = None
        self.draw()

    # --- the regions --------------------------------------------------------

    def close(self):
        """Closed, asking first where regions were drawn and not saved."""
        from tkinter import messagebox

        if (
            self.regions.changed
            and not self.saved
            and not messagebox.askyesno(
                "mxi ui",
                "The untrusted regions drawn are not saved. Close anyway?",
                parent=self.top,
            )
        ):
            return
        self.top.destroy()

    def _imageset(self):
        return self.sweeps[self.sweep.get()]["index"]

    def delete_chosen(self):
        if self.chosen is not None and self.regions.remove(
            self._imageset(), self.chosen
        ):
            self.saved = False
        self.chosen = None
        self.draw()

    def undo(self):
        self.regions.undo(self._imageset())
        self.chosen = None
        self.draw()

    def clear(self):
        if self.regions.clear(self._imageset()):
            self.saved = False
        self.chosen = None
        self.draw()

    def same_as_previous(self):
        if self._previous_sweep is None:
            self.readout.set("No other sweep shown yet")
            return
        other = self.sweeps[self._previous_sweep]["index"]
        if self.regions.copy_from(self._imageset(), other):
            self.saved = False
        self.draw()

    def save(self):
        """Every sweep's regions into the output, masked.rflx, beside the
        file shown (mxi ui) or where mxi view runs, through mxi untrusted,
        and mxi_mask.json; the file shown as it was."""
        where = self.where
        try:
            report = mxi_tools.save_untrusted(
                self.path, self.regions.by_imageset, self.output, where
            )
        except (RuntimeError, OSError, SystemExit) as e:
            # SystemExit too: mxi_tools says so where mxi cannot be found,
            # which here would close the interface.
            self.readout.set(
                f"Not saved: {str(e).strip().splitlines()[-1] if str(e).strip() else e}"
            )
            return None
        self.saved = True
        self.regions.changed = False
        last = report.strip().splitlines()[-1] if report.strip() else "saved"
        self.readout.set(f"{last}; and mxi_mask.json")
        masked = (
            self.output
            if os.path.isabs(self.output)
            else os.path.join(where, self.output)
        )
        if self.on_saved is not None:
            self.on_saved(masked)
        self.draw()  # the count no longer "not saved"
        return masked

    def _draw_regions(self):
        for i, region in enumerate(self.regions.of(self._imageset())):
            style = dict(
                outline=CHOSEN if i == self.chosen else REGION,
                width=3 if i == self.chosen else 2,
                tags="region",
            )
            if "rectangle" in region:
                x0, y0, x1, y1 = region["rectangle"]
                self.canvas.create_rectangle(
                    *self.view.to_canvas(x0, y0), *self.view.to_canvas(x1, y1), **style
                )
            elif "circle" in region:
                cx, cy, r = region["circle"]
                (a, b), q = self.view.to_canvas(cx, cy), r / self.view.scale
                self.canvas.create_oval(a - q, b - q, a + q, b + q, **style)
            elif "polygon" in region:
                points = [c for p in region["polygon"] for c in self.view.to_canvas(*p)]
                self.canvas.create_polygon(*points, fill="", **style)
        if self._vertices:
            points = [c for p in self._vertices for c in self.view.to_canvas(*p)]
            if len(points) >= 4:
                self.canvas.create_line(
                    *points, fill=REGION, dash=(4, 2), tags="region"
                )
            for x, y in zip(points[::2], points[1::2]):
                self.canvas.create_oval(
                    x - 3, y - 3, x + 3, y + 3, outline=REGION, tags="region"
                )
        n = len(self.regions.of(self._imageset()))
        unsaved = self.regions.changed and not self.saved
        self.count.configure(
            text=f"Untrusted regions: {n}" + (", not saved" if unsaved else "")
        )

    def zoom(self, direction, at=None, at_pointer=False):
        """A level in (direction > 0) or out, about canvas point `at`; at the
        pointer where asked and it is over the image; else the middle."""
        if self.view is None:
            return
        if at is None and at_pointer:
            px = self.canvas.winfo_pointerx() - self.canvas.winfo_rootx()
            py = self.canvas.winfo_pointery() - self.canvas.winfo_rooty()
            w, h = self._canvas_size()
            if 0 <= px < w and 0 <= py < h:
                at = (px, py)
        if at is None:
            w, h = self._canvas_size()
            at = (w / 2, h / 2)
        self._auto_fit = False
        self.view.zoom(
            1.25 if direction > 0 else 0.8,
            at[0],
            at[1],
            largest=4 * max(self.view.width, self.view.height),
        )
        self.draw()

    def _wheel(self, e, direction, now=None):
        """A wheel's notch, or a trackpad's scroll: one level a tenth of a
        second at most. A trackpad sends many small scroll events, and at a
        level each the view had raced through every zoom."""
        now = time.monotonic() if now is None else now
        if now - self._last_wheel < WHEEL_INTERVAL:
            return False
        self._last_wheel = now
        self.zoom(direction, at=(e.x, e.y))
        return True

    def _hover(self, e):
        if self.view is None or self.frame is None:
            return
        x, y = self.view.to_detector(e.x, e.y)
        h, w = self.frame.shape
        if not (0 <= x < w and 0 <= y < h):
            self.readout.set("")
            return
        counts = self.frame[int(y), int(x)]
        experiments = self.sweeps[self.sweep.get()]["experiments"]
        d = (
            render.resolution(self.document, experiments[0], x, y)
            if experiments
            else None
        )
        self.readout.set(
            f"{self._which()}   x {x:8.1f}   y {y:8.1f}   counts {counts:8.0f}"
            + (f"   d {d:6.2f} Å" if d else "")
        )
