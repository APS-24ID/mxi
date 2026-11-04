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
import tkinter as tk
from tkinter import ttk

import numpy as np

import mxi_tools

from mxi_ui import render
from mxi_ui.tooltip import Tooltip

RADIUS = 7  # canvas pixels, a marker's
LAYERS = (
    # key, label, the table's column, colour
    ("indexed", "Indexed spots", "xyzobs.px.value", "#2e7d32"),
    ("unindexed", "Spots", "xyzobs.px.value", "#d32f2f"),
    ("predicted", "Predictions", "xyzcal.px", "#1565c0"),
)
COLUMNS = ["id", "xyzobs.px.value", "xyzcal.px", "bbox", "miller_index"]


def masters(document, path):
    """Each sweep (mxi_tools.sweeps), its master made absolute: a relative one
    is the file's directory's."""
    here = os.path.dirname(os.path.abspath(path))
    out = mxi_tools.sweeps(document)
    for s in out:
        if s["master"] and not os.path.isabs(s["master"]):
            s["master"] = os.path.join(here, s["master"])
    return out


def layers_for(table, experiments):
    """The marker layers a sweep's rows make: spots indexed and not where the
    table has indices -- else all its spots as one -- and predictions where it
    has them; each a mxi_tools.Markers, None where it would be empty."""
    if "id" not in table:
        return {}
    ids = np.asarray(table["id"])
    mine = np.isin(ids, list(experiments))
    out = {}
    if "xyzobs.px.value" in table:
        if "miller_index" in table:
            hkl = np.asarray(table["miller_index"]).reshape(len(ids), -1)
            indexed = mine & np.any(hkl != 0, axis=1)
            groups = {"indexed": indexed, "unindexed": mine & ~indexed}
        else:
            groups = {"unindexed": mine}
        for key, rows in groups.items():
            sub = {k: np.asarray(v)[rows] for k, v in table.items()}
            markers = mxi_tools.Markers(sub, "xyzobs.px.value", experiments)
            out[key] = markers if len(markers) else None
    if "xyzcal.px" in table:
        sub = {k: np.asarray(v)[mine] for k, v in table.items()}
        markers = mxi_tools.Markers(sub, "xyzcal.px", experiments)
        out["predicted"] = markers if len(markers) else None
    return out


class Viewer:
    def __init__(self, parent, path, title=None):
        self.path = path
        self.document = mxi_tools.read_experiments(path)
        self.sweeps = masters(self.document, path)
        if not self.sweeps:
            raise ValueError(f"{os.path.basename(path)} names no images")
        self.table = (
            mxi_tools.read_reflections(path, COLUMNS) if path.endswith(".rflx") else {}
        )
        self.images = {}
        self.layers = {}
        self.top = tk.Toplevel(parent)
        self.top.title(title or f"mxi ui: {os.path.basename(path)}")
        self.top.geometry("1000x900")
        self.sweep = tk.IntVar(value=0)
        self.image = tk.IntVar(value=1)
        self.ceiling = tk.DoubleVar(value=10.0)
        self.scheme = tk.StringVar(value="grey")
        self.showing = {key: tk.BooleanVar(value=True) for key, *_ in LAYERS}
        self.readout = tk.StringVar(value="")
        self.view = None
        self.photo = None
        self.frame = None
        self._drag = None
        self._reduced = None
        self._controls()
        self.canvas = tk.Canvas(self.top, background="#808080", highlightthickness=0)
        self.canvas.pack(fill=tk.BOTH, expand=True)
        ttk.Label(
            self.top, textvariable=self.readout, padding=(8, 3), font="TkFixedFont"
        ).pack(fill=tk.X)
        self._bindings()
        self.choose_sweep(0)

    # --- building ---------------------------------------------------------

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
        self.of.pack(side=tk.LEFT, padx=(4, 12))
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
        ttk.Button(bar, text="Fit", width=4, command=self.fit).pack(side=tk.LEFT)
        self.toggles = ttk.Frame(self.top, padding=(8, 0, 8, 4))
        self.toggles.pack(fill=tk.X)

    def _bindings(self):
        c = self.canvas
        c.bind("<Configure>", lambda _e: self.draw())
        c.bind("<ButtonPress-1>", self._press)
        c.bind("<B1-Motion>", self._move)
        c.bind("<Double-Button-1>", lambda _e: self.fit())
        c.bind("<Motion>", self._hover)
        c.bind("<MouseWheel>", lambda e: self._wheel(e, 1 if e.delta > 0 else -1))
        c.bind("<Button-4>", lambda e: self._wheel(e, 1))  # X11's wheel
        c.bind("<Button-5>", lambda e: self._wheel(e, -1))
        self.top.bind("<Left>", lambda _e: self.step(-1))
        self.top.bind("<Right>", lambda _e: self.step(1))

    # --- the sweep and image ----------------------------------------------

    def choose_sweep(self, k):
        self.sweep.set(k)
        s = self.sweeps[k]
        if s["master"] not in self.images:
            self.images[s["master"]] = mxi_tools.Images(s["master"])
        images = self.images[s["master"]]
        last = s["last"] if s["last"] is not None else len(images)
        self.first, self.last = s["first"], last
        self.spin.configure(from_=self.first, to=self.last)
        self.of.configure(text=f"of {self.first} to {self.last}")
        if k not in self.layers:
            self.layers[k] = layers_for(self.table, s["experiments"])
        for child in self.toggles.winfo_children():
            child.destroy()
        for key, label, _column, colour in LAYERS:
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
        self.image.set(self.first)
        self.show(fit_contrast=True)

    def step(self, by):
        self.image.set(min(self.last, max(self.first, self._image() + by)))
        self.show()

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
        s = self.sweeps[self.sweep.get()]
        self.frame = self.images[s["master"]].frame(n - 1)
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

    def _markers(self):
        layers = self.layers.get(self.sweep.get(), {})
        frame = self._image() - 1
        width, height = self._canvas_size()
        for key, _label, _column, colour in LAYERS:
            markers = layers.get(key)
            if markers is None or not self.showing[key].get():
                continue
            yx, centred = markers.on(frame)
            for (y, x), own in zip(yx, centred):
                # Markers are napari's, half a pixel less than mxi's
                cx, cy = self.view.to_canvas(x + 0.5, y + 0.5)
                if -RADIUS <= cx <= width + RADIUS and -RADIUS <= cy <= height + RADIUS:
                    self.canvas.create_oval(
                        cx - RADIUS,
                        cy - RADIUS,
                        cx + RADIUS,
                        cy + RADIUS,
                        outline=colour,
                        width=2 if own else 1,
                        dash=() if own else (2, 2),
                        tags=key,
                    )

    # --- the mouse ----------------------------------------------------------

    def _press(self, e):
        self._drag = (e.x, e.y)

    def _move(self, e):
        if self._drag is None or self.view is None:
            return
        self.view.pan(e.x - self._drag[0], e.y - self._drag[1])
        self._drag = (e.x, e.y)
        self.draw()

    def _wheel(self, e, direction):
        if self.view is None:
            return
        self.view.zoom(
            1.25 if direction > 0 else 0.8,
            e.x,
            e.y,
            largest=4 * max(self.view.width, self.view.height),
        )
        self.draw()

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
            f"image {self._image()}   x {x:8.1f}   y {y:8.1f}   counts {counts:8.0f}"
            + (f"   d {d:6.2f} Å" if d else "")
        )
