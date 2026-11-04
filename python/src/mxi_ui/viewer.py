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
import time
import tkinter as tk
from tkinter import ttk

import numpy as np

import mxi_tools

from mxi_ui import render
from mxi_ui.tooltip import Tooltip

MIN_RADIUS = 5  # canvas pixels: the least a marker is, however small its spot
MARGIN = 3  # canvas pixels between a spot's box and its marker, to see the spot
DRAG_LEVEL = 30  # canvas pixels dragged a zoom level, with the right button
WHEEL_INTERVAL = 0.1  # seconds: at most a zoom level each, from a scroll
LAYERS = (
    # key, label, the table's column, colour
    ("indexed", "Indexed spots", "xyzobs.px.value", "#2e7d32"),
    ("unindexed", "Spots", "xyzobs.px.value", "#d32f2f"),
    ("predicted", "Predictions", "xyzcal.px", "#1565c0"),
)
COLUMNS = ["id", "xyzobs.px.value", "xyzcal.px", "bbox", "miller_index"]


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
        self._last_wheel = -1e9
        self._zooming = None
        self._auto_fit = True
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
        # Buttons and keys to zoom, not a wheel alone: there need be no wheel
        # (Graeme, a trackpad).
        zoom_in = ttk.Button(bar, text="+", width=3, command=lambda: self.zoom(1))
        zoom_in.pack(side=tk.LEFT)
        zoom_out = ttk.Button(bar, text="−", width=3, command=lambda: self.zoom(-1))
        zoom_out.pack(side=tk.LEFT, padx=(2, 4))
        fit = ttk.Button(bar, text="Fit", width=4, command=self.fit)
        fit.pack(side=tk.LEFT)
        Tooltip(
            zoom_in,
            "Zoom in: + or =; or drag up with the right button (or Control-drag)",
        )
        Tooltip(zoom_out, "Zoom out: -; or drag down with the right button")
        Tooltip(fit, "The whole image: 0, or double-click")
        self.toggles = ttk.Frame(self.top, padding=(8, 0, 8, 4))
        self.toggles.pack(fill=tk.X)

    def _bindings(self):
        c = self.canvas
        c.bind("<Configure>", lambda _e: self._resized())
        c.bind("<ButtonPress-1>", self._press)
        c.bind("<B1-Motion>", self._move)
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
        c.bind("<Double-Button-1>", lambda _e: self.fit())
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
        for key, action in keys.items():
            self.top.bind(key, lambda e, act=action: None if typing(e) else act())

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
        self._auto_fit = True
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

    def _markers(self):
        """Each marker the size of its spot on screen -- its box's half-width
        and half-height at the view's scale -- and never less than MIN_RADIUS,
        so that a spot shows zoomed out: drawn all one size, as at first, a
        marker said nothing of its spot (Graeme)."""
        layers = self.layers.get(self.sweep.get(), {})
        frame = self._image() - 1
        width, height = self._canvas_size()
        scale = self.view.scale
        for key, _label, _column, colour in LAYERS:
            markers = layers.get(key)
            if markers is None or not self.showing[key].get():
                continue
            yx, centred, half = markers.on_sized(frame)
            for (y, x), own, (hx, hy) in zip(yx, centred, half):
                # Markers are napari's, half a pixel less than mxi's
                cx, cy = self.view.to_canvas(x + 0.5, y + 0.5)
                rx = max(MIN_RADIUS, hx / scale + MARGIN)
                ry = max(MIN_RADIUS, hy / scale + MARGIN)
                if -rx <= cx <= width + rx and -ry <= cy <= height + ry:
                    self.canvas.create_oval(
                        cx - rx,
                        cy - ry,
                        cx + rx,
                        cy + ry,
                        outline=colour,
                        width=2 if own else 1,
                        dash=() if own else (2, 2),
                        tags=key,
                    )

    # --- the mouse ----------------------------------------------------------

    def _press(self, e):
        self._drag = (e.x, e.y)
        self.canvas.focus_set()  # the keys the image's, not the image box's

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
        try:
            _dx, dy = (
                int(v)
                for v in self.top.tk.splitlist(
                    self.top.tk.call("tk::PreciseScrollDeltas", e.delta)
                )
            )
        except (tk.TclError, ValueError):
            dy = e.delta & 0xFFFF
            dy = dy - 0x10000 if dy >= 0x8000 else dy
        if dy:
            self._wheel(e, 1 if dy > 0 else -1)

    def _move(self, e):
        if self._drag is None or self.view is None:
            return
        self.view.pan(e.x - self._drag[0], e.y - self._drag[1])
        self._auto_fit = False
        self._drag = (e.x, e.y)
        self.draw()

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
            f"image {self._image()}   x {x:8.1f}   y {y:8.1f}   counts {counts:8.0f}"
            + (f"   d {d:6.2f} Å" if d else "")
        )
