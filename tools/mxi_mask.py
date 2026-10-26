#!/usr/bin/env python3
"""mxi_mask: draw the detector's untrusted regions -- a backstop, its arm, a
shadow -- over a sum of each sweep's images, and write them into the experiment
list, where spot finding and integration leave them out as they leave out the
detector's own bad pixels.

    mxi_mask imported.rflx

One sweep at a time: circles, rectangles and polygons; Undo, Clear, Same as
last (the previous sweep's shapes, or the last saved mask's); Previous and
Next; Save, which writes every sweep's regions into the file (through
mxi_untrusted) and the last mask beside it, mxi_mask.json, to apply again --
`mxi_untrusted other.rflx --set mxi_mask.json` -- with no window.

The image is a sum of images, so that a shadow shows against the background:
the first image and how many, and the contrast, are set at the bottom. Pixels
as DIALS counts them: pixel i spans i to i + 1. While the toolbar's zoom or pan
is on, clicks are the toolbar's.

Needs numpy, h5py, hdf5plugin (for the images' compression) and matplotlib.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

VERSION = "@PROJECT_VERSION@"
LAST = "mxi_mask.json"


# --- the experiment list and the images -------------------------------------


def program(name: str) -> str:
    """An mxi program beside this one, or on the PATH."""
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), name)
    if os.path.exists(here):
        return here
    found = shutil.which(name)
    if not found:
        raise SystemExit(f"mxi_mask: {name} is neither beside mxi_mask nor on the PATH")
    return found


def read_experiments(path: str) -> dict:
    """The experiment list as JSON: an .expt read as it is, a .rflx through
    mxi_convert, which knows its encoding."""
    if path.endswith(".expt"):
        with open(path) as f:
            return json.load(f)
    with tempfile.TemporaryDirectory() as tmp:
        expt = os.path.join(tmp, "list.expt")
        run = subprocess.run(
            [program("mxi_convert"), path, "--output-expt", expt],
            capture_output=True,
            text=True,
            cwd=tmp,
        )
        if run.returncode != 0 or not os.path.exists(expt):
            raise SystemExit(f"mxi_mask: cannot read {path}:\n{run.stderr}{run.stdout}")
        with open(expt) as f:
            return json.load(f)


def sweeps(document: dict) -> list[dict]:
    """Each image set: its index, its master file, its images (first, last)
    and its untrusted regions so far."""
    out = []
    for k, s in enumerate(document.get("imageset", [])):
        first, last = 1, None
        for e in document.get("experiment", []):
            if e.get("imageset") == k and e.get("scan") is not None:
                scan = document["scan"][e["scan"]]
                first, last = scan["image_range"]
                break
        out.append(
            {
                "index": k,
                "master": s.get("template", ""),
                "first": int(first),
                "last": int(last) if last is not None else None,
                "shapes": list(s.get("untrusted", [])),
            }
        )
    return out


class Images:
    """An NXmx series' frames, summed: the master's /entry/data, one dataset
    or several (data_000001, ...), each read as it is stored."""

    def __init__(self, master: str):
        import h5py
        import numpy as np

        try:
            import hdf5plugin  # noqa: F401 -- bitshuffle and LZ4, registered
        except ImportError:
            pass
        self.np = np
        self.file = h5py.File(master, "r")
        data = self.file["/entry/data"]
        names = sorted(n for n in data if n == "data" or n.startswith("data_"))
        self.blocks = []  # (dataset, first frame in it)
        start = 0
        for n in names:
            try:
                d = data[n]
                count = d.shape[0]
            except (KeyError, OSError):
                continue  # a link to a file not there
            self.blocks.append((d, start))
            start += count
        if not self.blocks:
            raise SystemExit(f"mxi_mask: {master} has no frames under /entry/data")
        self.frames = start
        self.shape = self.blocks[0][0].shape[1:]
        self.bad = np.zeros(self.shape, bool)
        for where in (
            "/entry/instrument/detector/pixel_mask",
            "/entry/instrument/detector/detectorSpecific/pixel_mask",
        ):
            if where in self.file:
                m = self.file[where][()]
                if m.shape == self.shape:
                    self.bad |= m != 0
                break

    def sum(self, first: int, count: int):
        """Frames first .. first + count - 1, counting from 0, summed; bad and
        overloaded-marker pixels as 0."""
        np = self.np
        total = np.zeros(self.shape, np.float64)
        for i in range(max(0, first), min(self.frames, first + max(1, count))):
            for d, start in reversed(self.blocks):
                if i >= start:
                    frame = d[i - start].astype(np.float64)
                    break
            limit = (
                float(np.iinfo(d.dtype).max) - 2 if d.dtype.kind in "ui" else math.inf
            )
            frame[(frame >= limit) | (frame < 0)] = 0
            total += frame
        total[self.bad] = 0
        return total


# --- the window ---------------------------------------------------------------


class Editor:
    def __init__(self, path: str, document: dict, start_sum: int):
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Button, RadioButtons, Slider, TextBox

        self.plt = plt
        self.path = path
        self.document = document
        self.sweeps = sweeps(document)
        if not self.sweeps:
            raise SystemExit(f"mxi_mask: {path} has no image sets")
        self.current = 0
        self.history = {}  # sweep -> list of shapes as they were before each change
        self.last = self.load_last()
        self.previous_edited = None
        self.images = {}
        self.tool = "rectangle"
        self.press = None
        self.polygon = []
        self.preview = None
        self.count = start_sum
        self.first = 0

        self.fig = plt.figure(figsize=(11, 9))
        self.ax = self.fig.add_axes([0.05, 0.22, 0.75, 0.73])
        self.status = self.fig.text(0.05, 0.965, "", fontsize=10)

        def button(rect, label, action):
            b = Button(self.fig.add_axes(rect), label)
            b.on_clicked(lambda _: action())
            return b

        self.buttons = [
            button([0.83, 0.86, 0.14, 0.05], "Undo", self.undo),
            button([0.83, 0.80, 0.14, 0.05], "Clear", self.clear),
            button([0.83, 0.74, 0.14, 0.05], "Same as last", self.same_as_last),
            button([0.83, 0.62, 0.065, 0.05], "Previous", lambda: self.go(-1)),
            button([0.905, 0.62, 0.065, 0.05], "Next", lambda: self.go(+1)),
            button([0.83, 0.50, 0.14, 0.06], "Save", self.save),
        ]
        self.radio = RadioButtons(
            self.fig.add_axes([0.83, 0.28, 0.14, 0.14]),
            ("rectangle", "circle", "polygon"),
        )
        self.radio.on_clicked(self.choose)
        self.first_box = TextBox(
            self.fig.add_axes([0.16, 0.12, 0.1, 0.04]), "first image ", "1"
        )
        self.first_box.on_submit(self.set_first)
        self.count_box = TextBox(
            self.fig.add_axes([0.42, 0.12, 0.1, 0.04]), "images summed ", str(start_sum)
        )
        self.count_box.on_submit(self.set_count)
        self.contrast = Slider(
            self.fig.add_axes([0.16, 0.06, 0.5, 0.03]),
            "contrast",
            90.0,
            99.99,
            valinit=99.5,
        )
        self.contrast.on_changed(lambda _: self.redraw_image())
        self.fig.text(
            0.05,
            0.015,
            "rectangle: drag from a corner; circle: drag from the centre; polygon: click "
            "the vertices, right-click or Enter to close",
            fontsize=8,
        )
        self.fig.canvas.mpl_connect("button_press_event", self.on_press)
        self.fig.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.fig.canvas.mpl_connect("button_release_event", self.on_release)
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.image_artist = None
        self.patches = []
        self.show_sweep()

    # The last saved mask, for Same as last on the first sweep.
    def load_last(self):
        if os.path.exists(LAST):
            try:
                with open(LAST) as f:
                    saved = json.load(f)
                return saved.get("last", []) if isinstance(saved, dict) else saved
            except (OSError, ValueError):
                return []
        return []

    @property
    def sweep(self):
        return self.sweeps[self.current]

    def remember(self):
        self.history.setdefault(self.current, []).append(list(self.sweep["shapes"]))

    def images_of(self, sweep):
        master = sweep["master"]
        if master not in self.images:
            self.images[master] = Images(master)
        return self.images[master]

    # --- drawing the image and the shapes
    def show_sweep(self):
        # Each sweep from its own first image: image n of the file is its
        # frame n - 1, whatever image the scan begins at.
        self.polygon = []
        self.first = max(0, self.sweep["first"] - 1)
        self.first_box.set_val(str(self.first + 1))
        self.redraw_image(reset=True)

    def redraw_image(self, reset=False):
        import numpy as np

        s = self.sweep
        images = self.images_of(s)
        total = images.sum(self.first, self.count)
        shown = np.sqrt(np.clip(total, 0, None))  # a shadow against background
        top = (
            np.percentile(shown[shown > 0], self.contrast.val)
            if np.any(shown > 0)
            else 1
        )
        h, w = total.shape
        if self.image_artist is None or reset:
            self.ax.clear()
            self.image_artist = self.ax.imshow(
                shown,
                cmap="gray_r",
                vmin=0,
                vmax=top,
                extent=(0, w, h, 0),
                interpolation="nearest",
            )
            self.patches = []
        else:
            self.image_artist.set_data(shown)
            self.image_artist.set_clim(0, top)
        self.draw_shapes()

    def draw_shapes(self):
        from matplotlib.patches import Circle, Polygon, Rectangle

        for p in self.patches:
            p.remove()
        self.patches = []
        style = dict(fill=True, facecolor=(1, 0, 0, 0.25), edgecolor="red", lw=1.5)
        for shape in self.sweep["shapes"]:
            if "circle" in shape:
                x, y, r = shape["circle"]
                p = Circle((x, y), r, **style)
            elif "rectangle" in shape:
                x0, y0, x1, y1 = shape["rectangle"]
                p = Rectangle(
                    (min(x0, x1), min(y0, y1)), abs(x1 - x0), abs(y1 - y0), **style
                )
            elif "polygon" in shape:
                p = Polygon(shape["polygon"], closed=True, **style)
            else:
                continue
            self.ax.add_patch(p)
            self.patches.append(p)
        if self.polygon:
            xs, ys = zip(*self.polygon)
            (line,) = self.ax.plot(xs, ys, "r.-", lw=1)
            self.patches.append(line)
        s = self.sweep
        last = f"{s['last']}" if s["last"] else "?"
        self.status.set_text(
            f"sweep {self.current + 1} of {len(self.sweeps)}: "
            f"{os.path.basename(s['master'])}, images {s['first']}-{last}; summed "
            f"{self.count} from image {self.first + 1}; "
            f"{len(s['shapes'])} regions; drawing {self.tool}s"
        )
        self.fig.canvas.draw_idle()

    # --- the controls
    def choose(self, label):
        self.tool = label
        self.polygon = []
        self.draw_shapes()

    def set_first(self, text):
        try:
            self.first = max(0, int(text) - 1)
        except ValueError:
            return
        self.redraw_image()

    def set_count(self, text):
        try:
            self.count = max(1, int(text))
        except ValueError:
            return
        self.redraw_image()

    def undo(self):
        if self.polygon:
            self.polygon = []
        elif self.history.get(self.current):
            self.sweep["shapes"] = self.history[self.current].pop()
        self.draw_shapes()

    def clear(self):
        if self.sweep["shapes"]:
            self.remember()
            self.sweep["shapes"] = []
        self.polygon = []
        self.draw_shapes()

    def same_as_last(self):
        source = (
            self.sweeps[self.previous_edited]["shapes"]
            if self.previous_edited is not None and self.previous_edited != self.current
            else self.last
        )
        if not source:
            self.status.set_text(
                "no last mask: nothing drawn on another sweep, none saved"
            )
            self.fig.canvas.draw_idle()
            return
        self.remember()
        self.sweep["shapes"] = json.loads(json.dumps(source))
        self.draw_shapes()

    def go(self, step):
        if self.sweep["shapes"]:
            self.previous_edited = self.current
        self.current = (self.current + step) % len(self.sweeps)
        self.show_sweep()

    def save(self):
        mask = {"imagesets": {str(s["index"]): s["shapes"] for s in self.sweeps}}
        last = self.sweep["shapes"] or (
            self.sweeps[self.previous_edited]["shapes"]
            if self.previous_edited is not None
            else []
        )
        with open(LAST, "w") as f:
            json.dump({"imagesets": mask["imagesets"], "last": last}, f, indent=1)
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
            json.dump(mask, f)
            name = f.name
        try:
            run = subprocess.run(
                [program("mxi_untrusted"), self.path, "--set", name],
                capture_output=True,
                text=True,
            )
        finally:
            os.unlink(name)
        message = (
            run.stdout.strip().splitlines()[-1] if run.returncode == 0 else run.stderr
        )
        self.status.set_text(f"{message}; and {LAST}")
        print(run.stdout + run.stderr, end="")
        self.fig.canvas.draw_idle()

    # --- the mouse and the keys
    def busy(self):
        toolbar = getattr(self.fig.canvas, "toolbar", None)
        return toolbar is not None and getattr(toolbar, "mode", "") not in ("", None)

    def on_press(self, event):
        if event.inaxes is not self.ax or self.busy() or event.xdata is None:
            return
        x, y = float(event.xdata), float(event.ydata)
        if self.tool == "polygon":
            if event.button == 3 or event.dblclick:
                self.close_polygon()
            else:
                self.polygon.append([round(x, 1), round(y, 1)])
                self.draw_shapes()
            return
        self.press = (x, y)

    def on_motion(self, event):
        if self.press is None or event.inaxes is not self.ax or event.xdata is None:
            return
        from matplotlib.patches import Circle, Rectangle

        if self.preview is not None:
            self.preview.remove()
        x0, y0 = self.press
        x, y = float(event.xdata), float(event.ydata)
        if self.tool == "circle":
            self.preview = Circle(
                (x0, y0), math.hypot(x - x0, y - y0), fill=False, ec="red", ls="--"
            )
        else:
            self.preview = Rectangle(
                (min(x0, x), min(y0, y)),
                abs(x - x0),
                abs(y - y0),
                fill=False,
                ec="red",
                ls="--",
            )
        self.ax.add_patch(self.preview)
        self.fig.canvas.draw_idle()

    def on_release(self, event):
        if self.press is None:
            return
        if self.preview is not None:
            self.preview.remove()
            self.preview = None
        x0, y0 = self.press
        self.press = None
        if event.xdata is None:
            self.draw_shapes()
            return
        x, y = float(event.xdata), float(event.ydata)
        if self.tool == "circle":
            r = math.hypot(x - x0, y - y0)
            if r >= 1:
                self.remember()
                self.sweep["shapes"].append(
                    {"panel": 0, "circle": [round(x0, 1), round(y0, 1), round(r, 1)]}
                )
        elif abs(x - x0) >= 1 and abs(y - y0) >= 1:
            self.remember()
            self.sweep["shapes"].append(
                {
                    "panel": 0,
                    "rectangle": [
                        round(min(x0, x), 1),
                        round(min(y0, y), 1),
                        round(max(x0, x), 1),
                        round(max(y0, y), 1),
                    ],
                }
            )
        self.draw_shapes()

    def on_key(self, event):
        if event.key == "enter" and self.tool == "polygon":
            self.close_polygon()
        elif event.key in ("ctrl+z", "cmd+z", "u"):
            self.undo()

    def close_polygon(self):
        if len(self.polygon) >= 3:
            self.remember()
            self.sweep["shapes"].append({"panel": 0, "polygon": self.polygon})
        self.polygon = []
        self.draw_shapes()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="mxi_mask",
        description="Draw the detector's untrusted regions over a sum of each sweep's "
        "images, and write them into the experiment list (through mxi_untrusted), "
        f"where spot finding and integration leave them out; the last mask kept in "
        f"{LAST} to apply again: mxi_untrusted other.rflx --set {LAST}.",
    )
    parser.add_argument(
        "experiments", help="the experiment list, as mxi_import writes it"
    )
    parser.add_argument(
        "--sum", type=int, default=10, help="images summed at first (10)"
    )
    parser.add_argument(
        "--snapshot",
        default=None,
        help="no window: the first sweep's sum and regions drawn to this PNG",
    )
    parser.add_argument("--version", action="version", version=f"mxi_mask {VERSION}")
    args = parser.parse_args(argv)

    document = read_experiments(args.experiments)
    if args.snapshot:
        import matplotlib

        matplotlib.use("Agg")
    editor = Editor(args.experiments, document, args.sum)
    if args.snapshot:
        editor.fig.savefig(args.snapshot, dpi=90)
        print(f"Wrote {args.snapshot}")
        return 0
    editor.plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
