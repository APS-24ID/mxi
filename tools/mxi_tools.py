"""What mxi's Python tools share: finding mxi's programs, reading an experiment
list, its reflections and its images, and the regions a user marks untrusted.

Copied beside the programs by the build, and found there by the tools, which
put their own directory on sys.path. Independent of mxeq, the referee, as mxi
must be: it reads the files itself.

Pixels, two conventions: DIALS's, which mxi's files use, where pixel i spans i
to i + 1 and a spot's centroid at x lies in pixel floor(x); and napari's,
where pixel i is centred on i. So napari = DIALS - 0.5, in x, y and z alike --
image n, frame n - 1, spanning z from n - 1 to n, is napari's slice n - 1.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import tempfile

import numpy as np

# --- mxi's programs ----------------------------------------------------------


def program(name: str) -> str:
    """An mxi program beside this module, or on the PATH."""
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), name)
    if os.path.exists(here):
        return here
    found = shutil.which(name)
    if not found:
        raise SystemExit(f"{name} is neither beside the mxi tools nor on the PATH")
    return found


# --- the experiment list -----------------------------------------------------


def read_experiments(path: str) -> dict:
    """The experiment list as JSON: an .expt as it is, a .rflx through
    mxi_convert, which knows its encoding -- run in a scratch directory, the
    path made absolute first."""
    path = os.path.abspath(path)
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
            raise SystemExit(f"cannot read {path}:\n{run.stderr}{run.stdout}")
        with open(expt) as f:
            return json.load(f)


def sweeps(document: dict) -> list[dict]:
    """Each image set: its index, its master, its images (first, last), its
    untrusted regions, and the experiments on it."""
    out = []
    for k, s in enumerate(document.get("imageset", [])):
        first, last, members = 1, None, []
        for i, e in enumerate(document.get("experiment", [])):
            if e.get("imageset") != k:
                continue
            members.append(i)
            if e.get("scan") is not None and last is None:
                first, last = document["scan"][e["scan"]]["image_range"]
        out.append(
            {
                "index": k,
                "master": s.get("template", ""),
                "first": int(first),
                "last": int(last) if last is not None else None,
                "shapes": list(s.get("untrusted", [])),
                "experiments": members,
            }
        )
    return out


# --- reflections -------------------------------------------------------------


def read_reflections(path: str, columns=None) -> dict:
    """A .rflx's reflection tables, joined in the order mxi writes them
    (shortest name first, then by name), as {column: array}: the columns
    asked for, or all; an empty dict where the file holds none."""
    import h5py

    out: dict = {}
    with h5py.File(path, "r") as f:
        if "reflections" not in f:
            return out
        group = f["reflections"]
        names = sorted(
            (n for n in group if isinstance(group[n], h5py.Group)),
            key=lambda n: (len(n), n),
        )
        for name in names:
            table = group[name]
            for key, dataset in table.items():
                if not isinstance(dataset, h5py.Dataset):
                    continue
                if columns is not None and key not in columns:
                    continue
                array = dataset[()]
                if array.ndim == 3:
                    array = array.reshape(array.shape[0], -1)
                out.setdefault(key, []).append(array)
    return {k: np.concatenate(v) for k, v in out.items()}


def points_of(table: dict, column: str, experiments) -> np.ndarray:
    """The rows of `experiments` as napari points, (z, y, x) in napari's
    pixels: from a position column in DIALS's, (x, y, z)."""
    if column not in table or "id" not in table:
        return np.zeros((0, 3))
    xyz = np.asarray(table[column], float)
    rows = np.isin(np.asarray(table["id"]), list(experiments)) & np.all(
        np.isfinite(xyz), axis=1
    )
    xyz = xyz[rows]
    rows = np.any(xyz != 0, axis=1)  # a position of zero: none
    xyz = xyz[rows]
    return np.column_stack([xyz[:, 2], xyz[:, 1], xyz[:, 0]]) - 0.5


# --- images ------------------------------------------------------------------


class Images:
    """An NXmx series' frames, read as they are asked for: /entry/data, one
    dataset or several (data_000001, ...), the detector's pixel mask applied
    and overloaded markers zeroed. Looks like an array to napari -- shape,
    dtype, ndim and slicing -- and sums `window` frames from each one asked
    for, so that a shadow shows against the background."""

    def __init__(self, master: str, window: int = 1):
        import h5py

        try:
            import hdf5plugin  # noqa: F401 -- bitshuffle and LZ4, registered
        except ImportError:
            pass
        self.file = h5py.File(master, "r")
        data = self.file["/entry/data"]
        names = sorted(n for n in data if n == "data" or n.startswith("data_"))
        self.blocks = []
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
            raise SystemExit(f"{master} has no frames under /entry/data")
        self.frames = start
        first = self.blocks[0][0]
        self.frame_shape = tuple(first.shape[1:])
        self.limit = (
            float(np.iinfo(first.dtype).max) - 2
            if first.dtype.kind in "ui"
            else math.inf
        )
        self.bad = np.zeros(self.frame_shape, bool)
        for where in (
            "/entry/instrument/detector/pixel_mask",
            "/entry/instrument/detector/detectorSpecific/pixel_mask",
        ):
            if where in self.file:
                m = self.file[where][()]
                if m.shape == self.frame_shape:
                    self.bad |= m != 0
                break
        self.window = max(1, int(window))
        self._cache: dict = {}

    # What napari asks of an array.
    @property
    def shape(self):
        return (self.frames,) + self.frame_shape

    @property
    def dtype(self):
        return np.dtype(np.float32)

    @property
    def ndim(self):
        return 3

    def __len__(self):
        return self.frames

    def frame(self, i: int) -> np.ndarray:
        """Frame i, from 0, its bad and overloaded pixels zero."""
        for d, start in reversed(self.blocks):
            if i >= start:
                f = d[i - start].astype(np.float32)
                break
        f[(f >= self.limit) | (f < 0)] = 0
        f[self.bad] = 0
        return f

    def summed(self, i: int) -> np.ndarray:
        """Frames i to i + window - 1 summed; a few kept, as napari asks for
        the same slice again as it draws."""
        key = (i, self.window)
        if key not in self._cache:
            if len(self._cache) > 8:
                self._cache.clear()
            total = np.zeros(self.frame_shape, np.float32)
            for k in range(max(0, i), min(self.frames, i + self.window)):
                total += self.frame(k)
            self._cache[key] = total
        return self._cache[key]

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        first, rest = key[0], key[1:]
        if isinstance(first, (int, np.integer)):
            return self.summed(int(first))[rest] if rest else self.summed(int(first))
        indices = (
            range(*first.indices(self.frames)) if isinstance(first, slice) else first
        )
        stack = np.stack([self.summed(int(i)) for i in indices])
        return stack[(slice(None),) + rest] if rest else stack

    def __array__(self, dtype=None):  # only for small series: napari slices
        return np.asarray(self[:], dtype=dtype)


def display_top(frame: np.ndarray) -> float:
    """Where a diffraction image's grey scale tops out: well above its
    background, so that the background shows light and the spots dark -- the
    largest of the 99.9th percentile, five times the median and ten counts. A
    diffraction image is nearly all background: the 99.5th percentile of its
    counting pixels was the background itself, and every pixel at or above it
    black (23 October)."""
    counting = frame[frame > 0]
    if counting.size == 0:
        return 10.0
    return float(max(np.percentile(counting, 99.9), 5.0 * np.median(counting), 10.0))


# --- untrusted regions -------------------------------------------------------


def shape_to_untrusted(kind: str, vertices) -> dict | None:
    """A napari shape -- its type and its vertices, (y, x) in napari's pixels
    (a shape drawn over a stack may carry the slice first, which is dropped:
    a region is the whole sweep's) -- as an untrusted region in DIALS's
    pixels; None for a line or a path, which mark no area."""
    v = np.asarray(vertices, float)[:, -2:] + 0.5
    xy = v[:, ::-1]
    if kind == "rectangle":
        x0, y0 = xy.min(axis=0)
        x1, y1 = xy.max(axis=0)
        return {"panel": 0, "rectangle": [round(float(c), 2) for c in (x0, y0, x1, y1)]}
    if kind == "ellipse":
        x0, y0 = xy.min(axis=0)
        x1, y1 = xy.max(axis=0)
        cx, cy, rx, ry = (x0 + x1) / 2, (y0 + y1) / 2, (x1 - x0) / 2, (y1 - y0) / 2
        if abs(rx - ry) <= 0.01 * max(rx, ry, 1e-9):
            return {
                "panel": 0,
                "circle": [round(float(c), 2) for c in (cx, cy, (rx + ry) / 2)],
            }
        t = np.linspace(0, 2 * math.pi, 48, endpoint=False)
        return {
            "panel": 0,
            "polygon": [
                [
                    round(float(cx + rx * math.cos(a)), 2),
                    round(float(cy + ry * math.sin(a)), 2),
                ]
                for a in t
            ],
        }
    if kind == "polygon" and len(xy) >= 3:
        return {
            "panel": 0,
            "polygon": [[round(float(x), 2), round(float(y), 2)] for x, y in xy],
        }
    return None


def untrusted_to_shape(region: dict):
    """An untrusted region as a napari shape: (type, vertices (y, x) in
    napari's pixels)."""
    if "circle" in region:
        x, y, r = region["circle"]
        return (
            "ellipse",
            np.array([[y - r, x - r], [y - r, x + r], [y + r, x + r], [y + r, x - r]])
            - 0.5,
        )
    if "rectangle" in region:
        x0, y0, x1, y1 = region["rectangle"]
        return (
            "rectangle",
            np.array([[y0, x0], [y0, x1], [y1, x1], [y1, x0]], float) - 0.5,
        )
    if "polygon" in region:
        return "polygon", np.array([[y, x] for x, y in region["polygon"]], float) - 0.5
    raise ValueError(f"not an untrusted region: {region}")


def save_untrusted(path: str, by_imageset: dict, output: str = "masked.rflx") -> str:
    """Each image set's regions written by mxi_untrusted into a new file,
    `output`, the input left as it was; and kept in mxi_mask.json, to apply
    again with no window. Returns mxi_untrusted's report."""
    mask = {"imagesets": {str(k): v for k, v in by_imageset.items()}}
    with open("mxi_mask.json", "w") as f:
        json.dump(mask, f, indent=1)
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(mask, f)
        name = f.name
    try:
        run = subprocess.run(
            [
                program("mxi_untrusted"),
                os.path.abspath(path),
                "--set",
                name,
                "-o",
                output,
            ],
            capture_output=True,
            text=True,
        )
    finally:
        os.unlink(name)
    if run.returncode != 0:
        raise RuntimeError(run.stderr or run.stdout)
    return run.stdout
