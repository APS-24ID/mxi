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


# The .rflx's encoding of JSON in HDF5 (docs/rflx.md): an object a group, its
# members attributes or members of it, names escaped; a list a group with
# __type__ "list" and __len__; a null an empty dataspace.
_HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"
_RESERVED = {"__type__", "__keys__", "__len__"}


def _text(value) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def _escape(key: str) -> str:
    if key == "":
        return "%00"
    out = key.replace("%", "%25").replace("/", "%2F")
    return "%5F" + out[1:] if out.startswith("_") else out


def _unescape(name: str) -> str:
    if name == "%00":
        return ""
    if name.startswith("%5F"):
        name = "_" + name[3:]
    return name.replace("%2F", "/").replace("%25", "%")


def _leaf(h5py, value, null: bool):
    if null or isinstance(value, h5py.Empty):
        return None
    if isinstance(value, np.ndarray):
        if value.dtype.kind in "OSU":
            return [_text(v) for v in value.ravel()] if value.size else []
        return value.tolist()
    if isinstance(value, (bytes, str)):
        return _text(value)
    if isinstance(value, np.generic):
        return value.item()
    return value


def _member(h5py, group, name: str):
    if name in group.attrs:
        a = group.attrs.get_id(name)
        return _leaf(
            h5py,
            group.attrs[name],
            a.get_space().get_simple_extent_type() == h5py.h5s.NULL,
        )
    obj = group[name]
    if isinstance(obj, h5py.Group):
        return _decode(h5py, obj)
    null = obj.id.get_space().get_simple_extent_type() == h5py.h5s.NULL
    return _leaf(h5py, None if null else obj[()], null)


def _decode(h5py, group):
    kind = _text(group.attrs.get("__type__", "dict"))
    if kind == "list":
        return [
            _member(h5py, group, str(i)) for i in range(int(group.attrs["__len__"]))
        ]
    if "__keys__" in group.attrs:
        return {
            _text(k): _member(h5py, group, _escape(_text(k)))
            for k in group.attrs["__keys__"]
        }
    names = sorted(set(group.keys()) | {n for n in group.attrs if n not in _RESERVED})
    return {_unescape(n): _member(h5py, group, n) for n in names}


def read_experiments(path: str) -> dict:
    """The experiment list as JSON, numbered as the file numbers it -- each
    sweep's crystal-less experiment kept -- so that a reflection's id names its
    own experiment: an .expt as it is; a .rflx decoded from its own encoding,
    read directly. (Read through mxi_convert, the list had been DIALS's,
    lattices only and renumbered from 0, so that an indexed file's unindexed
    spots were the only ones shown -- Graeme's 3798 of 164173, 23 October --
    and with --keep-sweeps mxi_convert writes the whole table besides.)"""
    path = os.path.abspath(path)
    with open(path, "rb") as f:
        hdf5 = f.read(8) == _HDF5_SIGNATURE
    if not hdf5:
        with open(path) as f:
            return json.load(f)
    h5py = h5py_with_filters()
    with h5py.File(path, "r") as f:
        if "experiments" not in f:
            raise SystemExit(f"{path} holds no experiment list")
        return _decode(h5py, f["experiments"])


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


# --- HDF5 --------------------------------------------------------------------


def h5py_with_filters():
    """h5py, with hdf5plugin's filters registered first: mxi writes every
    numeric column of 256 elements or more bitshuffle-LZ4 (HDF5 filter 32008),
    and DECTRIS's frames are compressed the same, so neither reads without it
    -- HDF5 then looks for a plugin directory and fails with that instead. Its
    absence said plainly. (Reading reflections had not imported it: a fresh
    process failed on Graeme's indexed.rflx, 23 October.)"""
    import h5py

    try:
        import hdf5plugin  # noqa: F401 -- bitshuffle and LZ4, registered
    except ImportError:
        raise SystemExit(
            "mxi's tools need hdf5plugin to read compressed files: "
            "pip install hdf5plugin"
        )
    return h5py


# --- reflections -------------------------------------------------------------


def read_reflections(path: str, columns=None) -> dict:
    """A .rflx's reflection tables, joined in the order mxi writes them
    (shortest name first, then by name), as {column: array}: the columns
    asked for, or all; an empty dict where the file holds none."""
    h5py = h5py_with_filters()

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


class Markers:
    """A sweep's spots (or predictions) as markers on the images they are
    recorded on -- not only their centroid's: a spot spans several images, a
    box's median 13 on a small molecule -- for one image at a time: the
    extents worked out once, vectorised, and the markers for the image shown
    picked out as the slider moves. Built for the whole sweep at once they had
    been 6 million points for an integrated 16M sweep's spots and 26 million
    for its predictions, some 15 s to start (Graeme, 26 October); for one image
    they are a few hundred, a few milliseconds.

    A spot's extent: the images its box spans, z0 to z1 - 1 -- the spot
    finder's box, its pixels' images, exactly, in a strong or indexed file;
    integration's, the reflection's extent by the mosaicity model, in an
    integrated one -- else its own image. Not from xyzobs.px.variance, which
    is the centroid's uncertainty, not the spot's spread: so used, with at
    least an image either side, a spot recorded on one image was marked on its
    neighbours too, and against the spot finder's boxes on a 16M sweep only 15
    per cent of spots were marked on just their images -- 119000 images marked
    with no spot, 211000 of a spot's unmarked (Graeme, 26 October).
    Positions are napari's: DIALS's less half a pixel."""

    def __init__(self, table: dict, column: str, experiments):
        if column not in table or "id" not in table:
            xyz = np.zeros((0, 3))
            index = np.zeros(0, np.int64)
        else:
            xyz = np.asarray(table[column], float)
            rows = np.isin(np.asarray(table["id"]), list(experiments)) & np.all(
                np.isfinite(xyz), axis=1
            )
            rows &= np.any(xyz != 0, axis=1)  # a position of zero: none
            index = np.where(rows)[0]
        z = xyz[index, 2]
        own = np.floor(z)
        if len(index) and "bbox" in table:
            box = np.asarray(table["bbox"])[index]
            # Its own image always among them, whatever the box says.
            self.lo = np.minimum(box[:, 4].astype(float), own)
            self.hi = np.maximum(box[:, 5].astype(float) - 1, own)
        else:
            self.lo = self.hi = own
        self.centre = np.floor(z)
        self.yx = np.column_stack([xyz[index, 1], xyz[index, 0]]) - 0.5

    def __len__(self):
        return len(self.yx)

    def on(self, frame: int, window: int = 1):
        """The markers on image `frame` (napari's slice: image frame + 1) --
        or, with `window` images summed for each shown, on any of frame to
        frame + window - 1: their (y, x), and whether each is centred there."""
        last = frame + max(1, int(window)) - 1
        shown = (self.lo <= last) & (self.hi >= frame)
        centre = self.centre[shown]
        return self.yx[shown], (centre >= frame) & (centre <= last)


# --- images ------------------------------------------------------------------


class Images:
    """An NXmx series' frames, read as they are asked for: /entry/data, one
    dataset or several (data_000001, ...), the detector's pixel mask applied
    and overloaded markers zeroed. Looks like an array to napari -- shape,
    dtype, ndim and slicing -- and sums `window` frames from each one asked
    for, so that a shadow shows against the background."""

    def __init__(self, master: str, window: int = 1):
        h5py = h5py_with_filters()
        self.file = h5py.File(master, "r")
        data = self.file["/entry/data"]
        # Every frame in "data" where there is one (a virtual dataset over the
        # blocks, or the frames themselves); the blocks data_000001 ... only
        # where there is not -- both counted had shown every image twice, a
        # slider to 7199 for 3600 (23 October).
        if "data" in data:
            names = ["data"]
        else:
            names = sorted(n for n in data if n.startswith("data_"))
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
