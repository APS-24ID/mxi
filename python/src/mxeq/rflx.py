""".rflx: one HDF5 file holding an experiment list, reflection tables or both,
in dxtbx-h5's layout (github.com/graeme-winter/dxtbx-h5, docs/HANDOVER.md;
mxi's docs/rflx.md). Read here from the specification with h5py alone -- not
through dxtbx-h5's own reader, so that agreeing with mxi's C++ is a check on
both.

``experiments_dict`` gives the experiment list's JSON dictionary from an .expt
or a .rflx alike; ``load_table`` the reflections, for ``refl.load``.
"""

from __future__ import annotations

import json

import numpy as np

HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"
RESERVED = {"__type__", "__keys__", "__len__"}


def _h5py():
    """h5py, with hdf5plugin loaded: .rflx numbers are stored bitshuffle-LZ4,
    HDF5 filter 32008, which h5py reads only with the plugin -- without it a
    read fails deep in HDF5 with a missing plugin directory."""
    try:
        import h5py
        import hdf5plugin  # noqa: F401
    except ImportError as e:  # pragma: no cover - a plain environment
        raise RuntimeError(
            "reading a .rflx needs h5py and hdf5plugin (bitshuffle-LZ4, HDF5 "
            f"filter 32008): pip install h5py hdf5plugin ({e})"
        ) from e
    return h5py


def is_hdf5(path: str) -> bool:
    with open(path, "rb") as f:
        return f.read(8) == HDF5_SIGNATURE


def _escape(key: str) -> str:
    """HANDOVER 2.2: % first, then /, then a leading underscore; the empty key."""
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


def _text(value) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def _leaf(value, null: bool):
    """A value's type from its dtype and dataspace alone (HANDOVER 2.2)."""
    h5py = _h5py()

    if null or isinstance(value, h5py.Empty):
        return None
    if isinstance(value, np.ndarray):
        if value.dtype.kind in "OSU":
            return (
                np.vectorize(_text, otypes=[object])(value).tolist()
                if value.size
                else []
            )
        return value.tolist()
    if isinstance(value, (bytes, str)):
        return _text(value)
    if isinstance(value, np.generic):
        return value.item()
    return value


def _member(group, name: str):
    h5py = _h5py()

    if name in group.attrs:
        a = group.attrs.get_id(name)
        return _leaf(
            group.attrs[name], a.get_space().get_simple_extent_type() == h5py.h5s.NULL
        )
    obj = group[name]
    if isinstance(obj, h5py.Group):
        return _decode(obj)
    null = obj.id.get_space().get_simple_extent_type() == h5py.h5s.NULL
    return _leaf(None if null else obj[()], null)


def _decode(group):
    kind = _text(group.attrs.get("__type__", "dict"))
    if kind == "list":
        return [_member(group, str(i)) for i in range(int(group.attrs["__len__"]))]
    if "__keys__" in group.attrs:
        return {
            _text(k): _member(group, _escape(_text(k))) for k in group.attrs["__keys__"]
        }
    names = sorted(set(group.keys()) | {n for n in group.attrs if n not in RESERVED})
    return {_unescape(n): _member(group, n) for n in names}


def experiments_dict(path: str) -> dict:
    """The experiment list's JSON dictionary, from an .expt or a .rflx."""
    if not is_hdf5(path):
        with open(path) as f:
            return json.load(f)
    h5py = _h5py()

    with h5py.File(path, "r") as f:
        if "experiments" not in f:
            raise ValueError(f"{path} holds no experiment list")
        return _decode(f["experiments"])


def summary_dict(path: str) -> dict | None:
    """The processing's summary a .rflx carries -- {"steps": [...]}, each
    program's figures and tables in the order they ran (docs/output.md) -- or
    None if it carries none."""
    h5py = _h5py()

    with h5py.File(path, "r") as f:
        if "summary" not in f:
            return None
        return _decode(f["summary"])


_TAGS = {
    "bool": "bool",
    "int": "int",
    "size_t": "std::size_t",
    "double": "double",
    "vec2_double": "vec2<double>",
    "vec3_double": "vec3<double>",
    "mat3_double": "mat3<double>",
    "miller_index": "cctbx::miller::index<>",
    "int6": "int6",
    "tiny_int_2": "tiny<int,2>",
}


def _column_type(dataset) -> str | None:
    """By dials_type where given, else by dtype and shape (HANDOVER 3.3);
    None for a type mxi's tables have no place for."""
    tag = dataset.attrs.get("dials_type")
    if tag is not None:
        return _TAGS.get(_text(tag))
    shape, dtype = dataset.shape, dataset.dtype
    width = int(np.prod(shape[1:])) if len(shape) > 1 else 1
    if dtype == np.bool_:
        return "bool"
    if dtype == np.float64:
        return {
            1: "double",
            2: "vec2<double>",
            3: "vec3<double>",
            9: "mat3<double>",
        }.get(width)
    if dtype.kind == "i" and width == 3 and len(shape) == 2:
        return "cctbx::miller::index<>"
    if dtype.kind == "i" and width == 6:
        return "int6"
    if dtype == np.uint64 and width == 1:
        return "std::size_t"
    if dtype == np.int32 and width == 1:
        return "int"
    return None


def load_table(path: str):
    """The tables under /reflections, joined in order; shoeboxes left aside,
    which mxeq has no use for."""
    h5py = _h5py()

    from mxeq.refl import ReflectionTable, ReflFormatError

    columns: dict[str, list] = {}
    types: dict[str, str] = {}
    identifiers: dict[int, str] = {}
    nrows = 0
    with h5py.File(path, "r") as f:
        if "reflections" not in f:
            raise ReflFormatError(f"{path} holds no reflections")
        group = f["reflections"]
        names = sorted(
            (n for n in group if isinstance(group[n], h5py.Group)),
            key=lambda n: (len(n), n),
        )
        for name in names:
            table = group[name]
            for i, identifier in zip(
                table.attrs.get("experiment_ids", []),
                table.attrs.get("identifiers", []),
            ):
                identifiers[int(i)] = _text(identifier)
            rows = int(table.attrs["nrows"]) if "nrows" in table.attrs else None
            for key, dataset in table.items():
                if not isinstance(dataset, h5py.Dataset):
                    continue
                kind = _column_type(dataset)
                if kind is None:
                    continue
                array = dataset[()]
                if array.ndim == 3:
                    array = array.reshape(array.shape[0], -1)
                rows = len(array) if rows is None else rows
                types[key] = kind
                columns.setdefault(key, []).append(array)
            nrows += rows or 0
    out = ReflectionTable(nrows=nrows, identifiers=identifiers)
    for key, parts in columns.items():
        out.columns[key] = np.concatenate(parts)
        out.types[key] = types[key]
    out.validate()
    return out


def load_shoeboxes(path: str, rows) -> dict:
    """The shoeboxes of the rows asked for, numbered as load_table numbers
    them (its tables joined in the same order): each box's bbox (x0, x1, y0,
    y1, z0, z1) and its data, background and mask, z then y then x, as mxi
    writes them (docs/rflx.md). Returns {"bbox": (n, 6), "data": [...],
    "background": [...], "mask": [...]} in the order of `rows`; a row with no
    box gets None."""
    import numpy as np

    h5py = _h5py()
    rows = [int(r) for r in rows]
    out = {"bbox": [], "data": [], "background": [], "mask": []}
    wanted: dict[int, list[int]] = {}
    for at, r in enumerate(rows):
        wanted.setdefault(r, []).append(at)
    found = [None] * len(rows)
    with h5py.File(path, "r") as f:
        group = f["reflections"]
        names = sorted(
            (n for n in group if isinstance(group[n], h5py.Group)),
            key=lambda n: (len(n), n),
        )
        first = 0
        for name in names:
            table = group[name]
            count = None
            for k in table:
                if isinstance(table[k], h5py.Dataset):
                    count = table[k].shape[0]
                    break
            if "shoebox" in table:
                g = table["shoebox"]
                bbox = g["bbox"][()]
                count = len(bbox) if count is None else count
                sizes = (
                    (bbox[:, 1] - bbox[:, 0])
                    * (bbox[:, 3] - bbox[:, 2])
                    * (bbox[:, 5] - bbox[:, 4])
                ).astype(np.int64)
                start = np.concatenate([[0], np.cumsum(sizes)])
                for r in [r for r in wanted if first <= r < first + len(bbox)]:
                    i = r - first
                    a, b = int(start[i]), int(start[i + 1])
                    x0, x1, y0, y1, z0, z1 = (int(v) for v in bbox[i])
                    shape = (z1 - z0, y1 - y0, x1 - x0)
                    box = (
                        bbox[i].copy(),
                        g["shoebox_data"][a:b].reshape(shape),
                        g["shoebox_background"][a:b].reshape(shape),
                        g["shoebox_mask"][a:b].reshape(shape),
                    )
                    for at in wanted[r]:
                        found[at] = box
            first += count or 0
    for box in found:
        if box is None:
            for k in out:
                out[k].append(None)
        else:
            for k, v in zip(("bbox", "data", "background", "mask"), box):
                out[k].append(v)
    return out
