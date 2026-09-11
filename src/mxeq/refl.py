"""Read and write DIALS reflection tables (``.refl``) without cctbx.

A ``.refl`` file is a msgpack document.  At the top level it is a two-element
array whose first element is the tag ``dials::af::reflection_table`` and whose
second is a map with the row count, the experiment identifiers, and the
columns.  Each column is itself a two-element array of a C++ type name and a
payload; for every numeric type the payload is a single ``bin`` holding the
values packed little-endian, with the components of a compound type adjacent
rather than in separate arrays.

Nothing here knows about cctbx.  That is the point: the checker has to run in a
container with no DIALS build in it, against files written by a DIALS build
somewhere else.

STATUS
------
The structure above is taken from DIALS' own msgpack adapter and round-trips
against :func:`write` in ``tests/``.  It has *not* been validated against a
``.refl`` written by a real DIALS.  ``mxeq inspect`` exists so that the first
contact with a real file says exactly what the difference is: it walks the
document without assuming any of the key names below.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field

import msgpack
import numpy as np

TAG = "dials::af::reflection_table"

#: The tag accepted for the map key holding the experiment identifiers.  DIALS
#: has spelled this more than one way; both are read, the first is written.
IDENTIFIER_KEYS = ("identifiers", "experiment_identifiers")

# C++ type name -> (numpy scalar dtype, components per row).  A column of
# ``vec3<double>`` is 3N doubles, not N of anything, so the shape is imposed
# here and not inferred from the byte count alone -- inferring it would silently
# accept a truncated payload.
NUMERIC_TYPES = {
    "int": ("<i4", 1),
    "std::size_t": ("<u8", 1),
    "double": ("<f8", 1),
    "bool": ("|u1", 1),
    "vec2<double>": ("<f8", 2),
    "vec3<double>": ("<f8", 3),
    "mat3<double>": ("<f8", 9),
    "int6": ("<i4", 6),
    "cctbx::miller::index<>": ("<i4", 3),
    "tiny<int,2>": ("<i4", 2),
}

#: Types carried through without being decoded.  A shoebox is a nested record
#: with its own variable-length pixel arrays, and no equivalence check here
#: looks inside one; decoding it would be work spent to produce something
#: nothing reads.
OPAQUE_TYPES = ("Shoebox<>",)


class ReflFormatError(Exception):
    """The file is msgpack, but not a reflection table this can read."""


@dataclass
class ReflectionTable:
    """Columns as numpy arrays, plus the row count and the identifiers.

    ``nrows`` is read from the file rather than derived from the columns, and
    :meth:`validate` compares the two.  A file whose header and payload
    disagree is a corrupt file, and saying so is more useful than quietly
    believing one of them.
    """

    nrows: int
    columns: dict[str, np.ndarray] = field(default_factory=dict)
    identifiers: dict[int, str] = field(default_factory=dict)
    types: dict[str, str] = field(default_factory=dict)
    opaque: dict[str, str] = field(default_factory=dict)

    def __len__(self) -> int:
        return self.nrows

    def __contains__(self, key: str) -> bool:
        return key in self.columns

    def __getitem__(self, key: str) -> np.ndarray:
        try:
            return self.columns[key]
        except KeyError:
            raise KeyError(
                f"no column {key!r}; have {', '.join(sorted(self.columns))}"
            ) from None

    def validate(self) -> None:
        for name, values in self.columns.items():
            if len(values) != self.nrows:
                raise ReflFormatError(
                    f"column {name!r} has {len(values)} rows, header says {self.nrows}"
                )

    def select(self, mask: np.ndarray) -> ReflectionTable:
        """A new table holding the selected rows.  Opaque columns are dropped."""
        mask = np.asarray(mask)
        if mask.dtype == bool:
            n = int(mask.sum())
        else:
            n = len(mask)
        return ReflectionTable(
            nrows=n,
            columns={k: v[mask] for k, v in self.columns.items()},
            identifiers=dict(self.identifiers),
            types=dict(self.types),
        )


def _decode_column(name: str, type_name: str, payload: object) -> np.ndarray:
    if type_name == "std::string":
        if not isinstance(payload, (list, tuple)):
            raise ReflFormatError(f"column {name!r}: string payload is not an array")
        return np.array(
            [p.decode() if isinstance(p, bytes) else p for p in payload], dtype=object
        )

    dtype, width = NUMERIC_TYPES[type_name]
    if not isinstance(payload, (bytes, bytearray)):
        raise ReflFormatError(
            f"column {name!r} of type {type_name}: payload is "
            f"{type(payload).__name__}, expected a binary blob"
        )
    values = np.frombuffer(bytes(payload), dtype=dtype)
    if width == 1:
        out = values
    else:
        if len(values) % width:
            raise ReflFormatError(
                f"column {name!r} of type {type_name}: {len(values)} scalars "
                f"is not a multiple of {width}"
            )
        out = values.reshape(-1, width)
    if type_name == "bool":
        out = out.astype(bool)
    return out


def loads(raw: bytes) -> ReflectionTable:
    """Decode a reflection table from the bytes of a ``.refl`` file."""
    doc = msgpack.unpackb(raw, raw=True, strict_map_key=False)
    return _from_document(doc)


def load(path: str) -> ReflectionTable:
    with open(path, "rb") as f:
        return loads(f.read())


def _text(value: object) -> object:
    return value.decode() if isinstance(value, bytes) else value


def _unraw(obj: object) -> object:
    """Recursively turn msgpack ``raw`` bytes back into ``str`` for map keys."""
    if isinstance(obj, dict):
        return {_text(k): _unraw(v) for k, v in obj.items()}
    return obj


def _from_document(doc: object) -> ReflectionTable:
    if not isinstance(doc, (list, tuple)) or len(doc) != 2:
        raise ReflFormatError(
            "not a reflection table: the document is not a two-element array. "
            "Run `mxeq inspect` on it."
        )
    tag = _text(doc[0])
    if tag != TAG:
        raise ReflFormatError(f"unexpected tag {tag!r}, expected {TAG!r}")

    body = _unraw(doc[1])
    if not isinstance(body, dict):
        raise ReflFormatError("the reflection table body is not a map")

    if "nrows" not in body:
        raise ReflFormatError(
            f"no 'nrows' in the table body; keys are {', '.join(map(str, body))}"
        )
    nrows = int(body["nrows"])

    identifiers: dict[int, str] = {}
    for key in IDENTIFIER_KEYS:
        if key in body and isinstance(body[key], dict):
            identifiers = {int(k): _text(v) for k, v in body[key].items()}
            break

    data = body.get("data")
    if not isinstance(data, dict):
        raise ReflFormatError("no 'data' map in the table body")

    table = ReflectionTable(nrows=nrows, identifiers=identifiers)
    for raw_name, column in data.items():
        name = str(_text(raw_name))
        if not isinstance(column, (list, tuple)) or len(column) != 2:
            raise ReflFormatError(f"column {name!r} is not a [type, payload] pair")
        type_name = str(_text(column[0]))
        if type_name in OPAQUE_TYPES:
            table.opaque[name] = type_name
            continue
        if type_name not in NUMERIC_TYPES and type_name != "std::string":
            print(
                f"mxeq: skipping column {name!r}: unknown type {type_name!r}",
                file=sys.stderr,
            )
            table.opaque[name] = type_name
            continue
        table.columns[name] = _decode_column(name, type_name, column[1])
        table.types[name] = type_name

    table.validate()
    return table


def _encode_column(values: np.ndarray, type_name: str) -> object:
    if type_name == "std::string":
        return [str(v) for v in values]
    dtype, _ = NUMERIC_TYPES[type_name]
    if type_name == "bool":
        return np.asarray(values, dtype=bool).astype("|u1").tobytes()
    return np.ascontiguousarray(values, dtype=dtype).tobytes()


def dumps(table: ReflectionTable) -> bytes:
    """Encode a reflection table.  Used for fixtures, and to round-trip the reader."""
    table.validate()
    data = {}
    for name, values in table.columns.items():
        type_name = table.types.get(name)
        if type_name is None:
            raise ReflFormatError(f"column {name!r} has no recorded C++ type")
        data[name] = [type_name, _encode_column(values, type_name)]
    body = {
        "nrows": table.nrows,
        IDENTIFIER_KEYS[0]: {int(k): v for k, v in table.identifiers.items()},
        "data": data,
    }
    return msgpack.packb([TAG, body], use_bin_type=True)


def write(path: str, table: ReflectionTable) -> None:
    with open(path, "wb") as f:
        f.write(dumps(table))


def structure(raw: bytes, max_items: int = 40) -> list[str]:
    """Describe an unknown msgpack document, assuming nothing about its keys.

    This is what to reach for when :func:`loads` refuses a real file.  It never
    raises on structure, so its output is the diagnostic.
    """
    lines: list[str] = []

    def walk(obj: object, indent: int, label: str) -> None:
        pad = "  " * indent
        if isinstance(obj, dict):
            lines.append(f"{pad}{label}map, {len(obj)} keys")
            for i, (k, v) in enumerate(obj.items()):
                if i >= max_items:
                    lines.append(f"{pad}  ... {len(obj) - max_items} more")
                    break
                walk(v, indent + 1, f"{_text(k)!r}: ")
        elif isinstance(obj, (list, tuple)):
            lines.append(f"{pad}{label}array, {len(obj)} items")
            for i, v in enumerate(obj):
                if i >= max_items:
                    lines.append(f"{pad}  ... {len(obj) - max_items} more")
                    break
                walk(v, indent + 1, f"[{i}] ")
        elif isinstance(obj, (bytes, bytearray)):
            text = _text(obj)
            if isinstance(text, str) and text.isprintable() and len(text) < 60:
                lines.append(f"{pad}{label}{text!r}")
            else:
                lines.append(f"{pad}{label}binary, {len(obj)} bytes")
        else:
            lines.append(f"{pad}{label}{obj!r}")

    try:
        doc = msgpack.unpackb(raw, raw=True, strict_map_key=False)
    except Exception as exc:  # noqa: BLE001 - the exception *is* the report
        return [f"not decodable as msgpack: {exc}"]
    walk(doc, 0, "")
    return lines
