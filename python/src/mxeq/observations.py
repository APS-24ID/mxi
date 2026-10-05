"""Every observation of given reflections, side by side.

For chasing a reflection two programs disagree on: every observation of it --
all its symmetry equivalents, grouped by the experiment list's own space group,
so a scaled table's (after mxi_symmetry or dials.symmetry) gives them all where
an integrated table's gives only what its indexing setting relates -- with
whatever the table records of each: where it was predicted, the summed and
fitted intensities, the background, the profile's correlation and mxi's reason
for not fitting it, LP and partiality, the flags that bear on it, and, for a
scaled table, the scale factor and the scaled intensity. A column a table does
not have is left out, so DIALS's tables and mxi's both list.

    mxeq observations scaled.expt scaled.refl 1,1,1 2,2,2 3,1,1
"""

from __future__ import annotations

import numpy as np

from . import equivalents as eq
from . import expt as expt_module
from . import refl

FLAGS = (
    (1 << 8, "sum"),
    (1 << 9, "prf"),
    (1 << 10, "overloaded"),
    (1 << 14, "fg-masked"),
    (1 << 15, "bg-masked"),
    (eq.SCALED, "scaled"),
    (eq.OUTLIER_IN_SCALING, "outlier"),
    (eq.EXCLUDED_FOR_SCALING, "excluded"),
)


def parse_index(text: str) -> list[int]:
    parts = text.replace(" ", "").split(",")
    if len(parts) != 3:
        raise ValueError(f"{text!r} is not h,k,l")
    return [int(p) for p in parts]


def listing(
    expt_path: str, refl_path: str, indices: list[list[int]], limit: int = 0
) -> str:
    e = expt_module.load(expt_path)
    if not e.experiments or e[0].crystal is None or not e[0].crystal.hall:
        raise ValueError(f"{expt_path} has no crystal with a space group")
    hall = e[0].crystal.hall
    t = refl.load(refl_path)
    c = t.columns
    hkl = np.asarray(c["miller_index"]).reshape(-1, 3)
    keys = eq.asu_keys(hkl, hall)
    px = np.asarray(c["xyzcal.px"]).reshape(-1, 3) if "xyzcal.px" in c else None
    flags = (
        np.asarray(c["flags"]).astype(np.int64)
        if "flags" in c
        else np.zeros(t.nrows, np.int64)
    )
    isf = (
        np.asarray(c["inverse_scale_factor"], float)
        if "inverse_scale_factor" in c
        else None
    )

    def col(name):
        return np.asarray(c[name]).ravel() if name in c else None

    columns = [
        ("sum", col("intensity.sum.value"), "{:10.1f}"),
        ("prf", col("intensity.prf.value"), "{:10.1f}"),
        ("bg", col("background.mean"), "{:7.2f}"),
        ("disp", col("background.dispersion"), "{:7.2f}"),
        ("disp t", col("background.dispersion_trimmed"), "{:7.2f}"),
        ("bg z", col("background.z"), "{:7.1f}"),
        ("bg exp", col("background.expected"), "{:7.2f}"),
        ("cc", col("profile.correlation"), "{:5.2f}"),
        ("why", col("profile.failure"), "{:3d}"),
        ("lp", col("lp"), "{:7.4f}"),
        ("part", col("partiality"), "{:4.2f}"),
    ]
    if isf is not None and "intensity.scale.value" in c:
        scaled = np.asarray(c["intensity.scale.value"], float) / np.where(
            isf > 0, isf, np.nan
        )
        columns += [("isf", isf, "{:6.3f}"), ("scaled I", scaled, "{:10.2f}")]
    columns = [(name, v, fmt) for name, v, fmt in columns if v is not None]

    lines = [f"{refl_path}, grouped by {expt_path}'s space group"]
    for target in indices:
        k = eq.asu_keys(np.array([target]), hall)[0]
        rows = np.nonzero(keys == k)[0]
        lines.append(f"== {target}: {len(rows)} observations")
        header = (
            f"  {'hkl':>13} {'x':>7} {'y':>7} {'z':>7}"
            if px is not None
            else f"  {'hkl':>13}"
        )
        for name, _, fmt in columns:
            width = len(fmt.format(0 if "d}" in fmt else 0.0))
            header += f" {name:>{width}}"
        lines.append(header + "  flags")
        shown = rows if not limit else rows[:limit]
        for r in shown:
            text = f"  {str(hkl[r].tolist()):>13}"
            if px is not None:
                text += f" {px[r, 0]:7.1f} {px[r, 1]:7.1f} {px[r, 2]:7.1f}"
            for _, v, fmt in columns:
                value = v[r]
                text += " " + (
                    fmt.format(int(value)) if "d}" in fmt else fmt.format(float(value))
                )
            names = [name for bit, name in FLAGS if flags[r] & bit]
            lines.append(text + "  " + ",".join(names))
        if limit and len(rows) > limit:
            lines.append(f"  ... and {len(rows) - limit} more (--limit 0 for all)")
    return "\n".join(lines)
