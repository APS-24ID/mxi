"""Where indexed spots lie around their lattice points, in reciprocal space.

`mxeq residuals` measures positions in the images' frame -- pixels and images,
after integration. This measures them where indexing judges them: each
indexed spot's reciprocal-lattice point (the rlp column) against its lattice
point A h, the crystal refined, scan-varying where it is, taken at the spot's
image. Two frames:

* the CRYSTAL'S AXES, as fractional index offsets dh, dk, dl -- whether one
  direction is worse than the others, as a wrong cell or a second lattice
  offset along one axis would make it;
* each REFLECTION'S OWN FRAME, as offsets relative to |r| (angles, in
  degrees, or for the radial one a fraction): RADIAL, along r -- the cell,
  the wavelength or the distance, an error that grows with resolution if the
  cell is off; ALONG THE ROTATION, the direction the rotation moves r -- the
  rocking width and the image centroids; and ACROSS, out of the rotation
  plane -- the beam's divergence. As fractions of |r| they pool across
  resolution.

For each: a robust width (1.4826 median absolute deviation), the RMS, the
fraction beyond three widths, and a histogram. A clean lattice is one narrow
peak; chance coincidences a floor rising to the indexing tolerance; a split
crystal a shoulder or a second peak (Graeme, 26 November: on data where the
first lattice is a good solution -- the map presumes it).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import expt, refl


@dataclass
class Offsets:
    """One indexed spot each."""

    dhkl: np.ndarray  # (n, 3) fractional index offsets
    radial: np.ndarray  # (r_obs - r_pred) . e_r / |r_pred|
    along: np.ndarray  # degrees, along the rotation
    across: np.ndarray  # degrees, out of the rotation plane
    d: np.ndarray  # resolution, A
    image: np.ndarray  # the spot's image, from xyzobs.px
    indexed: int
    total: int


def _settings(crystal: expt.Crystal, scan: expt.Scan | None, z: np.ndarray):
    """The setting matrix at each spot: interpolated between scan points for a
    scan-varying crystal (one per image boundary, the first at the scan's
    start), the static one otherwise."""
    if not crystal.scan_varying or scan is None:
        return np.broadcast_to(crystal.setting, (len(z), 3, 3))
    points = crystal.a_at_scan_points
    t = np.clip(z - (scan.image_range[0] - 1), 0.0, len(points) - 1.0)
    lo = np.floor(t).astype(int)
    hi = np.minimum(lo + 1, len(points) - 1)
    f = (t - lo)[:, None, None]
    return points[lo] * (1.0 - f) + points[hi] * f


def gather(table: refl.ReflectionTable, experiments: expt.ExperimentList) -> Offsets:
    """Every indexed spot's offset from its lattice point."""
    if "rlp" not in table.columns or "miller_index" not in table.columns:
        raise ValueError(
            "needs the rlp and miller_index columns: an indexed or refined file"
        )
    hkl = np.asarray(table["miller_index"]).astype(float)
    r = np.asarray(table["rlp"], dtype=float)
    ids = (
        np.asarray(table["id"]).astype(int)
        if "id" in table.columns
        else np.zeros(len(hkl), dtype=int)
    )
    z = np.asarray(table["xyzobs.px.value"], dtype=float)[:, 2]
    indexed = np.any(hkl != 0, axis=1) & (ids >= 0) & np.all(np.isfinite(r), axis=1)

    parts = []
    for i, e in enumerate(experiments.experiments):
        if e.crystal is None:
            continue
        sel = np.nonzero(indexed & (ids == i))[0]
        if sel.size == 0:
            continue
        A = _settings(e.crystal, e.scan, z[sel])
        pred = np.einsum("nij,nj->ni", A, hkl[sel])
        dhkl = np.einsum("nij,nj->ni", np.linalg.inv(A), r[sel]) - hkl[sel]
        delta = r[sel] - pred
        size = np.linalg.norm(pred, axis=1)
        e_r = pred / size[:, None]
        axis = (
            e.goniometer.axis_in_crystal_frame
            if e.goniometer is not None
            else np.array([1.0, 0.0, 0.0])
        )
        e_along = np.cross(axis, e_r)
        reach = np.linalg.norm(e_along, axis=1)
        # Near the rotation axis the direction of rotation is undefined.
        far = reach > 0.05
        e_along = e_along / np.where(far, reach, 1.0)[:, None]
        e_across = np.cross(e_r, e_along)
        deg = 180.0 / np.pi
        along = np.where(far, np.sum(delta * e_along, axis=1) / size * deg, np.nan)
        across = np.where(far, np.sum(delta * e_across, axis=1) / size * deg, np.nan)
        parts.append(
            (
                dhkl,
                np.sum(delta * e_r, axis=1) / size,
                along,
                across,
                1.0 / size,
                z[sel],
            )
        )
    if not parts:
        raise ValueError("no indexed spots with a crystal")
    cat = [np.concatenate([p[k] for p in parts]) for k in range(6)]
    return Offsets(
        dhkl=cat[0],
        radial=cat[1],
        along=cat[2],
        across=cat[3],
        d=cat[4],
        image=cat[5],
        indexed=int(indexed.sum()),
        total=len(hkl),
    )


def _width(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def summary(x: np.ndarray) -> dict:
    """Median, robust width, RMS and the fraction beyond three widths."""
    x = x[np.isfinite(x)]
    w = _width(x)
    m = float(np.median(x)) if x.size else float("nan")
    return {
        "n": int(x.size),
        "median": m,
        "width": w,
        "rms": float(np.sqrt(np.mean(x * x))) if x.size else float("nan"),
        "beyond_3": float(np.mean(np.abs(x - m) > 3 * w)) if x.size and w > 0 else 0.0,
    }


def histogram(x: np.ndarray, width: float, bins: int = 25, reach: float = 6.0):
    """Counts in bins across +-reach widths, and how many lie beyond."""
    x = x[np.isfinite(x)]
    edges = np.linspace(-reach * width, reach * width, bins + 1)
    counts, _ = np.histogram(x, bins=edges)
    return edges, counts, int(np.sum(x < edges[0])), int(np.sum(x > edges[-1]))


def _bars(edges, counts, below, above, unit: str, scale: float) -> list[str]:
    most = max(int(counts.max()), 1)
    lines = [f"      {'below':>10s} {below:8d}"]
    for lo, hi, c in zip(edges[:-1], edges[1:], counts):
        bar = "#" * int(round(40 * c / most))
        lines.append(f"      {0.5 * (lo + hi) * scale:10.4f} {c:8d} {bar}")
    lines.append(f"      {'above':>10s} {above:8d}   ({unit})")
    return lines


QUANTITIES = (
    ("dh", "index units", 1.0),
    ("dk", "index units", 1.0),
    ("dl", "index units", 1.0),
    ("radial", "relative, x 1e-4", 1e4),
    ("along the rotation", "degrees", 1.0),
    ("across", "degrees", 1.0),
)


def _values(o: Offsets):
    return [o.dhkl[:, 0], o.dhkl[:, 1], o.dhkl[:, 2], o.radial, o.along, o.across]


def report(o: Offsets, shells: int = 8, histograms: bool = True) -> str:
    lines = [
        f"Indexed {o.indexed} of {o.total} spots "
        f"({100.0 * o.indexed / max(o.total, 1):.1f}%); offsets from the refined "
        "lattice points",
        "",
        f"  {'':22s} {'median':>10s} {'width':>10s} {'RMS':>10s} "
        f"{'beyond 3 widths':>16s}",
    ]
    values = _values(o)
    stats = [summary(v) for v in values]
    for (name, unit, scale), s in zip(QUANTITIES, stats):
        lines.append(
            f"  {name:22s} {s['median'] * scale:10.4f} {s['width'] * scale:10.4f} "
            f"{s['rms'] * scale:10.4f} {100.0 * s['beyond_3']:15.2f}%   {unit}"
        )
    lines += [
        "",
        "  (width: 1.4826 x the median absolute deviation; a Gaussian has 0.27%",
        "  beyond three widths -- more is a tail: chance coincidences, a second",
        "  lattice, outliers)",
    ]

    ok = np.isfinite(o.d)
    if ok.sum() >= shells * 20:
        edges = np.quantile(1.0 / o.d[ok] ** 2, np.linspace(0, 1, shells + 1))
        s2 = 1.0 / o.d**2
        lines += [
            "",
            "By resolution, equal numbers of spots: the radial median trending",
            "outward is a cell, wavelength or distance error",
            f"  {'d_max':>7s} {'d_min':>7s} {'n':>7s} {'radial med':>11s} "
            f"{'radial w':>10s} {'along w':>9s} {'across w':>9s}",
        ]
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = ok & (s2 >= lo) & (s2 <= hi)
            lines.append(
                f"  {1 / np.sqrt(lo):7.2f} {1 / np.sqrt(hi):7.2f} {int(m.sum()):7d} "
                f"{np.median(o.radial[m]) * 1e4:11.2f} {_width(o.radial[m]) * 1e4:10.2f} "
                f"{_width(o.along[m]):9.4f} {_width(o.across[m]):9.4f}"
            )
        lines.append("  (radial x 1e-4; along and across in degrees)")

    if histograms:
        for (name, unit, scale), v, s in zip(QUANTITIES, values, stats):
            if not s["width"] or not np.isfinite(s["width"]):
                continue
            lines += ["", f"  {name}: histogram across +-6 widths"]
            lines += _bars(*histogram(v, s["width"]), unit=unit, scale=scale)
    return "\n".join(lines)


def run(path: str, shells: int = 8, histograms: bool = True) -> str:
    return report(
        gather(refl.load(path), expt.load(path)), shells=shells, histograms=histograms
    )
