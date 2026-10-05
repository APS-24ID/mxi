"""Each reflection's background against resolution, for looking at before
modelling it.

The background under the reflections comes from air and from what holds the
sample -- water, nylon -- scattering smoothly with angle, so it should be a
smooth function of resolution, and a reflection whose background does not
follow it is suspect: in the backstop shadow, low; in its flare, high. Before
such a model is built, this shows what the data say: the background against
1/d^2 with its running median and spread, and each reflection's background over
that median -- its ratio -- against the azimuth around the beam (where
polarisation, the flare and the backstop would show), against the image (where
a drift through the scan would), and on the detector.

From an integrated table's background.mean, which every mxi or DIALS
integrated table has, and the experiment list for the beam centre. Only numpy:
the figure is mxeq.plots.background's.

    mxeq background integrated.expt integrated.refl -o background.png
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import expt as expt_module
from . import refl


@dataclass
class Backgrounds:
    d: np.ndarray
    background: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    azimuth: np.ndarray  # degrees, 0 along +fast, 90 along +slow
    radius: np.ndarray  # pixels from the beam centre
    centre: tuple[float, float]


def beam_centre(detector, beam) -> tuple[float, float]:
    """Where the beam meets panel 0, in pixels, fast and slow."""
    o = detector.origins[0]
    f = detector.fast_axes[0] / np.linalg.norm(detector.fast_axes[0])
    s = detector.slow_axes[0] / np.linalg.norm(detector.slow_axes[0])
    normal = np.cross(f, s)
    s0 = -np.asarray(beam.direction, float)
    s0 = s0 / np.linalg.norm(s0)
    t = np.dot(o, normal) / np.dot(s0, normal)
    hit = t * s0 - o
    return float(np.dot(hit, f) / detector.pixel_size[0]), float(
        np.dot(hit, s) / detector.pixel_size[1]
    )


def load(expt_path: str, refl_path: str) -> Backgrounds:
    e = expt_module.load(expt_path)
    if not e.experiments or e[0].detector is None or e[0].beam is None:
        raise ValueError(f"{expt_path} has no detector and beam")
    t = refl.load(refl_path)
    c = t.columns
    for name in ("background.mean", "d", "xyzcal.px"):
        if name not in c:
            raise ValueError(f"{refl_path} has no {name}: is it integrated?")
    bg = np.asarray(c["background.mean"], float).ravel()
    d = np.asarray(c["d"], float).ravel()
    px = np.asarray(c["xyzcal.px"], float).reshape(-1, 3)
    keep = np.isfinite(bg) & np.isfinite(d) & (d > 0)
    if "num_pixels.background" in c:
        keep &= (
            np.asarray(c["num_pixels.background"]).ravel() > 0
        )  # a box with no background
    cx, cy = beam_centre(e[0].detector, e[0].beam)
    x, y, z = px[keep, 0], px[keep, 1], px[keep, 2]
    return Backgrounds(
        d[keep],
        bg[keep],
        x,
        y,
        z,
        np.degrees(np.arctan2(y - cy, x - cx)),
        np.hypot(x - cx, y - cy),
        (cx, cy),
    )


@dataclass
class Curve:
    """The running median of the background against 1/d^2, with its spread."""

    s2: np.ndarray  # bin centres, 1/d^2
    median: np.ndarray
    low: np.ndarray  # 5th percentile
    high: np.ndarray  # 95th

    def at(self, s2: np.ndarray) -> np.ndarray:
        return np.interp(s2, self.s2, self.median)


def curve(b: Backgrounds, bins: int = 60) -> Curve:
    s2 = 1.0 / b.d**2
    edges = np.unique(np.quantile(s2, np.linspace(0, 1, bins + 1)))
    which = np.clip(np.searchsorted(edges, s2, side="right") - 1, 0, len(edges) - 2)
    centres, med, lo, hi = [], [], [], []
    for k in range(len(edges) - 1):
        sel = which == k
        if sel.sum() < 5:
            continue
        centres.append(np.median(s2[sel]))
        med.append(np.median(b.background[sel]))
        lo.append(np.percentile(b.background[sel], 5))
        hi.append(np.percentile(b.background[sel], 95))
    return Curve(np.array(centres), np.array(med), np.array(lo), np.array(hi))


def ratios(b: Backgrounds, c: Curve) -> np.ndarray:
    """Each reflection's background over the median at its resolution."""
    expected = c.at(1.0 / b.d**2)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(expected > 0, b.background / expected, np.nan)


def _subset(b: Backgrounds, r: np.ndarray, d_range):
    """The reflections in a resolution range, and their ratios -- measured
    against the curve from all of them, so a range does not move its own
    reference."""
    if d_range is None:
        return b, r, "all resolutions"
    hi, lo = max(d_range), min(d_range)
    sel = (b.d <= hi) & (b.d >= lo)
    part = Backgrounds(
        b.d[sel],
        b.background[sel],
        b.x[sel],
        b.y[sel],
        b.z[sel],
        b.azimuth[sel],
        b.radius[sel],
        b.centre,
    )
    return part, r[sel], f"{hi:g} to {lo:g} A"


def table(
    b: Backgrounds,
    shells: int = 20,
    d_range=None,
    azimuth_bins: int = 24,
    image_bins: int = 12,
) -> str:
    """The background by resolution shell -- its median, its spread about it,
    and how many reflections lie far below or above -- and, over a resolution
    range if one is given, the same against the azimuth and against the image."""
    c = curve(b)
    r_all = ratios(b, c)
    part, r, label = _subset(b, r_all, d_range)
    s2 = 1.0 / part.d**2
    edges = np.quantile(s2, np.linspace(0, 1, shells + 1))
    lines = [
        f"{len(b.d)} reflections with a background; beam centre at pixel {b.centre[0]:.1f}, {b.centre[1]:.1f}",
        f"{len(part.d)} of them in {label}; each ratio against the median at its resolution over all",
        "",
        "by resolution",
        f"  {'d (A)':>15} {'n':>8} {'median bg':>10} {'ratio':>7} {'spread':>7} {'below 0.5':>10} {'above 2':>8}",
    ]
    for k in range(shells):
        sel = (s2 >= edges[k]) & (
            (s2 < edges[k + 1]) if k < shells - 1 else (s2 <= edges[k + 1])
        )
        if not sel.any():
            continue
        dd = part.d[sel]
        med = np.median(part.background[sel])
        rr = r[sel]
        rmed = float(np.nanmedian(rr))
        mad = float(np.nanmedian(np.abs(rr - rmed)))
        lines.append(
            f"  {dd.max():6.2f} -{dd.min():7.2f} {sel.sum():8d} {med:10.3f} {rmed:7.3f} "
            f"{mad / rmed if rmed > 0 else np.nan:7.3f} {np.mean(rr < 0.5) * 100:9.2f}% "
            f"{np.mean(rr > 2.0) * 100:7.2f}%"
        )
    for title, values, bins, unit in (
        ("by azimuth around the beam", part.azimuth, azimuth_bins, "deg"),
        ("by image", part.z, image_bins, ""),
    ):
        lines += [
            "",
            title + f", {label}",
            f"  {('from - to ' + unit).strip():>15} {'n':>8} {'ratio':>7} {'spread':>7} {'below 0.5':>10} "
            f"{'above 2':>8}",
        ]
        lo, hi = (
            (-180.0, 180.0)
            if unit == "deg"
            else (float(values.min()), float(values.max()))
        )
        e = np.linspace(lo, hi, bins + 1)
        for k in range(bins):
            last = k == bins - 1
            sel = (values >= e[k]) & (
                (values <= e[k + 1]) if last else (values < e[k + 1])
            )
            if not sel.any():
                continue
            rr = r[sel]
            rmed = float(np.nanmedian(rr))
            mad = float(np.nanmedian(np.abs(rr - rmed)))
            lines.append(
                f"  {e[k]:7.1f} -{e[k + 1]:7.1f} {sel.sum():8d} {rmed:7.3f} "
                f"{mad / rmed if rmed > 0 else np.nan:7.3f} {np.mean(rr < 0.5) * 100:9.2f}% "
                f"{np.mean(rr > 2.0) * 100:7.2f}%"
            )
    lines += [
        "",
        "  ratio is the median of each reflection's background over the median at its",
        "  resolution, over all reflections; spread its median absolute deviation over it;",
        "  below 0.5 and above 2 the reflections with less than half, or more than twice, that",
        "  median: the backstop shadow, and its flare, would be among them.",
    ]
    return "\n".join(lines)
