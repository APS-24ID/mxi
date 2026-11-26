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


# The HTML page: binned data only, so that it stays small whatever the
# number of spots; drawn with plotly from its CDN, as mxeq residuals is.

NAMES = ["dh", "dk", "dl", "radial", "along the rotation", "across"]
UNITS = ["index units", "index units", "index units", "x 1e-4", "deg", "deg"]
SCALES = [1.0, 1.0, 1.0, 1e4, 1.0, 1.0]
MAPS = [(4, 5), (0, 1), (1, 2), (0, 2)]  # along x across, dh x dk, dk x dl, dh x dl


def _binned(x, y, n_bins, equal_count=True):
    """Median and robust width of y in bins of x."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 2 * n_bins:
        return None
    edges = (
        np.quantile(x, np.linspace(0, 1, n_bins + 1))
        if equal_count
        else np.linspace(x.min(), x.max(), n_bins + 1)
    )
    centre, median, width = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x <= hi)
        if m.sum() < 5:
            continue
        centre.append(float(np.median(x[m])))
        median.append(float(np.median(y[m])))
        width.append(_width(y[m]))
    return {"x": centre, "median": median, "width": width}


def build(o: Offsets, bins: int = 61, map_bins: int = 61, n_trend: int = 20) -> dict:
    values = _values(o)
    stats = [summary(v) for v in values]
    out = {
        "indexed": o.indexed,
        "total": o.total,
        "names": NAMES,
        "units": UNITS,
        "summary": [
            {
                "name": name,
                "unit": unit,
                "median": s["median"] * scale,
                "width": s["width"] * scale,
                "rms": s["rms"] * scale,
                "beyond_3": s["beyond_3"],
                "n": s["n"],
            }
            for name, unit, scale, s in zip(NAMES, UNITS, SCALES, stats)
        ],
    }
    # Histograms across +-8 widths about nought, counts per bin; the page
    # draws the Gaussian of the robust width over each.
    hists = []
    for v, s, scale in zip(values, stats, SCALES):
        w = s["width"] if s["width"] and np.isfinite(s["width"]) else 1.0
        edges = np.linspace(-8 * w, 8 * w, bins + 1)
        counts, _ = np.histogram(v[np.isfinite(v)], bins=edges)
        hists.append(
            {
                "centre": (0.5 * (edges[:-1] + edges[1:]) * scale).tolist(),
                "counts": counts.tolist(),
                "step": float((edges[1] - edges[0]) * scale),
                "below": int(np.sum(v < edges[0])),
                "above": int(np.sum(v > edges[-1])),
            }
        )
    out["histograms"] = hists
    # The reciprocal-space maps: two offsets at once, +-6 widths each.
    maps = []
    for i, j in MAPS:
        a, b = values[i], values[j]
        ok = np.isfinite(a) & np.isfinite(b)
        wa, wb = stats[i]["width"], stats[j]["width"]
        ea = np.linspace(-6 * wa, 6 * wa, map_bins + 1)
        eb = np.linspace(-6 * wb, 6 * wb, map_bins + 1)
        h, _, _ = np.histogram2d(a[ok], b[ok], bins=[ea, eb])
        maps.append(
            {
                "x": (0.5 * (ea[:-1] + ea[1:]) * SCALES[i]).tolist(),
                "y": (0.5 * (eb[:-1] + eb[1:]) * SCALES[j]).tolist(),
                "z": h.T.astype(int).tolist(),  # rows are y
                "xname": NAMES[i],
                "yname": NAMES[j],
                "xunit": UNITS[i],
                "yunit": UNITS[j],
            }
        )
    out["maps"] = maps
    # Against resolution (1/d^2, equal counts) and against image.
    s2 = 1.0 / o.d**2
    out["against_resolution"] = []
    out["against_image"] = []
    for k in (3, 4, 5):
        r = _binned(s2, values[k], n_trend)
        if r:
            r["x"] = [1.0 / np.sqrt(v) for v in r["x"]]
            r["median"] = [m * SCALES[k] for m in r["median"]]
            r["width"] = [w * SCALES[k] for w in r["width"]]
        out["against_resolution"].append(r)
        r = _binned(o.image, values[k], n_trend, equal_count=False)
        if r:
            r["median"] = [m * SCALES[k] for m in r["median"]]
            r["width"] = [w * SCALES[k] for w in r["width"]]
        out["against_image"].append(r)
    return out


PAGE = """<!doctype html>
<meta charset="utf-8">
<title>__TITLE__</title>
<script src="__PLOTLY__"></script>
<style>
 body { font: 14px/1.5 system-ui, sans-serif; margin: 2rem auto; max-width: 1100px; color: #222; }
 h1 { font-size: 1.4rem; } h2 { font-size: 1.1rem; margin-top: 2.2rem; border-top: 1px solid #ddd; padding-top: 1rem; }
 table { border-collapse: collapse; margin: 0.5rem 0 1rem; }
 td, th { padding: 0.2rem 0.8rem; text-align: right; border-bottom: 1px solid #eee; }
 th { font-weight: 600; } .note { color: #666; font-size: 0.9rem; }
 .grid { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 0.4rem; }
 .grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 0.4rem; }
</style>
<h1>__TITLE__</h1>
<p class="note">__SUBTITLE__</p>
<div id="body"></div>
<script>
const R = __DATA__;
const body = document.getElementById("body");
const cfg = {displayModeBar: false, responsive: true};
const margin = {t: 36, l: 55, r: 10, b: 45};
function el(tag, cls, text) { const e = document.createElement(tag); if (cls) e.className = cls; if (text) e.textContent = text; return e; }
function plot(where, traces, layout) { const d = el("div"); where.appendChild(d); Plotly.newPlot(d, traces, layout, cfg); }
function fmt(v) { return Math.abs(v) >= 0.1 ? v.toFixed(3) : v.toPrecision(3); }

body.appendChild(el("h2", null, "Summary"));
const t = el("table");
t.innerHTML = "<tr><th></th><th>median</th><th>width</th><th>RMS</th><th>beyond 3 widths</th></tr>" +
  R.summary.map(s => `<tr><td style="text-align:left">${s.name} (${s.unit})</td>` +
    `<td>${fmt(s.median)}</td><td>${fmt(s.width)}</td><td>${fmt(s.rms)}</td>` +
    `<td>${(100 * s.beyond_3).toFixed(2)}%</td></tr>`).join("");
body.appendChild(t);
body.appendChild(el("p", "note",
  "Width is 1.4826 x the median absolute deviation, the sigma a Gaussian core would have. A Gaussian has 0.27% " +
  "beyond three widths; more is a tail -- chance coincidences, a second lattice, outliers. dh, dk, dl are fractional " +
  "index offsets along the crystal's axes; radial, along the rotation and across are each reflection's own frame, " +
  "relative to |r|: radial for the cell, wavelength and distance, along for the rocking width and image centroids, " +
  "across for the beam's divergence."));

body.appendChild(el("h2", null, "Histograms"));
body.appendChild(el("p", "note",
  "Counts on a logarithmic scale, so that tails and asymmetry show; the red line is the Gaussian of the robust " +
  "width, scaled to the same number of spots. Spots beyond +-8 widths are counted in the titles."));
const h = el("div", "grid"); body.appendChild(h);
R.histograms.forEach((g, k) => {
  const s = R.summary[k], w = s.width, n = s.n;
  const gauss = g.centre.map(x => n * g.step / (w * Math.sqrt(2 * Math.PI)) * Math.exp(-0.5 * ((x - s.median) / w) ** 2));
  plot(h, [
    {x: g.centre, y: g.counts, type: "bar", marker: {color: "rgba(31,119,180,0.65)"}, name: "spots"},
    {x: g.centre, y: gauss, mode: "lines", line: {color: "#d62728", width: 1.5}, name: "Gaussian"}
  ], {height: 280, title: `${s.name} (${s.unit}); ${g.below} below, ${g.above} above`, showlegend: false,
      bargap: 0, yaxis: {type: "log", title: "spots"}, xaxis: {title: s.unit}, margin});
});

body.appendChild(el("h2", null, "Reciprocal-space maps"));
body.appendChild(el("p", "note",
  "Two offsets at once, counts on a logarithmic colour scale, +-6 widths each way. One round blob is a clean " +
  "lattice; a second blob a second lattice or domain; a streak a spread -- a bent crystal, mosaicity -- in that " +
  "direction; an even floor chance coincidences."));
const m = el("div", "grid2"); body.appendChild(m);
R.maps.forEach(map => {
  const z = map.z.map(row => row.map(v => v > 0 ? Math.log10(v) : null));
  plot(m, [{x: map.x, y: map.y, z, type: "heatmap", colorscale: "Viridis", colorbar: {title: "log10 n"}}],
       {height: 420, title: `${map.yname} against ${map.xname}`,
        xaxis: {title: `${map.xname} (${map.xunit})`}, yaxis: {title: `${map.yname} (${map.yunit})`}, margin});
});

function against(title, blocks, xtitle, note, reversed) {
  body.appendChild(el("h2", null, title));
  body.appendChild(el("p", "note", note));
  const g = el("div", "grid"); body.appendChild(g);
  blocks.forEach((b, i) => {
    if (!b) return;
    const k = 3 + i;
    const upper = b.median.map((v, j) => v + b.width[j]), lower = b.median.map((v, j) => v - b.width[j]);
    plot(g, [
      {x: b.x.concat(b.x.slice().reverse()), y: upper.concat(lower.reverse()), fill: "toself",
       fillcolor: "rgba(31,119,180,0.15)", line: {width: 0}, hoverinfo: "skip", showlegend: false},
      {x: b.x, y: b.median, mode: "lines+markers", line: {color: "#1f77b4"}, marker: {size: 4}, showlegend: false},
      {x: [Math.min(...b.x), Math.max(...b.x)], y: [0, 0], mode: "lines", line: {dash: "dot", color: "#999", width: 1}, showlegend: false}
    ], {height: 280, title: `${R.names[k]} (${R.units[k]})`, margin,
        xaxis: Object.assign({title: xtitle}, reversed ? {autorange: "reversed"} : {})});
  });
}
against("Against resolution", R.against_resolution, "d (A)",
  "Median and robust width in bins of equal count. A radial median that is constant but not nought is a cell " +
  "error; one that trends with resolution is not -- the detector model, the distance, or centroids biased " +
  "by resolution.", true);
against("Against image", R.against_image, "image",
  "Median and width through the scan: a median that wanders is the crystal or the beam moving in a way the " +
  "refinement did not follow.", false);
</script>
"""


def write(
    o: Offsets, path: str, title: str = "lattice map", plotly_src: str | None = None
) -> str:
    import html
    import json

    from .htmlreport import PLOTLY_CDN

    data = build(o)
    subtitle = (
        f"{o.indexed} of {o.total} spots indexed "
        f"({100.0 * o.indexed / max(o.total, 1):.1f}%): each against its lattice point, "
        "the crystal as refined, scan-varying where it is."
    )
    page = (
        PAGE.replace("__TITLE__", html.escape(title))
        .replace("__SUBTITLE__", html.escape(subtitle))
        .replace("__PLOTLY__", html.escape(plotly_src or PLOTLY_CDN, quote=True))
        .replace("__DATA__", json.dumps(data))
    )
    with open(path, "w") as f:
        f.write(page)
    return report(o, histograms=False) + f"\n\nWrote {path}"
