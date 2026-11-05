"""What a step produced, from the summary its file carries (docs/output.md):
its facts, its tables, its plots -- as data, drawn by the window with
matplotlib. No tkinter or matplotlib here, so that it is tested anywhere.

Each step's output file carries the summary of every step before it, its own
last; a step's record is its own program's, the last of them. The plots are
mxi report's, which reads the same records: a figure of the record this does
not know is still shown, as a fact, so that one added to a command later
appears without a change here.
"""

import math
import os
from dataclasses import dataclass, field

import mxi_tools

from mxi_ui.steps import by_key


@dataclass
class Series:
    name: str
    x: list
    y: list
    thin: bool = False  # one of many, drawn grey and fine, in no legend


@dataclass
class Plot:
    """One plot: its series, and how to draw them. On a resolution axis x is
    1/d^2, ticked in d; `right` series, if any, share x on an axis of their
    own (multiplicity beside completeness)."""

    title: str
    xlabel: str
    ylabel: str
    series: list
    resolution: bool = False
    log: bool = False
    bars: bool = False
    right: list = field(default_factory=list)
    right_label: str = ""
    integer: bool = False  # x in whole numbers: cycles, images
    # "line" (bars as above), or "matrix" -- `matrix` drawn as an image,
    # its scale from `low` to 1 -- or "points", the series as points alone.
    kind: str = "line"
    matrix: list = field(default_factory=list)
    low: float = 0.0


@dataclass
class Table:
    title: str
    columns: list
    rows: list


# --- reading --------------------------------------------------------------


def steps_of(path):
    """Every step the file's summary records, in the order they ran. A file
    that cannot be read raises, for the window to say why: caught here and
    taken for a file with no summary, a missing filter had passed for one."""
    return mxi_tools.read_summary(path).get("steps", [])


def record(directory, key):
    """(path, record) for a step: its output file in the directory, and its
    own program's record in that file's summary -- record None where the file
    is missing or records no such step. mxi correlation's is the JSON it
    writes (--json), its program said mxi_correlation_json."""
    step = by_key(key)
    if key == "correlation":
        return correlation_record(os.path.join(directory, step.output))
    if not step.output or not step.output.endswith(".rflx"):
        return None, None
    path = os.path.join(directory, step.output)
    if not os.path.isfile(path):
        return path, None
    mine = [s for s in steps_of(path) if s.get("program") == f"mxi_{key}"]
    return path, (mine[-1] if mine else None)


def correlation_record(path):
    """mxi correlation's --json, as a record: the matrix, each data set's
    cluster, OPTICS's ordering and reachability, the coordinates' principal
    components, the dimension functional."""
    import json
    import time

    if not os.path.isfile(path):
        return path, None
    with open(path) as f:
        doc = json.load(f)
    doc["program"] = "mxi_correlation_json"
    doc["time"] = time.strftime(
        "%Y-%m-%dT%H:%M:%S", time.localtime(os.path.getmtime(path))
    )
    return path, doc


def is_table(value):
    return isinstance(value, dict) and isinstance(value.get("columns"), list)


def column(table, name):
    return list(table.get(name, [])) if is_table(table) else []


def _number(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


# --- facts ----------------------------------------------------------------

LABELS = {
    "spots": "Spots",
    "signal_pixels": "Signal pixels",
    "indexed": "Indexed",
    "indexed_fraction": "Indexed fraction",
    "lattices": "Lattices",
    "reflections": "Reflections",
    "refined_on": "Refined on",
    "rejected": "Rejected as outliers",
    "rmsd_px_px_images": "RMSD (x px, y px, z images)",
    "scan_varying": "Scan-varying",
    "cell": "Cell",
    "volume": "Volume (Å³)",
    "d_min": "Resolution limit (Å)",
    "max_cell": "Longest cell edge (Å)",
    "laue_group": "Laue group",
    "space_group": "Space group",
    "reindex": "Reindexed by",
    "observations": "Observations",
    "unique": "Unique",
    "outliers": "Outliers",
    "resolution_cc_half": "Resolution at CC½ 0.3 (Å)",
    "resolution_cc_half_significance": "Resolution, CC½ significant (Å)",
    "error_model_a": "Error model a",
    "error_model_b": "Error model b",
    "two_theta_cell": "Cell from 2θ",
    "two_theta_cell_esd": "Its e.s.d.s",
    "ambiguous": "Ambiguous spots",
    # A cluster's (mxi correlation --write-clusters, summary_cut.hh):
    "cluster": "Cluster",
    "data_sets": "Data sets",
    "of_data_sets": "Of",
    "cut_from": "Cut from",
    "whole_run_left_out": "Left out (the whole run's)",
}

SKIP = {
    "program",
    "time",
    "first_image",
    "spots_per_image",
    "cell_over_scan",
    "each_sweep",
    "geometry",
    "geometry_given",
    "profile_model",
    "elements",
    "ice_rings",
}


def _one(v):
    """A number as a person reads it: whole where it is large or whole, else
    four figures -- not 2.361e+05 for a volume."""
    f = float(v)
    if abs(f) >= 1000 or (f.is_integer() and abs(f) < 1e15):
        return f"{f:.0f}"
    return f"{f:.4g}"


def _text(value, key=""):
    if isinstance(value, bool):
        return "yes" if value else "no"
    if key in ("cell", "two_theta_cell", "two_theta_cell_esd") and isinstance(
        value, list
    ):
        # A cell as mxi show gives it: three decimals, lengths and angles.
        return "  ".join(f"{float(v):.3f}" for v in value if _number(v) is not None)
    if isinstance(value, (int, float)) and _number(value) is not None:
        return _one(value)
    if isinstance(value, list) and all(_number(v) is not None for v in value):
        return "  ".join(_one(v) for v in value)
    if isinstance(value, list) and all(isinstance(v, str) for v in value):
        return ", ".join(value)
    return str(value)


def facts(rec):
    """[(label, text)] for the record's figures: those it knows first, in
    order, then any other figure it holds that is not a table or a series."""
    if not rec:
        return []
    if rec.get("program") == "mxi_correlation_json":
        return correlation_facts(rec)
    out = []
    for key, label in LABELS.items():
        if key in rec and not is_table(rec[key]):
            value = rec[key]
            if key == "indexed_fraction" and _number(value) is not None:
                out.append((label, f"{100 * float(value):.1f}%"))
            elif key == "scan_varying":
                out.append((label, "yes" if value else "no"))
            else:
                out.append((label, _text(value, key)))
    for key in sorted(rec):
        if key in LABELS or key in SKIP or is_table(rec[key]):
            continue
        value = rec[key]
        if isinstance(value, (dict, list)) and not (
            isinstance(value, list)
            and len(value) <= 8
            and all(_number(v) is not None for v in value)
        ):
            continue
        out.append((key.replace("_", " ").capitalize(), _text(value)))
    return out


def correlation_facts(rec):
    labels = [int(v) for v in rec.get("labels", [])]
    clusters = sorted({v for v in labels if v >= 0})
    out = [("Data sets", str(len(rec.get("experiments", labels))))]
    out.append(("Clusters", str(len(clusters))))
    sizes = ", ".join(f"{c}: {labels.count(c)}" for c in clusters)
    if sizes:
        out.append(("Their data sets", sizes))
    out.append(("In no cluster", str(labels.count(-1))))
    if "dimension" in rec:
        out.append(("Dimensions", str(rec["dimension"])))
    if "min_samples" in rec:
        out.append(("OPTICS min_samples", str(rec["min_samples"])))
    return out


# --- tables ---------------------------------------------------------------

MERGING_ROWS = (
    ("d_min", "High resolution limit", "{:.2f}"),
    ("d_max", "Low resolution limit", "{:.2f}"),
    ("completeness", "Completeness", "{:.1f}"),
    ("multiplicity", "Multiplicity", "{:.1f}"),
    ("i_over_sigma", "I/sigma", "{:.1f}"),
    ("rmerge", "Rmerge(I)", "{:.3f}"),
    ("rmeas", "Rmeas(I)", "{:.3f}"),
    ("rpim", "Rpim(I)", "{:.3f}"),
    ("cc_half", "CC half", "{:.3f}"),
    ("observations", "Observations", "{:.0f}"),
    ("unique", "Unique", "{:.0f}"),
    ("anom_completeness", "Anomalous completeness", "{:.1f}"),
    ("anom_multiplicity", "Anomalous multiplicity", "{:.1f}"),
    ("cc_anom", "CC anomalous", "{:.3f}"),
)


def _formatted(fmt, v):
    n = _number(v)
    return fmt.format(n) if n is not None else ("" if v is None else str(v))


def tables(rec):
    """The record's tables worth reading as tables: scaling's merging summary
    as DIALS lays it out (a row a figure, a column overall, low, high), the
    subgroups symmetry weighed, integration's overall figures."""
    if not rec:
        return []
    out = []
    summary = rec.get("merging_summary")
    if is_table(summary):
        shells = column(summary, "shell") or ["overall", "low", "high"]
        rows = [
            [label] + [_formatted(fmt, v) for v in column(summary, key)]
            for key, label, fmt in MERGING_ROWS
            if column(summary, key)
            and any(_number(v) is not None for v in column(summary, key))
        ]
        out.append(Table("Merging statistics", [""] + [str(s) for s in shells], rows))
    for name, title in (
        ("subgroups", "Subgroups weighed"),
        ("overall", "Overall"),
        ("candidates", "Candidate cells"),
    ):
        t = rec.get(name)
        if is_table(t):
            names = t["columns"]
            n = max((len(column(t, c)) for c in names), default=0)
            rows = [
                [_cell(column(t, c)[i]) if i < len(column(t, c)) else "" for c in names]
                for i in range(n)
            ]
            out.append(Table(title, list(names), rows))
    return out


def _cell(v):
    if isinstance(v, list) and v and all(_number(x) is not None for x in v):
        # A table's cell that is itself a list -- each candidate's unit cell --
        # compact: printed whole, at full precision, it had pushed the other
        # columns out of sight (Graeme's indexed.rflx).
        return " ".join(f"{float(x):.2f}" for x in v)
    n = _number(v)
    if n is None:
        return "" if v is None else str(v)
    return _one(n)  # as the facts: large or whole numbers whole


# --- plots ----------------------------------------------------------------


def _resolution_x(table):
    """Each shell's middle in 1/d^2, from its d_max and d_min."""
    out = []
    for hi, lo in zip(column(table, "d_max"), column(table, "d_min")):
        a, b = _number(hi), _number(lo)
        inv_hi = 1.0 / (a * a) if a else 0.0
        inv_lo = 1.0 / (b * b) if b else None
        out.append(None if inv_lo is None else 0.5 * (inv_hi + inv_lo))
    return out


def _by_resolution(table, names, title, ylabel, log=False):
    x = _resolution_x(table)
    series = [
        Series(label, x, [_number(v) for v in column(table, key)])
        for key, label in names
        if column(table, key)
    ]
    return (
        Plot(title, "Resolution (Å)", ylabel, series, resolution=True, log=log)
        if series
        else None
    )


def plots(rec):
    """The record's plots, as mxi report draws them."""
    if not rec:
        return []
    program = rec.get("program", "")
    made = []
    if program == "mxi_find":
        per_image = rec.get("spots_per_image") or []
        firsts = rec.get("first_image") or []
        series = []
        for i, counts in enumerate(per_image):
            first = (
                int(firsts[i])
                if i < len(firsts) and _number(firsts[i]) is not None
                else 1
            )
            series.append(
                Series(
                    f"sweep {i}", list(range(first, first + len(counts))), list(counts)
                )
            )
        if series:
            made.append(Plot("Spots per image", "Image", "Spots", series, integer=True))
        t = rec.get("spots_by_resolution")
        if is_table(t):
            made.append(
                Plot(
                    "Spots by resolution",
                    "Resolution (Å)",
                    "Spots",
                    [Series("spots", _resolution_x(t), column(t, "spots"))],
                    resolution=True,
                    bars=True,
                )
            )
    elif program in ("mxi_index", "mxi_refine"):
        t = rec.get("cycles")
        if is_table(t):
            x = list(range(1, len(column(t, "cycle")) + 1))
            made.append(
                Plot(
                    "RMSDs, cycle by cycle",
                    "Cycle",
                    "RMSD",
                    integer=True,
                    series=[
                        Series("x (px)", x, column(t, "rmsd_x_px")),
                        Series("y (px)", x, column(t, "rmsd_y_px")),
                        Series("z (images)", x, column(t, "rmsd_z_images")),
                    ],
                )
            )
        for one in rec.get("cell_over_scan") or []:
            if isinstance(one, dict) and isinstance(one.get("angle"), list):
                tag = f"experiment {one.get('experiment', '')}"
                made.append(
                    Plot(
                        f"Cell lengths over the scan, {tag}",
                        "Rotation (°)",
                        "Å",
                        [
                            Series(k, one["angle"], one.get(k, []))
                            for k in ("a", "b", "c")
                        ],
                    )
                )
                made.append(
                    Plot(
                        f"Cell angles over the scan, {tag}",
                        "Rotation (°)",
                        "°",
                        [
                            Series(k, one["angle"], one.get(k, []))
                            for k in ("alpha", "beta", "gamma")
                        ],
                    )
                )
    elif program == "mxi_correlation_json":
        made.extend(correlation_plots(rec))
    elif program == "mxi_integrate":
        made.extend(each_sweep_plots(rec))
        t = rec.get("by_resolution")
        if is_table(t):
            for p in (
                _by_resolution(
                    t,
                    [("isigma_sum", "summed"), ("isigma_prf", "profile-fitted")],
                    "<I/σ> by resolution",
                    "<I/σ>",
                ),
                _by_resolution(
                    t,
                    [("cc_prf", "profile correlation")],
                    "Profile correlation by resolution",
                    "CC",
                ),
                _by_resolution(
                    t, [("ibg", "background")], "Background by resolution", "Counts"
                ),
            ):
                if p:
                    made.append(p)
    elif program == "mxi_scale":
        t = rec.get("merging_by_resolution")
        if is_table(t):
            for p in (
                _by_resolution(
                    t,
                    [("cc_half", "CC½"), ("cc_anom", "CC anomalous")],
                    "CC½ by resolution",
                    "CC",
                ),
                _by_resolution(
                    t, [("i_over_sigma", "<I/σ>")], "<I/σ> by resolution", "<I/σ>"
                ),
            ):
                if p:
                    made.append(p)
            x = _resolution_x(t)
            if column(t, "completeness"):
                made.append(
                    Plot(
                        "Completeness and multiplicity",
                        "Resolution (Å)",
                        "Completeness (%)",
                        [Series("completeness", x, column(t, "completeness"))],
                        resolution=True,
                        right=[Series("multiplicity", x, column(t, "multiplicity"))],
                        right_label="Multiplicity",
                    )
                )
            r = _by_resolution(
                t,
                [("rmerge", "Rmerge"), ("rmeas", "Rmeas"), ("rpim", "Rpim")],
                "R factors by resolution",
                "R",
                # Linear, as R factors are usually drawn: on a log axis,
                # as mxi report has them, they looked odd here (Graeme).
            )
            if r:
                made.append(r)
    return made


def cosine_angles(coordinates, n):
    """The cosine of the angle between every pair of data sets' cosym
    coordinates, as mxi correlation's report computes it: the coordinates
    dimension by dimension, x[k * n + i]; 0 where a vector has no length."""
    if n == 0 or len(coordinates) < n or len(coordinates) % n:
        return []
    dim = len(coordinates) // n
    x = [[float(coordinates[k * n + i]) for k in range(dim)] for i in range(n)]
    length = [math.sqrt(sum(v * v for v in row)) for row in x]
    return [
        [
            (
                sum(a * b for a, b in zip(x[i], x[j])) / (length[i] * length[j])
                if length[i] > 0 and length[j] > 0
                else 0.0
            )
            for j in range(n)
        ]
        for i in range(n)
    ]


def _matrix(title, full, order):
    """A matrix in the given order, as a plot: NaN where not a number, the
    scale from its least (or 0) to 1."""
    matrix = [
        [full[i][j] if full[i][j] is not None else float("nan") for j in order]
        for i in order
    ]
    values = [v for row in matrix for v in row if v == v]
    return Plot(
        title,
        "Data set, in order",
        "Data set, in order",
        [],
        kind="matrix",
        matrix=matrix,
        low=min(0.0, min(values, default=0.0)),
    )


MANY_SWEEPS = 8  # more, and a line each is a tangle: grey, a median, bars


def each_sweep_plots(rec):
    """Integration a sweep at a time, where there are several -- or where only
    they are left (a cluster's cut summary): so that a sweep unlike the
    others shows. Up to MANY_SWEEPS, a line each; more -- 36 crystals' -- each
    a thin grey line under their median, and a bar each for its overall
    <I/σ> and profile correlation, which shows which sweep is the odd one: a
    line each, with a legend of 36, had buried the plot (Graeme)."""
    each = [e for e in rec.get("each_sweep") or [] if isinstance(e, dict)]
    if not each or (len(each) < 2 and is_table(rec.get("by_resolution"))):
        return []
    if len(each) > MANY_SWEEPS:
        return _many_sweeps(each)
    made = []
    for key, title, ylabel in (
        (("isigma_prf", "isigma_sum"), "<I/σ> by resolution, each sweep", "<I/σ>"),
        (("cc_prf",), "Profile correlation by resolution, each sweep", "CC"),
    ):
        series = []
        for e in each:
            t = e.get("by_resolution")
            name = next((k for k in key if column(t, k)), None)
            if is_table(t) and name:
                series.append(
                    Series(
                        f"sweep {e.get('sweep', len(series))}",
                        _resolution_x(t),
                        [_number(v) for v in column(t, name)],
                    )
                )
        if series:
            made.append(Plot(title, "Resolution (Å)", ylabel, series, resolution=True))
    return made


def _many_sweeps(each):
    made = []
    lines = []
    for e in each:
        t = e.get("by_resolution")
        name = next((k for k in ("isigma_prf", "isigma_sum") if column(t, k)), None)
        if is_table(t) and name:
            lines.append(
                Series(
                    f"sweep {e.get('sweep')}",
                    _resolution_x(t),
                    [_number(v) for v in column(t, name)],
                    thin=True,
                )
            )
    if lines:
        series = list(lines)
        # The median, shell by shell, where the sweeps share their shells.
        if len({len(s.x) for s in lines}) == 1:
            xs, ys = [], []
            for k in range(len(lines[0].x)):
                at = [s.x[k] for s in lines if s.x[k] is not None]
                vs = sorted(s.y[k] for s in lines if s.y[k] is not None)
                if at and vs:
                    xs.append(sum(at) / len(at))
                    mid = len(vs) // 2
                    ys.append(vs[mid] if len(vs) % 2 else (vs[mid - 1] + vs[mid]) / 2)
            series.append(Series("median", xs, ys))
        made.append(
            Plot(
                f"<I/σ> by resolution, {len(lines)} sweeps",
                "Resolution (Å)",
                "<I/σ>",
                series,
                resolution=True,
            )
        )
    for keys, title, ylabel in (
        (("isigma_prf", "isigma_sum"), "<I/σ>, each sweep", "<I/σ>"),
        (("cc_prf",), "Profile correlation, each sweep", "CC"),
    ):
        xs, ys = [], []
        for k, e in enumerate(each):
            t = e.get("overall")
            name = next((n for n in keys if column(t, n)), None)
            if is_table(t) and name:
                xs.append(int(e.get("sweep", k)))
                ys.append(_number(column(t, name)[0]))
        if xs:
            made.append(
                Plot(
                    title,
                    "Sweep",
                    ylabel,
                    [Series(ylabel, xs, ys)],
                    bars=True,
                    integer=True,
                )
            )
    return made


def correlation_plots(rec):
    """mxi correlation's, as its HTML report draws them: the matrix in
    OPTICS's order, the reachability coloured by cluster, the data sets on
    their first two principal components, the functional by dimension."""
    made = []
    labels = [int(v) for v in rec.get("labels", [])]
    order = [int(v) for v in rec.get("ordering", [])] or list(range(len(labels)))
    r = rec.get("r") or []
    n = len(r)
    if r and len(order) == n:
        # The diagonal 1, as the report draws it; what is not a number, not
        # drawn (NaN), where matplotlib would refuse a None.
        full = [
            [1.0 if i == j else _number(r[i][j]) for j in range(n)] for i in range(n)
        ]
        made.append(_matrix("Correlation, in OPTICS's order", full, order))
        cosine = cosine_angles(rec.get("coordinates") or [], n)
        if cosine:
            made.append(
                _matrix(
                    "Cosine angles, in OPTICS's order",
                    cosine,
                    order,
                )
            )
    reach = [_number(v) for v in rec.get("reachability_in_order", [])]
    if reach and len(reach) == len(order) and labels:
        finite = [v for v in reach if v is not None]
        top = max(finite) * 1.1 if finite else 1.0
        series = []
        for c in sorted(set(labels)):
            xs = [k for k, i in enumerate(order) if labels[i] == c]
            ys = [reach[k] if reach[k] is not None else top for k in xs]
            series.append(Series("in no cluster" if c < 0 else f"cluster {c}", xs, ys))
        made.append(
            Plot(
                "Reachability, in OPTICS's order",
                "Data set, in order",
                "Reachability",
                series,
                bars=True,
            )
        )
    comps = rec.get("components") or []
    if comps and labels and all(isinstance(c, list) and len(c) >= 2 for c in comps):
        ratio = rec.get("explained_variance_ratio") or []
        axis = [
            f"PC{k + 1}" + (f" ({100 * ratio[k]:.0f}%)" if k < len(ratio) else "")
            for k in range(2)
        ]
        series = []
        for c in sorted(set(labels)):
            pts = [
                comps[i] for i in range(min(len(comps), len(labels))) if labels[i] == c
            ]
            series.append(
                Series(
                    "in no cluster" if c < 0 else f"cluster {c}",
                    [p[0] for p in pts],
                    [p[1] for p in pts],
                )
            )
        made.append(
            Plot(
                "The data sets, principal components",
                axis[0],
                axis[1],
                series,
                kind="points",
            )
        )
    dims, func = rec.get("scan_dimensions") or [], rec.get("scan_functional") or []
    if dims and len(dims) == len(func):
        made.append(
            Plot(
                "Functional by dimension",
                "Dimensions",
                "Functional",
                [Series("functional", dims, func)],
                integer=True,
            )
        )
    return made


def d_ticks(x_values, count=6):
    """Ticks for a resolution axis: positions in 1/d^2, labels in d."""
    xs = [x for x in x_values if x is not None and x > 0]
    if not xs:
        return [], []
    lo, hi = min(xs), max(xs)
    positions = (
        [lo + (hi - lo) * i / (count - 1) for i in range(count)] if hi > lo else [lo]
    )
    return positions, [f"{1 / math.sqrt(p):.2f}" for p in positions]
