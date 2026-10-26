"""Each scaled observation against a refined structure's Fc^2: the processing
judged against something close to the truth, not against its own equivalents
or another program.

A small molecule refined to an R1 of a few per cent has Fc^2 that are, to that
accuracy, what the intensities should have been. Unmerged observations from a
scaled .rflx are matched to SHELXL's .fcf (LIST 4: Fc^2, Fo^2, sigma) by their
symmetry-unique reflection, Friedel mates apart unless asked, one scale k fitted
between I and Fc^2, and each observation's z = (I - k Fc^2) / sigma studied:

* its calibration: z's variance should be 1, and its tails Gaussian's;
* binned by what processing could get wrong -- intensity, resolution, image,
  sweep, partiality, detector position, background, profile against summation;
* split into what all observations of a reflection share (the model's part:
  extinction, anomalous scattering, hydrogens) and what one observation does
  apart from its own equivalents (the processing's part);
* the worst observations, each with what is known of it.

Independent of mxi, as all of mxeq: it reads files.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from itertools import permutations

import numpy as np

from mxeq.equivalents import asu_keys

# --- reading -----------------------------------------------------------------


@dataclass
class Fcf:
    """A SHELXL .fcf: its cell, space group and reflections."""

    cell: tuple  # a, b, c, alpha, beta, gamma
    hall: str
    hkl: np.ndarray  # (n, 3) int
    fc2: np.ndarray
    fo2: np.ndarray
    sigma: np.ndarray
    name: str = ""


def read_fcf(path: str) -> Fcf:
    """A SHELXL .fcf written with LIST 4 (h k l, Fc^2, Fo^2, sigma(Fo^2)); the
    space group from its own symmetry operations where it lists them, its name
    otherwise."""
    import gemmi

    block = gemmi.cif.read_file(path).sole_block()

    def number(tag):
        v = block.find_value(tag)
        if v is None:
            raise ValueError(f"{path}: no {tag}")
        return float(re.sub(r"\(.*\)", "", v))

    cell = tuple(
        number(t)
        for t in (
            "_cell_length_a",
            "_cell_length_b",
            "_cell_length_c",
            "_cell_angle_alpha",
            "_cell_angle_beta",
            "_cell_angle_gamma",
        )
    )
    ops = []
    for tag in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"):
        col = block.find_values(tag)
        if col:
            ops = [gemmi.cif.as_string(v) for v in col]
            break
    if ops:
        group = gemmi.find_spacegroup_by_ops(gemmi.GroupOps([gemmi.Op(o) for o in ops]))
        if group is None:
            raise ValueError(f"{path}: its symmetry operations name no space group")
    else:
        name = None
        for tag in ("_space_group_name_H-M_alt", "_symmetry_space_group_name_H-M"):
            v = block.find_value(tag)
            if v:
                name = gemmi.cif.as_string(v)
                break
        if not name:
            raise ValueError(f"{path}: neither symmetry operations nor a space group")
        group = gemmi.SpaceGroup(name)
    table = block.find(
        "_refln_",
        [
            "index_h",
            "index_k",
            "index_l",
            "F_squared_calc",
            "F_squared_meas",
            "F_squared_sigma",
        ],
    )
    if not table:
        raise ValueError(
            f"{path}: no _refln_ loop with F_squared_calc, F_squared_meas and "
            "F_squared_sigma -- written with LIST 4?"
        )
    rows = np.array([[float(x) for x in row] for row in table], dtype=float)
    return Fcf(
        cell=cell,
        hall=group.hall,
        hkl=rows[:, :3].astype(np.int64),
        fc2=rows[:, 3],
        fo2=rows[:, 4],
        sigma=rows[:, 5],
        name=group.hm,
    )


def read_res_extinction(path: str) -> float | None:
    """EXTI from a SHELXL .res, if refined: the strong reflections' Fc^2 in
    the .fcf may or may not carry its correction, so it is reported."""
    with open(path) as f:
        for line in f:
            words = line.split()
            if words and words[0].upper() == "EXTI" and len(words) > 1:
                return float(words[1])
    return None


@dataclass
class Observations:
    """Scaled, unmerged observations, as read from a .rflx: what the analysis
    needs, and whatever else is known of each (None where a column is
    missing)."""

    hkl: np.ndarray
    intensity: np.ndarray
    variance: np.ndarray
    cell: tuple
    hall: str
    d: np.ndarray | None = None
    image: np.ndarray | None = None
    x: np.ndarray | None = None
    y: np.ndarray | None = None
    experiment: np.ndarray | None = None
    partiality: np.ndarray | None = None
    background: np.ndarray | None = None
    prf: np.ndarray | None = None
    sum: np.ndarray | None = None
    scale: np.ndarray | None = None
    rejected: np.ndarray | None = None  # rejected as an outlier by scaling
    names: dict = field(default_factory=dict)  # experiment -> its images' file
    zeta: np.ndarray | None = None  # nearness to the rotation axis, 0 on it


# DIALS's scaling flags, as mxi writes them (src/scale.hh: kOutlierInScaling,
# kExcludedForScaling, kScaled) and dials.correlation_matrix reads them back.
# Positions 36 to 40 had been guessed here at first, so that neither the
# outliers were marked nor the excluded rows left out (21 October).
_OUTLIER_IN_SCALING = 1 << 23
_EXCLUDED_FOR_SCALING = 1 << 24
_SCALED = 1 << 26


def _cell_of(crystal: dict) -> tuple:
    a = np.array(crystal["real_space_a"], float)
    b = np.array(crystal["real_space_b"], float)
    c = np.array(crystal["real_space_c"], float)

    def angle(u, v):
        return math.degrees(
            math.acos(float(u @ v / np.linalg.norm(u) / np.linalg.norm(v)))
        )

    return (
        float(np.linalg.norm(a)),
        float(np.linalg.norm(b)),
        float(np.linalg.norm(c)),
        angle(b, c),
        angle(a, c),
        angle(a, b),
    )


def read_scaled(path: str) -> Observations:
    """A scaled .rflx: the scaled intensities and their variances, rows
    scaling excluded left out, rows it rejected as outliers kept and marked."""
    from mxeq import rflx

    table = rflx.load_table(path)
    experiments = rflx.experiments_dict(path)
    crystals = experiments.get("crystal", [])
    if not crystals:
        raise ValueError(f"{path}: no crystal in its experiment list")
    first = crystals[0]
    # Each experiment named by its images' file: an experiment's number is not
    # its sweep's -- joint indexing puts the lattices after all the sweeps, so
    # l-cyst_01's lattice is experiment 4 (Graeme, 21 October).
    imagesets = experiments.get("imageset", [])
    names = {}
    for i, e in enumerate(experiments.get("experiment", [])):
        k = e.get("imageset")
        if k is not None and k < len(imagesets):
            path = imagesets[k].get("template") or imagesets[k].get("master") or ""
            names[i] = str(path).split("/")[-1]
    hall = first.get("space_group_hall_symbol", " P 1")
    for name in ("intensity.scale.value", "intensity.scale.variance", "miller_index"):
        if name not in table.columns:
            raise ValueError(f"{path}: no {name} -- is it scaled?")
    flags = (
        np.asarray(table["flags"], dtype=np.int64) if "flags" in table.columns else None
    )
    keep = np.ones(table.nrows, bool)
    if flags is not None:
        keep &= (flags & _EXCLUDED_FOR_SCALING) == 0
    hkl = np.asarray(table["miller_index"], dtype=np.int64)
    keep &= np.any(hkl != 0, axis=1)
    var = np.asarray(table["intensity.scale.variance"], float)
    keep &= np.isfinite(var) & (var > 0)

    def col(name, k=None):
        if name not in table.columns:
            return None
        v = np.asarray(table[name])
        v = v[:, k] if k is not None else v
        return v[keep]

    # DIALS's convention, which mxi keeps: intensity.scale.value is the
    # corrected intensity before scaling, and the scaled one is it divided by
    # inverse_scale_factor (as dials.export writes it). Read without the
    # division at first, the scales' own variation -- +/-20 per cent -- had
    # looked like the data's disagreement with Fc^2 (21 October).
    g = (
        np.asarray(table["inverse_scale_factor"], float)
        if "inverse_scale_factor" in table.columns
        else np.ones(table.nrows)
    )
    keep &= np.isfinite(g) & (g > 0)
    return Observations(
        hkl=hkl[keep],
        intensity=(np.asarray(table["intensity.scale.value"], float) / g)[keep],
        variance=(var / g**2)[keep],
        cell=_cell_of(first),
        hall=hall,
        d=col("d"),
        image=(
            col("xyzobs.px.value", 2)
            if "xyzobs.px.value" in table.columns
            else col("xyzcal.px", 2)
        ),
        x=(
            col("xyzobs.px.value", 0)
            if "xyzobs.px.value" in table.columns
            else col("xyzcal.px", 0)
        ),
        y=(
            col("xyzobs.px.value", 1)
            if "xyzobs.px.value" in table.columns
            else col("xyzcal.px", 1)
        ),
        experiment=col("id"),
        partiality=col("partiality"),
        background=col("background.mean"),
        prf=col("intensity.prf.value"),
        sum=col("intensity.sum.value"),
        scale=col("inverse_scale_factor"),
        rejected=(
            ((flags[keep] & _OUTLIER_IN_SCALING) != 0) if flags is not None else None
        ),
        names=names,
        zeta=col("zeta"),
    )


# --- matching ----------------------------------------------------------------


def parse_operator(text: str) -> np.ndarray:
    """A change of basis as dials writes it, "b,c,a" or "-a,c,b": the new axes
    in the old, as rows; Miller indices transform with it, h_new = M h."""
    rows = []
    for part in text.replace(" ", "").split(","):
        row = np.zeros(3, dtype=np.int64)
        for sign, letter in re.findall(r"([+-]?)([abc])", part):
            row["abc".index(letter)] += -1 if sign == "-" else 1
        rows.append(row)
    if len(rows) != 3:
        raise ValueError(f"not a change of basis: {text!r}")
    return np.array(rows)


def setting_from_cells(cell: tuple, target: tuple, tolerance: float = 0.02):
    """The permutation of axes taking `cell` to `target` (lengths within
    `tolerance`, relative; angles carried with them), handedness kept by a sign;
    None if none does. The identity when they already agree."""
    lengths, angles = np.array(cell[:3]), np.array(cell[3:])
    goal_l, goal_a = np.array(target[:3]), np.array(target[3:])
    best = None
    for perm in permutations(range(3)):
        l = lengths[list(perm)]
        a = angles[list(perm)]
        if np.all(np.abs(l - goal_l) <= tolerance * goal_l) and np.all(
            np.abs(a - goal_a) <= 2.0
        ):
            m = np.zeros((3, 3), dtype=np.int64)
            for new, old in enumerate(perm):
                m[new, old] = 1
            if round(np.linalg.det(m)) < 0:
                m[2] *= -1  # right-handed: an odd permutation takes a sign
            score = float(np.sum(np.abs(l - goal_l)))
            if best is None or score < best[0]:
                best = (score, m)
    return None if best is None else best[1]


# --- the analysis ------------------------------------------------------------


@dataclass
class Bin:
    label: str
    n: int
    mean_z: float
    var_z: float
    robust_var_z: float
    median_rel: float  # median of (I - k Fc^2) / (k Fc^2), where k Fc^2 > 0
    over_5: float  # fraction with |z| > 5


@dataclass
class Result:
    matched: int
    unmatched: int
    k: float
    b: float
    r_unmerged: float  # sum |I - k Fc^2| / sum |I|
    var_z: float
    robust_var_z: float
    tails: dict  # threshold -> (observed fraction, Gaussian fraction)
    between: float  # share of sum z^2 from reflections' mean offsets
    within: float  # share from observations about their own reflection's mean
    # z about its own reflection's mean, scaled by n / (n - 1), reflections
    # seen three times or more: the sigmas against the processing alone, the
    # model's share cancelled. Against Fc^2 the variance also holds the
    # model's errors -- 55 to 65 per cent of the disagreement on a small
    # molecule -- so that sigmas right for the data read as too small there
    # (23 October). This is the calibration to judge an error model by.
    robust_var_within: float = float("nan")
    bins: dict = field(default_factory=dict)  # name -> list[Bin]
    worst: list = field(default_factory=list)
    rejected: dict = field(default_factory=dict)
    z: np.ndarray | None = None
    zbar: np.ndarray | None = None
    model: np.ndarray | None = None
    rows: np.ndarray | None = None  # matched observations' indices


def _robust_var(z):
    z = z[np.isfinite(z)]
    if len(z) == 0:
        return float("nan")
    return float((1.4826 * np.median(np.abs(z - np.median(z)))) ** 2)


def _bin(label, z, rel):
    ok = np.isfinite(z)
    zz, rr = z[ok], rel[ok]
    rr = rr[np.isfinite(rr)]
    return Bin(
        label=label,
        n=int(len(zz)),
        mean_z=float(np.mean(zz)) if len(zz) else float("nan"),
        var_z=float(np.var(zz)) if len(zz) else float("nan"),
        robust_var_z=_robust_var(zz),
        median_rel=float(np.median(rr)) if len(rr) else float("nan"),
        over_5=float(np.mean(np.abs(zz) > 5)) if len(zz) else float("nan"),
    )


def binned(values, z, rel, bins=10, fmt="{:.3g}"):
    """Equal-count bins of `values`, each one's z and relative difference."""
    values = np.asarray(values, float)
    ok = np.isfinite(values)
    if ok.sum() < bins:
        return []
    order = np.argsort(values[ok])
    idx = np.where(ok)[0][order]
    out = []
    for part in np.array_split(idx, bins):
        lo, hi = values[part[0]], values[part[-1]]
        out.append(_bin(f"{fmt.format(lo)} to {fmt.format(hi)}", z[part], rel[part]))
    return out


def fit_scale(intensity, sigma, model, s2=None, fit_b=False, clip=10.0, rounds=4):
    """k (and B, I = k exp(-2 B s^2) Fc^2, s = 1 / 2d) by weighted least squares,
    rows more than `clip` sigma off left out after the first round -- the
    outliers are what is being looked for, and must not set the scale."""
    w = 1.0 / sigma**2
    use = np.isfinite(model) & (model > 0) & np.isfinite(intensity)
    k, b = 1.0, 0.0
    for _ in range(rounds):
        f = model * (np.exp(-2.0 * b * s2) if (fit_b and s2 is not None) else 1.0)
        k = float(np.sum((w * intensity * f)[use]) / np.sum((w * f * f)[use]))
        if fit_b and s2 is not None:
            # B from the log of the ratio, against s^2, the strong rows
            strong = use & (intensity > 3 * sigma) & (model > 0)
            if strong.sum() > 20:
                y = np.log(np.clip(intensity[strong] / (k * model[strong]), 1e-6, None))
                slope = np.polyfit(s2[strong], y, 1)[0]
                b += -slope / 2.0
        f = model * (np.exp(-2.0 * b * s2) if (fit_b and s2 is not None) else 1.0)
        z = (intensity - k * f) / sigma
        use = np.isfinite(model) & (model > 0) & (np.abs(z) < clip)
    return k, b


def analyse(
    obs: Observations,
    fcf: Fcf,
    reindex: np.ndarray | None = None,
    anomalous: bool = True,
    fit_b: bool = False,
    bins: int = 10,
    worst: int = 30,
) -> Result:
    """The observations against the .fcf's Fc^2 -- see the module's note."""
    hkl = obs.hkl if reindex is None else obs.hkl @ reindex.T
    keys_obs = asu_keys(hkl, fcf.hall, anomalous=anomalous)
    keys_fcf = asu_keys(fcf.hkl, fcf.hall, anomalous=anomalous)
    lookup = {int(k): i for i, k in enumerate(keys_fcf)}
    where = np.array([lookup.get(int(k), -1) for k in keys_obs])
    rows = np.where(where >= 0)[0]
    fc2 = fcf.fc2[where[rows]]
    i_obs = obs.intensity[rows]
    sigma = np.sqrt(obs.variance[rows])
    s2 = None
    if obs.d is not None:
        s2 = 1.0 / (4.0 * obs.d[rows] ** 2)
    k, b = fit_scale(i_obs, sigma, fc2, s2, fit_b)
    model = k * fc2 * (np.exp(-2.0 * b * s2) if (fit_b and s2 is not None) else 1.0)
    z = (i_obs - model) / sigma
    rel = np.where(model > 0, (i_obs - model) / np.where(model > 0, model, 1.0), np.nan)

    # What a reflection's observations share against what each does alone.
    keys = keys_obs[rows]
    uniq, inverse = np.unique(keys, return_inverse=True)
    count = np.bincount(inverse)
    zbar = (np.bincount(inverse, weights=z) / count)[inverse]
    multi = count[inverse] >= 2
    total = float(np.sum(z[multi] ** 2))
    between = float(np.sum(zbar[multi] ** 2)) / total if total > 0 else float("nan")
    within = (
        float(np.sum((z[multi] - zbar[multi]) ** 2)) / total
        if total > 0
        else float("nan")
    )

    nk = count[inverse]
    many = nk >= 3
    within_z = (z - zbar) * np.sqrt(nk / np.maximum(nk - 1, 1))
    robust_within = _robust_var(within_z[many]) if many.any() else float("nan")

    tails = {}
    for t in (3, 5, 10):
        gauss = math.erfc(t / math.sqrt(2.0))
        tails[t] = (float(np.mean(np.abs(z) > t)), gauss)

    res = Result(
        matched=int(len(rows)),
        unmatched=int(len(obs.hkl) - len(rows)),
        k=k,
        b=b,
        r_unmerged=float(np.sum(np.abs(i_obs - model)) / np.sum(np.abs(i_obs))),
        var_z=float(np.var(z)),
        robust_var_z=_robust_var(z),
        tails=tails,
        between=between,
        within=within,
        robust_var_within=robust_within,
        z=z,
        zbar=zbar,
        model=model,
        rows=rows,
    )

    def sub(v):
        return None if v is None else np.asarray(v)[rows]

    res.bins["k Fc^2 (intensity)"] = binned(model, z, rel, bins)
    if obs.d is not None:
        res.bins["d (A)"] = binned(sub(obs.d), z, rel, bins)
    if obs.image is not None:
        res.bins["image (z, px)"] = binned(sub(obs.image), z, rel, bins, "{:.0f}")
    if obs.experiment is not None:
        e = sub(obs.experiment)
        res.bins["experiment"] = [
            _bin(
                f"{int(x)} {obs.names.get(int(x), '')}".strip()[:24],
                z[e == x],
                rel[e == x],
            )
            for x in np.unique(e)
        ]
    if obs.partiality is not None:
        res.bins["partiality"] = binned(sub(obs.partiality), z, rel, min(bins, 5))
    if obs.x is not None and obs.y is not None:
        x, y = sub(obs.x), sub(obs.y)
        cx, cy = np.median(x), np.median(y)
        res.bins["distance from the spots' centre (px)"] = binned(
            np.hypot(x - cx, y - cy), z, rel, bins, "{:.0f}"
        )
    if obs.background is not None:
        res.bins["background (counts/px)"] = binned(sub(obs.background), z, rel, bins)
    if obs.prf is not None and obs.sum is not None:
        p, s_ = sub(obs.prf), sub(obs.sum)
        ratio = np.where(np.abs(s_) > 0, p / np.where(np.abs(s_) > 0, s_, 1.0), np.nan)
        res.bins["profile / summed"] = binned(ratio, z, rel, bins)
    if obs.scale is not None:
        res.bins["inverse scale factor"] = binned(sub(obs.scale), z, rel, bins)
    if obs.zeta is not None:
        res.bins["|zeta| (0 on the rotation axis)"] = binned(
            np.abs(sub(obs.zeta)), z, rel, bins
        )

    if obs.rejected is not None:
        rej = sub(obs.rejected).astype(bool)
        res.rejected = {
            "rejected": int(rej.sum()),
            "rejected, |z| > 5": int(np.sum(rej & (np.abs(z) > 5))),
            "kept, |z| > 5": int(np.sum(~rej & (np.abs(z) > 5))),
        }

    order = np.argsort(-np.abs(z))[:worst]
    for j in order:
        r = rows[j]
        entry = {
            "hkl": tuple(int(v) for v in hkl[r]),
            "I": float(i_obs[j]),
            "kFc2": float(model[j]),
            "z": float(z[j]),
            "reflection's mean z": float(zbar[j]),
            "observations of it": int(count[inverse[j]]),
        }
        for name, v in (
            ("d", obs.d),
            ("image", obs.image),
            ("x", obs.x),
            ("y", obs.y),
            ("experiment", obs.experiment),
            ("partiality", obs.partiality),
        ):
            if v is not None:
                entry[name] = float(v[r])
        if obs.experiment is not None and obs.names:
            entry["images"] = obs.names.get(int(obs.experiment[r]), "")
        if obs.rejected is not None:
            entry["rejected by scaling"] = bool(obs.rejected[r])
        res.worst.append(entry)
    return res


# --- reporting ---------------------------------------------------------------


def report(res: Result, fcf: Fcf, extinction: float | None = None) -> str:
    out = []
    out.append(
        f"{res.matched} observations matched to {fcf.name}'s Fc^2, "
        f"{res.unmatched} with no reflection in the .fcf"
    )
    out.append(
        f"Scale k {res.k:.6g}"
        + (f", B {res.b:.3f} A^2" if res.b else "")
        + f"; R(F^2) over the observations {res.r_unmerged:.4f}"
    )
    if extinction is not None:
        out.append(
            f"EXTI {extinction:g} in the .res: the strong low-angle reflections'"
            " Fc^2 may or may not carry it -- read their bins with that in mind"
        )
    out.append(
        f"z = (I - k Fc^2) / sigma: variance {res.var_z:.3f}, robust "
        f"{res.robust_var_z:.3f}; about each reflection's own mean (the model's "
        f"share cancelled, the processing's sigmas alone) {res.robust_var_within:.3f} "
        "-- 1 where the sigmas are right"
    )
    for t, (obs, gauss) in res.tails.items():
        out.append(
            f"  |z| > {t:2d}: {obs:.5f} of the observations, Gaussian {gauss:.2g}"
        )
    out.append(
        f"Of sum z^2 over reflections observed more than once: {res.between:.1%} what"
        f" a reflection's observations share (the model's part, or a bias common to"
        f" them), {res.within:.1%} what each does about its own reflection's mean"
        f" (the processing's part)"
    )
    if res.rejected:
        out.append(
            "Scaling's outlier rejection: "
            + ", ".join(f"{k} {v}" for k, v in res.rejected.items())
        )
    for name, rows in res.bins.items():
        if not rows:
            continue
        out.append("")
        out.append(f"By {name}")
        out.append(
            f"  {'':>24} {'n':>7} {'mean z':>8} {'var z':>8} {'robust':>8}"
            f" {'rel. diff':>10} {'|z|>5':>7}"
        )
        for b in rows:
            out.append(
                f"  {b.label:>24} {b.n:7d} {b.mean_z:8.3f} {b.var_z:8.3f}"
                f" {b.robust_var_z:8.3f} {b.median_rel:10.4f} {b.over_5:7.4f}"
            )
    if res.worst:
        out.append("")
        out.append("The worst observations")
        for e in res.worst:
            out.append(
                "  "
                + ", ".join(
                    f"{k} {v:.4g}" if isinstance(v, float) else f"{k} {v}"
                    for k, v in e.items()
                )
            )
    return "\n".join(out)


def write_csv(path: str, obs: Observations, res: Result, reindex=None) -> None:
    import csv

    hkl = obs.hkl if reindex is None else obs.hkl @ reindex.T
    names = ["h", "k", "l", "I", "sigma", "kFc2", "z", "reflection_mean_z"]
    extra = [
        (n, v)
        for n, v in (
            ("d", obs.d),
            ("image", obs.image),
            ("x", obs.x),
            ("y", obs.y),
            ("experiment", obs.experiment),
            ("partiality", obs.partiality),
            ("background", obs.background),
            ("prf", obs.prf),
            ("sum", obs.sum),
            ("inverse_scale", obs.scale),
            ("rejected", obs.rejected),
        )
        if v is not None
    ]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(names + [n for n, _ in extra])
        for j, r in enumerate(res.rows):
            w.writerow(
                [int(v) for v in hkl[r]]
                + [
                    obs.intensity[r],
                    math.sqrt(obs.variance[r]),
                    res.model[j],
                    res.z[j],
                    res.zbar[j],
                ]
                + [v[r] for _, v in extra]
            )


# --- figures -----------------------------------------------------------------


def _running(x, y, bins=50):
    """Equal-count bins of x: each bin's centre, median y, and 16th and 84th
    percentiles of y."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < bins * 5:
        bins = max(5, len(x) // 5)
    order = np.argsort(x)
    out = []
    for part in np.array_split(order, bins):
        if len(part):
            out.append(
                (
                    np.median(x[part]),
                    np.median(y[part]),
                    *np.percentile(y[part], [16, 84]),
                )
            )
    return np.array(out)


def plot(obs: Observations, res: Result, outdir: str, reindex=None) -> list:
    """Figures of the observations against Fc^2, as PNG files in `outdir`;
    their paths returned. z is drawn clipped to +/-15; points beyond 5 sigma
    are drawn larger, and those scaling rejected as crosses."""
    import os

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    rows = res.rows
    z = res.z
    zc = np.clip(z, -15, 15)
    bad = np.abs(z) > 5
    rej = obs.rejected[rows] if obs.rejected is not None else np.zeros(len(rows), bool)
    exp = obs.experiment[rows] if obs.experiment is not None else np.zeros(len(rows))
    sets = [int(x) for x in np.unique(exp)]

    def name(x):
        return obs.names.get(x, f"experiment {x}") if obs.names else f"experiment {x}"

    paths = []

    def save(fig, stem):
        path = os.path.join(outdir, stem + ".png")
        fig.savefig(path, dpi=130, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)

    def panels(n):
        cols = 2 if n > 1 else 1
        r = (n + cols - 1) // cols
        fig, axes = plt.subplots(r, cols, figsize=(6.5 * cols, 3.6 * r), squeeze=False)
        return fig, axes.ravel()

    def scatter(ax, x, y, sel):
        good = sel & ~bad
        ax.scatter(x[good], y[good], s=2, c="0.6", lw=0, rasterized=True)
        hi = sel & bad & ~rej
        ax.scatter(x[hi], y[hi], s=9, c="tab:red", lw=0, label="|z| > 5, kept")
        rr = sel & rej
        ax.scatter(
            x[rr],
            y[rr],
            s=14,
            c="tab:blue",
            marker="x",
            lw=0.8,
            label="rejected by scaling",
        )

    # 1. z against image, each data set.
    if obs.image is not None:
        img = obs.image[rows]
        fig, axes = panels(len(sets))
        for ax, x in zip(axes, sets):
            sel = exp == x
            scatter(ax, img, zc, sel)
            run = _running(img[sel], z[sel], 40)
            if len(run):
                ax.plot(run[:, 0], run[:, 1], c="k", lw=1.2, label="running median")
                ax.fill_between(
                    run[:, 0], run[:, 2], run[:, 3], color="k", alpha=0.12, lw=0
                )
            ax.axhline(0, c="k", lw=0.5)
            ax.set_title(
                f"{name(x)}: {np.mean(bad[sel]):.1%} beyond 5 sigma", fontsize=10
            )
            ax.set_xlabel("image")
            ax.set_ylabel("z = (I - k Fc^2) / sigma")
            ax.set_ylim(-15.5, 15.5)
        axes[0].legend(fontsize=7, loc="lower left", markerscale=2)
        for ax in axes[len(sets) :]:
            ax.axis("off")
        fig.suptitle("Each observation against Fc^2, along the rotation", fontsize=11)
        save(fig, "z_by_image")

    # 2. z on the detector, each data set.
    if obs.x is not None and obs.y is not None:
        x_, y_ = obs.x[rows], obs.y[rows]
        fig, axes = panels(len(sets))
        for ax, x in zip(axes, sets):
            sel = exp == x
            order = np.argsort(np.abs(zc[sel]))  # the worst drawn last, on top
            sc = ax.scatter(
                x_[sel][order],
                y_[sel][order],
                c=zc[sel][order],
                s=np.where(bad[sel][order], 14, 3),
                cmap="RdBu_r",
                vmin=-10,
                vmax=10,
                lw=0,
                rasterized=True,
            )
            ax.set_title(name(x), fontsize=10)
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set_xlabel("x (px)")
            ax.set_ylabel("y (px)")
        for ax in axes[len(sets) :]:
            ax.axis("off")
        fig.colorbar(sc, ax=axes.tolist(), shrink=0.7, label="z, clipped to +/-10")
        fig.suptitle("Where on the detector: z, the worst drawn larger", fontsize=11)
        save(fig, "z_on_the_detector")

    # 3. Observed against calculated.
    i_obs = obs.intensity[rows]
    fig, ax = plt.subplots(figsize=(6.5, 6))
    pos = (res.model > 0) & (i_obs > 0)
    for x in sets:
        sel = pos & (exp == x) & ~bad
        ax.scatter(
            res.model[sel],
            i_obs[sel],
            s=2,
            lw=0,
            alpha=0.5,
            label=name(x),
            rasterized=True,
        )
    sel = pos & bad
    ax.scatter(res.model[sel], i_obs[sel], s=10, c="tab:red", lw=0, label="|z| > 5")
    lo, hi = np.percentile(res.model[pos], [0.1, 100])
    ax.plot([lo, hi], [lo, hi], c="k", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("k Fc^2")
    ax.set_ylabel("I, scaled")
    ax.legend(fontsize=7, markerscale=4)
    ax.set_title("Observed against calculated (the line: equal)", fontsize=10)
    save(fig, "observed_against_calculated")

    # 4. z against resolution and against intensity.
    fig, axes = plt.subplots(1, 2, figsize=(13, 4))
    if obs.d is not None:
        d = obs.d[rows]
        scatter(axes[0], d, zc, np.ones(len(rows), bool))
        run = _running(d, z)
        axes[0].plot(run[:, 0], run[:, 1], c="k", lw=1.2)
        axes[0].fill_between(
            run[:, 0], run[:, 2], run[:, 3], color="k", alpha=0.12, lw=0
        )
        axes[0].set_xscale("log")
        axes[0].invert_xaxis()
        axes[0].set_xlabel("d (A)")
    scatter(axes[1], np.where(res.model > 0, res.model, np.nan), zc, res.model > 0)
    run = _running(np.where(res.model > 0, res.model, np.nan), z)
    axes[1].plot(run[:, 0], run[:, 1], c="k", lw=1.2)
    axes[1].fill_between(run[:, 0], run[:, 2], run[:, 3], color="k", alpha=0.12, lw=0)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("k Fc^2")
    for ax in axes:
        ax.axhline(0, c="k", lw=0.5)
        ax.set_ylabel("z")
        ax.set_ylim(-15.5, 15.5)
    axes[1].legend(fontsize=7, loc="lower left", markerscale=2)
    fig.suptitle(
        "z against resolution and intensity: the running median, and 16-84 per cent",
        fontsize=11,
    )
    save(fig, "z_by_resolution_and_intensity")

    # 5. Normal probability plot.
    from scipy.stats import norm

    fig, ax = plt.subplots(figsize=(6, 6))
    zs = np.sort(z)
    expected = norm.ppf((np.arange(len(zs)) + 0.5) / len(zs))
    ax.scatter(expected, zs, s=3, lw=0, c="k", rasterized=True)
    ax.plot([-4.5, 4.5], [-4.5, 4.5], c="tab:red", lw=0.8, label="sigmas right")
    sr = math.sqrt(res.robust_var_z)
    ax.plot(
        [-4.5, 4.5],
        [-4.5 * sr, 4.5 * sr],
        c="tab:blue",
        lw=0.8,
        ls="--",
        label=f"slope {sr:.2f}, the robust spread",
    )
    ax.set_xlabel("expected, Gaussian")
    ax.set_ylabel("z, observed")
    ax.set_ylim(-25, 25)
    ax.legend(fontsize=8)
    ax.set_title("Normal probability plot of z", fontsize=10)
    save(fig, "normal_probability")

    # 6. The model's part: each reflection's mean z against resolution.
    if obs.d is not None:
        fig, ax = plt.subplots(figsize=(7, 4))
        d = obs.d[rows]
        _, first = np.unique(res.zbar + 1e-9 * res.model, return_index=True)
        ax.scatter(
            d[first],
            np.clip(res.zbar[first], -15, 15),
            s=6,
            c=np.log10(np.clip(res.model[first], 1, None)),
            cmap="viridis",
            lw=0,
        )
        ax.axhline(0, c="k", lw=0.5)
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_xlabel("d (A)")
        ax.set_ylabel("the reflection's mean z")
        ax.set_title(
            "What a reflection's observations share (colour: log10 k Fc^2)", fontsize=10
        )
        save(fig, "reflection_mean_z")
    # 7. Where the measurement departs: (Fo - Fc) / Fc on the detector, the
    # reflections strong enough for it to mean something -- z can be many
    # sigma for a small difference, and points then look "too strong"
    # (Graeme, 21 October).
    if obs.x is not None and obs.y is not None:
        x_, y_ = obs.x[rows], obs.y[rows]
        sig = np.sqrt(obs.variance[rows])
        strong = (res.model > 0) & (res.model > 10 * sig)
        fo = np.sqrt(np.clip(i_obs, 0, None))
        fc = np.sqrt(np.clip(res.model, 1e-12, None))
        dF = (fo - fc) / fc
        fig, axes = panels(len(sets))
        for ax, x in zip(axes, sets):
            sel = (exp == x) & strong
            order = np.argsort(np.abs(dF[sel]))
            sc = ax.scatter(
                x_[sel][order],
                y_[sel][order],
                c=np.clip(dF[sel][order], -0.3, 0.3),
                s=np.where(np.abs(dF[sel][order]) > 0.1, 12, 3),
                cmap="RdBu_r",
                vmin=-0.3,
                vmax=0.3,
                lw=0,
                rasterized=True,
            )
            lowf = np.mean(dF[sel] < -0.1)
            highf = np.mean(dF[sel] > 0.1)
            ax.set_title(
                f"{name(x)}: {lowf:.1%} over 10% low, {highf:.1%} over 10% high",
                fontsize=9,
            )
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set_xlabel("x (px)")
            ax.set_ylabel("y (px)")
        for ax in axes[len(sets) :]:
            ax.axis("off")
        fig.colorbar(
            sc, ax=axes.tolist(), shrink=0.7, label="(Fo - Fc) / Fc, clipped to +/-0.3"
        )
        fig.suptitle(
            "Where Fo departs from Fc: reflections with k Fc^2 above 10 sigma",
            fontsize=11,
        )
        save(fig, "relative_difference_on_the_detector")

    # 8. Against nearness to the rotation axis.
    if obs.zeta is not None:
        zt = np.abs(obs.zeta[rows])
        sig = np.sqrt(obs.variance[rows])
        strong = (res.model > 0) & (res.model > 10 * sig)
        fo = np.sqrt(np.clip(i_obs, 0, None))
        fc = np.sqrt(np.clip(res.model, 1e-12, None))
        dF = (fo - fc) / fc
        fig, axes = plt.subplots(1, 2, figsize=(13, 4))
        scatter(axes[0], zt, zc, np.ones(len(rows), bool))
        axes[0].set_ylabel("z")
        axes[0].set_ylim(-15.5, 15.5)
        axes[1].scatter(
            zt[strong],
            np.clip(dF[strong], -0.6, 0.6),
            s=3,
            c="0.4",
            lw=0,
            rasterized=True,
        )
        run = _running(zt[strong], dF[strong], 30)
        axes[1].plot(run[:, 0], run[:, 1], c="k", lw=1.2)
        axes[1].fill_between(
            run[:, 0], run[:, 2], run[:, 3], color="k", alpha=0.12, lw=0
        )
        axes[1].set_ylabel("(Fo - Fc) / Fc, k Fc^2 above 10 sigma")
        for ax in axes:
            ax.axhline(0, c="k", lw=0.5)
            ax.set_xscale("log")
            ax.set_xlabel("|zeta| (0 on the rotation axis)")
        fig.suptitle("Near the rotation axis", fontsize=11)
        save(fig, "by_zeta")
    return paths


# --- the outliers' pixels -----------------------------------------------------


def export_shoeboxes(
    obs: Observations,
    res: Result,
    integrated: str,
    out: str,
    worst: int = 300,
    normal: int = 100,
    seed: int = 0,
) -> dict:
    """The shoeboxes of the `worst` observations and of `normal` strong ones
    agreeing with Fc^2 within 1 sigma, for comparison: each scaled observation
    matched to its row of an integrated .rflx written with --save-shoeboxes --
    by experiment and position (within 1.5 px and half an image), the symmetry
    step having reindexed since -- and written to one compressed .npz: the
    boxes' pixels, background and mask concatenated with their offsets, and
    each one's figures. Small enough to send where the images are not."""
    from scipy.spatial import cKDTree

    from mxeq import rflx

    table = rflx.load_table(integrated)
    pos_name = "xyzobs.px.value" if "xyzobs.px.value" in table.columns else "xyzcal.px"
    pos = np.asarray(table[pos_name], float)
    ids = np.asarray(table["id"]) if "id" in table.columns else np.zeros(table.nrows)

    rows = res.rows
    z = res.z
    order = np.argsort(-np.abs(z))[:worst]
    sig = np.sqrt(obs.variance[rows])
    rng = np.random.default_rng(seed)
    calm = np.where((np.abs(z) < 1) & (res.model > 10 * sig))[0]
    picked_normal = (
        rng.choice(calm, size=min(normal, len(calm)), replace=False)
        if len(calm)
        else []
    )
    chosen = [(j, "worst") for j in order] + [(int(j), "normal") for j in picked_normal]

    trees = {}
    for e in np.unique(ids):
        sel = np.where(ids == e)[0]
        trees[e] = (
            cKDTree(pos[sel] * [1, 1, 3]),
            sel,
        )  # half an image weighs as 1.5 px
    want, kept = [], []
    for j, group in chosen:
        r = rows[j]
        e = obs.experiment[r] if obs.experiment is not None else 0
        if e not in trees:
            continue
        tree, sel = trees[e]
        d, k = tree.query([obs.x[r], obs.y[r], obs.image[r] * 3])
        if d <= 1.5 * math.sqrt(2) + 1.5:
            want.append(int(sel[k]))
            kept.append((j, group))
    boxes = rflx.load_shoeboxes(integrated, want)
    have = [i for i, b in enumerate(boxes["bbox"]) if b is not None]
    if not have:
        raise ValueError(
            f"{integrated}: no shoeboxes for the observations -- integrated with "
            "--save-shoeboxes?"
        )
    sizes = np.array([boxes["data"][i].size for i in have])
    offsets = np.concatenate([[0], np.cumsum(sizes)])

    def meta(v, dtype=float):
        return np.array([v[rows[kept[i][0]]] for i in have], dtype=dtype)

    payload = dict(
        group=np.array([kept[i][1] for i in have]),
        integrated_row=np.array([want[i] for i in have]),
        bbox=np.array([boxes["bbox"][i] for i in have]),
        offsets=offsets,
        data=np.concatenate([boxes["data"][i].ravel() for i in have]).astype(
            np.float32
        ),
        background=np.concatenate(
            [boxes["background"][i].ravel() for i in have]
        ).astype(np.float32),
        mask=np.concatenate([boxes["mask"][i].ravel() for i in have]).astype(np.uint8),
        hkl=np.array([obs.hkl[rows[kept[i][0]]] for i in have]),
        z=np.array([z[kept[i][0]] for i in have]),
        I=meta(obs.intensity),
        sigma=np.sqrt(meta(obs.variance)),
        kFc2=np.array([res.model[kept[i][0]] for i in have]),
    )
    for name in (
        "x",
        "y",
        "image",
        "experiment",
        "partiality",
        "d",
        "zeta",
        "background",
    ):
        v = getattr(obs, name)
        if v is not None:
            payload[name if name != "background" else "background_mean"] = meta(v)
    np.savez_compressed(out, **payload)
    return {
        "worst": sum(1 for i in have if kept[i][1] == "worst"),
        "normal": sum(1 for i in have if kept[i][1] == "normal"),
        "unmatched": len(chosen) - len(have),
    }
