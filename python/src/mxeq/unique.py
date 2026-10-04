"""Two scaled data sets compared one unique reflection at a time.

Where two programs' merging statistics differ by a little -- DIALS's inner shell
a shade better in CC1/2 and I/sigma than mxi's -- the totals do not say why.
This matches the two data sets' unique reflections and splits the difference
into its parts.

A merged reflection's I/sigma is its mean intensity over its sigma, and that
sigma is the per-observation sigma over the root of the multiplicity. Written

    I/sigma = I_mean * sqrt(n) / sigma_obs,   sigma_obs = sigma_mean * sqrt(n)

the log of the ratio between the two data sets is exactly the sum of three:
the mean intensities' (a scale, or a bias), the multiplicities' (outlier
rejection, observations lost), and the per-observation sigmas' (variance
estimates, the error model). CC1/2 does not see sigma at all: it measures how
far observations scatter about their mean. So each reflection's scatter -- the
standard deviation of its observations, no sigma in it -- is compared too: less
scatter is intensities that agree better; the same scatter and a smaller sigma
is the error model.

Intensities are what merging uses: intensity.scale.value over
inverse_scale_factor, the variance likewise, observations scaling took
(flagged scaled, not outliers). The two are grouped by the first's space group,
and the second reindexed by whichever of the lattice's operators makes the two
data sets' merged intensities agree best -- the identity, when both programs
chose the same setting; the report says which.

    mxeq unique dials.expt dials.refl mxi.expt mxi.refl
    mxeq unique a.expt a.refl b.expt b.refl --shells 12 --csv matched.csv
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import equivalents as eq
from . import expt as expt_module
from . import refl
from . import reindex

#: DIALS's flag for a reflection with a pixel over the trusted range. mxi's
#: integrator does not set it: it does not apply the trusted range.
OVERLOADED = 1 << 10


@dataclass
class Observations:
    label: str
    hall: str
    hkl: np.ndarray
    intensity: np.ndarray
    variance: np.ndarray
    d: np.ndarray
    overloaded: np.ndarray | None = None  # per kept observation
    flagged_overloaded: int = 0  # rows the table flags, used or not


def load(expt_path: str, refl_path: str, label: str) -> Observations:
    e = expt_module.load(expt_path)
    if not e.experiments or e[0].crystal is None or not e[0].crystal.hall:
        raise ValueError(f"{expt_path} has no crystal with a space group")
    t = refl.load(refl_path)
    c = t.columns
    for name in (
        "intensity.scale.value",
        "intensity.scale.variance",
        "inverse_scale_factor",
    ):
        if name not in c:
            raise ValueError(f"{refl_path} has no {name}: is it scaled?")
    flags = np.asarray(c["flags"]).astype(np.int64)
    keep = (flags & (eq.OUTLIER_IN_SCALING | eq.EXCLUDED_FOR_SCALING)) == 0
    if np.any(flags & eq.SCALED):
        keep &= (flags & eq.SCALED) != 0
    isf = np.asarray(c["inverse_scale_factor"], float)
    value = np.asarray(c["intensity.scale.value"], float)
    var = np.asarray(c["intensity.scale.variance"], float)
    keep &= (isf > 0) & (var > 0) & np.isfinite(value)
    d = np.asarray(c["d"], float) if "d" in c else np.full(t.nrows, np.nan)
    return Observations(
        label,
        e[0].crystal.hall,
        np.asarray(c["miller_index"]).reshape(-1, 3)[keep],
        value[keep] / isf[keep],
        var[keep] / isf[keep] ** 2,
        d[keep],
        ((flags & OVERLOADED) != 0)[keep],
        int(((flags & OVERLOADED) != 0).sum()),
    )


@dataclass
class Merged:
    keys: np.ndarray  # one per unique reflection, sorted
    n: np.ndarray
    mean: np.ndarray  # weighted
    sigma: np.ndarray  # of the weighted mean
    sigma_obs: np.ndarray  # sigma * sqrt(n): the per-observation sigma
    chi2: np.ndarray  # per degree of freedom, NaN for one observation
    scatter: np.ndarray  # standard deviation of the observations, NaN for one
    half1: np.ndarray  # unweighted means of a random half each, NaN where empty
    half2: np.ndarray
    d: np.ndarray
    overloaded: np.ndarray  # observations flagged overloaded, per reflection


def merge(keys: np.ndarray, o: Observations, seed: int = 0) -> Merged:
    unique, inverse, n = np.unique(keys, return_inverse=True, return_counts=True)
    w = 1.0 / o.variance
    sw = np.bincount(inverse, weights=w)
    mean = np.bincount(inverse, weights=w * o.intensity) / sw
    sigma = 1.0 / np.sqrt(sw)
    resid = (o.intensity - mean[inverse]) ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        chi2 = np.where(
            n > 1, np.bincount(inverse, weights=resid * w) / (n - 1), np.nan
        )
        plain = np.bincount(inverse, weights=o.intensity) / n
        dev = (o.intensity - plain[inverse]) ** 2
        scatter = np.where(
            n > 1, np.sqrt(np.bincount(inverse, weights=dev) / (n - 1)), np.nan
        )
    # Random halves within each reflection, balanced: a random rank in each
    # group, even ranks one half and odd the other.
    rng = np.random.default_rng(seed)
    order = np.lexsort((rng.random(len(inverse)), inverse))
    rank = np.empty(len(inverse), dtype=np.int64)
    starts = np.r_[0, np.cumsum(n)[:-1]]
    rank[order] = np.arange(len(inverse)) - np.repeat(starts, n)
    halves = []
    for parity in (0, 1):
        sel = (rank % 2) == parity
        cnt = np.bincount(inverse[sel], minlength=len(unique))
        s = np.bincount(inverse[sel], weights=o.intensity[sel], minlength=len(unique))
        with np.errstate(invalid="ignore", divide="ignore"):
            halves.append(np.where(cnt > 0, s / np.maximum(cnt, 1), np.nan))
    d = np.bincount(inverse, weights=np.nan_to_num(o.d)) / n
    return Merged(
        unique,
        n,
        mean,
        sigma,
        sigma * np.sqrt(n),
        chi2,
        scatter,
        halves[0],
        halves[1],
        d,
        np.bincount(
            inverse,
            weights=(
                o.overloaded
                if o.overloaded is not None
                else np.zeros(len(inverse), bool)
            ).astype(float),
            minlength=len(unique),
        ).astype(int),
    )


def _cc(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan")
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def best_operator(a: Merged, b_obs: Observations, hall: str):
    """The lattice operator on the second data set's indices under which the
    two merged data sets agree best, its correlation, and the identity's."""
    results = []
    for m in reindex.candidate_operators():
        m = np.asarray(m)
        if abs(round(abs(np.linalg.det(m)))) != 1:
            continue
        hkl = np.rint(b_obs.hkl @ m).astype(np.int64)
        keys = eq.asu_keys(hkl, hall)
        unique, inverse = np.unique(keys, return_inverse=True)
        w = 1.0 / b_obs.variance
        mean = np.bincount(inverse, weights=w * b_obs.intensity) / np.bincount(
            inverse, weights=w
        )
        _, ia, ib = np.intersect1d(a.keys, unique, return_indices=True)
        if len(ia) < 10:
            continue
        results.append((_cc(a.mean[ia], mean[ib]), m))
    identity = next(
        (cc for cc, m in results if np.array_equal(m, np.eye(3, dtype=m.dtype))),
        float("nan"),
    )
    cc, m = max(
        results,
        key=lambda r: (np.nan_to_num(r[0], nan=-2.0), np.array_equal(r[1], np.eye(3))),
    )
    return m, cc, identity


def compare(
    a_obs: Observations, b_obs: Observations, shells: int = 10, csv: str | None = None
) -> str:
    hall = a_obs.hall
    lines = [
        f"{a_obs.label} against {b_obs.label}, grouped by {reindex.space_group_name(hall) or hall}"
    ]
    if b_obs.hall != hall:
        lines.append(
            f"  (their space groups differ: {b_obs.label}'s is "
            f"{reindex.space_group_name(b_obs.hall) or b_obs.hall}; grouped by the first's)"
        )
    a = merge(eq.asu_keys(a_obs.hkl, hall), a_obs)
    m, cc, identity = best_operator(a, b_obs, hall)
    is_identity = np.array_equal(m, np.eye(3, dtype=m.dtype))
    lines.append(
        f"  {b_obs.label}'s indices as they are"
        if is_identity
        else f"  {b_obs.label} reindexed by {m.astype(int).tolist()}: merged intensities CC "
        f"{cc:.4f}, against {identity:.4f} as they are"
    )
    b_hkl = np.rint(b_obs.hkl @ m).astype(np.int64)
    b = merge(eq.asu_keys(b_hkl, hall), b_obs)
    common, ia, ib = np.intersect1d(a.keys, b.keys, return_indices=True)
    lines.append(
        f"  unique reflections: {len(a.keys)} and {len(b.keys)}, {len(common)} in both; "
        f"observations {len(a_obs.intensity)} and {len(b_obs.intensity)}"
    )
    for o in (a_obs, b_obs):
        used = int(o.overloaded.sum()) if o.overloaded is not None else 0
        lines.append(
            f"  {o.label}: {o.flagged_overloaded} rows flagged overloaded, {used} of them used in scaling"
            + (
                ""
                if o.flagged_overloaded
                else " (none flagged; mxi_integrate never sets the flag)"
            )
        )
    lines.append("")

    d = a.d[ia]
    edges = np.quantile(d, np.linspace(0, 1, shells + 1))[::-1]
    shell_of = np.clip(np.searchsorted(-edges, -d, side="right") - 1, 0, shells - 1)

    def block(title, rows, header):
        lines.append(title)
        lines.append("    " + header)
        lines.extend("    " + r for r in rows)
        lines.append("")

    per = []
    for name, x, idx in ((a_obs.label, a, ia), (b_obs.label, b, ib)):
        rows = []
        for s in range(shells):
            sel = idx[shell_of == s]
            dd = d[shell_of == s]
            rows.append(
                f"{dd.min():6.2f}-{dd.max():6.2f} {len(sel):7d} {np.mean(x.n[sel]):7.2f} "
                f"{np.mean(x.mean[sel] / x.sigma[sel]):8.2f} {_cc(x.half1[sel], x.half2[sel]):7.4f} "
                f"{np.nanmedian(x.chi2[sel]):7.3f} {int(x.overloaded[sel].sum()):6d}"
            )
        sel = idx
        rows.append(
            f"{'overall':>13} {len(sel):7d} {np.mean(x.n[sel]):7.2f} {np.mean(x.mean[sel] / x.sigma[sel]):8.2f} "
            f"{_cc(x.half1[sel], x.half2[sel]):7.4f} {np.nanmedian(x.chi2[sel]):7.3f} "
            f"{int(x.overloaded[sel].sum()):6d}"
        )
        per.append(rows)
        block(
            f"  {name}, over the reflections both have",
            rows,
            f"{'d (A)':>13} {'unique':>7} {'mult':>7} {'I/sigma':>8} {'CC1/2':>7} {'chi2/nu':>7} {'ovld':>6}",
        )

    rows = []
    for s in list(range(shells)) + [None]:
        mask = np.ones(len(ia), bool) if s is None else shell_of == s
        sa, sb = ia[mask], ib[mask]
        pos = (a.mean[sa] > 0) & (b.mean[sb] > 0)
        lr_mean = np.log(a.mean[sa][pos] / b.mean[sb][pos])
        lr_n = 0.5 * np.log(a.n[sa][pos] / b.n[sb][pos])
        lr_sig = -np.log(a.sigma_obs[sa][pos] / b.sigma_obs[sb][pos])
        lr_total = lr_mean + lr_n + lr_sig
        both = (
            np.isfinite(a.scatter[sa])
            & np.isfinite(b.scatter[sb])
            & (b.scatter[sb] > 0)
        )
        scatter = (
            np.median(a.scatter[sa][both] / b.scatter[sb][both])
            if both.any()
            else np.nan
        )
        # The same, each scatter over its reflection's mean intensity: a
        # difference of overall scale between the two cancels, which it does
        # not in the scatter's ratio alone.
        rel = both & (a.mean[sa] > 0) & (b.mean[sb] > 0)
        relative = (
            np.median(
                (a.scatter[sa][rel] / a.mean[sa][rel])
                / (b.scatter[sb][rel] / b.mean[sb][rel])
            )
            if rel.any()
            else np.nan
        )
        label = "overall" if s is None else f"{d[mask].min():6.2f}-{d[mask].max():6.2f}"

        def pct(x):
            return f"{100 * (np.exp(np.median(x)) - 1):+7.2f}" if len(x) else f"{'':>7}"

        rows.append(
            f"{label:>13} {_cc(a.mean[sa], b.mean[sb]):7.4f} {pct(lr_total)} {pct(lr_mean)} "
            f"{pct(lr_n)} {pct(lr_sig)} {100 * (scatter - 1):+8.2f} {100 * (relative - 1):+9.2f}"
        )
    block(
        f"  {a_obs.label} against {b_obs.label}: medians, in per cent, of the ratio per reflection",
        rows,
        f"{'d (A)':>13} {'CC I':>7} {'I/sig':>7} {'= mean':>7} {'+ mult':>7} {'+ 1/sig':>7} "
        f"{'scatter':>8} {'scatter/I':>9}",
    )
    lines.append(
        "  I/sig's ratio is the product of the three after it -- the mean intensity's, the root"
    )
    lines.append(
        "  of the multiplicity's, and the inverse of the per-observation sigma's -- over"
    )
    lines.append(
        "  reflections positive in both; scatter is the observations' spread, no sigma in it,"
    )
    lines.append(
        "  which CC1/2 follows: less scatter, observations that agree better; scatter/I is it"
    )
    lines.append(
        "  over the reflection's mean intensity, so that a difference of overall scale cancels."
    )

    if csv:
        with open(csv, "w") as f:
            f.write(
                "h,k,l,d,n_a,n_b,i_a,i_b,sigma_a,sigma_b,sigma_obs_a,sigma_obs_b,chi2_a,chi2_b,"
                "scatter_a,scatter_b\n"
            )
            for j in range(len(common)):
                k = common[j]
                h = (
                    ((k >> 22) & 2047) - 1024,
                    ((k >> 11) & 2047) - 1024,
                    (k & 2047) - 1024,
                )
                p, q = ia[j], ib[j]
                f.write(
                    f"{h[0]},{h[1]},{h[2]},{a.d[p]:.4f},{a.n[p]},{b.n[q]},{a.mean[p]:.6g},{b.mean[q]:.6g},"
                    f"{a.sigma[p]:.6g},{b.sigma[q]:.6g},{a.sigma_obs[p]:.6g},{b.sigma_obs[q]:.6g},"
                    f"{a.chi2[p]:.6g},{b.chi2[q]:.6g},{a.scatter[p]:.6g},{b.scatter[q]:.6g}\n"
                )
        lines.append(f"\nWrote the {len(common)} matched reflections to {csv}")
    return "\n".join(lines)
