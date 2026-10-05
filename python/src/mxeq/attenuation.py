"""Attenuation judged in scaling: each observation good or attenuated, from its
equivalents and its background together -- a prototype (docs/backstop.md).

Beside a backstop an observation can be attenuated -- its signal and its
background dimmed together by the shadow's edge -- and an attenuated
observation is low. The equivalents say whether an observation is low; its
background, against its equivalents' backgrounds, says why. And attenuation only
lowers: an observation above the reflection's value cannot be attenuated, and
the truth lies with the upper consistent group.

Each observation of a reflection is good, I ~ N(Ibar, sigma_good^2), or
attenuated, I ~ N(T Ibar, sigma^2) with T uniform on (0, 1):

* **the prior** that it is attenuated, from its background ratio -- its
  background over the 75th percentile of its reflection's, each first divided by
  the background model's expectation where the table has one, which takes out
  the rotation's and the polarisation's variation between equivalents -- a
  smooth step, high where the ratio is low. Local: beside the backstop the
  flare raises the group's backgrounds and can bias the reference, which cannot
  be helped -- modelling it is hard;
* **the likelihood**, one-sided: the good with sigma_good the observation's
  sigma at Ibar -- its variance scaled by Ibar over its own intensity, counting
  variance growing with intensity -- never less than its own; the attenuated
  its own normal averaged over T, nil above Ibar;
* **expectation-maximisation**: each observation's probability of attenuation,
  then Ibar as the mean weighted by the probability of being good, iterated.

An observation more likely attenuated than good is flagged; rejected, not
corrected. A reflection too weak to judge by its equivalents (Ibar / sigma below
--min-i-sigma) -- perhaps weak because attenuated -- or observed fewer than three
times is judged by its background alone: the prior, its background against the
background model's expectation where the table has one, and kept where it has
none. --apply-to writes mxi_scale's next input with the flags applied, so
that scaling and judging can be iterated until the flags settle.

    mxeq attenuation scaled.expt scaled.refl --apply-to symmetrized.refl -o next.refl
    mxeq attenuation scaled.expt scaled.refl --show 3,1,1
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import ndtr

from . import equivalents as eq
from . import expt as expt_module
from . import refl

INTEGRATED_SUM = 1 << 8
INTEGRATED_PRF = 1 << 9


@dataclass
class Judged:
    rows: np.ndarray  # into the scaled table
    keys: np.ndarray  # each observation's reflection
    intensity: np.ndarray
    sigma: np.ndarray
    background: np.ndarray
    ratio: np.ndarray  # background over the reflection's 75th percentile
    prior: np.ndarray
    probability: np.ndarray  # of being attenuated
    ibar: np.ndarray  # the reflection's, per observation
    judged: np.ndarray  # the reflection strong enough to judge
    flagged: np.ndarray
    d: np.ndarray
    iterations: int


def _group_percentile(inverse, values, q, groups):
    """Each group's q-th percentile of values, per observation."""
    order = np.lexsort((values, inverse))
    counts = np.bincount(inverse, minlength=groups)
    starts = np.r_[0, np.cumsum(counts)[:-1]]
    pick = starts + np.floor(q * np.maximum(counts - 1, 0)).astype(np.int64)
    per_group = values[order][np.minimum(pick, len(values) - 1)]
    return per_group[inverse]


def judge(
    keys,
    intensity,
    variance,
    background,
    expected,
    d,
    prior_mid=0.6,
    prior_width=0.1,
    min_i_sigma=3.0,
    iterations=30,
    threshold=0.5,
) -> dict:
    unique, inverse = np.unique(keys, return_inverse=True)
    groups = len(unique)
    sigma = np.sqrt(np.maximum(variance, 1e-30))
    has_model = np.isfinite(expected) & (expected > 0)
    norm = np.where(
        has_model, background / np.where(has_model, expected, 1.0), background
    )
    reference = _group_percentile(inverse, norm, 0.75, groups)
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(reference > 0, norm / reference, 1.0)
    ratio = np.clip(np.nan_to_num(ratio, nan=1.0), 0.0, 10.0)
    # A reflection observed fewer than three times has no group to compare a
    # background with: the background model's expectation, where the table has
    # it, is the reference -- the global z the design keeps for this -- and
    # without it nothing can be said.
    count = np.bincount(inverse, minlength=groups)[inverse]
    few = count < 3
    ratio = np.where(few, np.where(has_model, np.clip(norm, 0.0, 10.0), 1.0), ratio)
    informed = ~few | has_model
    prior = 1.0 / (1.0 + np.exp((ratio - prior_mid) / prior_width))
    prior = np.clip(prior, 0.01, 0.99)

    # Start from the observations whose backgrounds are normal: near the upper group.
    w0 = (1.0 - prior) / variance
    ibar = np.bincount(inverse, weights=w0 * intensity, minlength=groups) / np.maximum(
        np.bincount(inverse, weights=w0, minlength=groups), 1e-300
    )
    p = prior.copy()
    for it in range(iterations):
        ib = ibar[inverse]
        scale = np.maximum(ib, 0.0) / np.maximum(intensity, sigma)
        var_good = variance * np.maximum(1.0, scale)
        good = np.exp(-0.5 * (intensity - ib) ** 2 / var_good) / np.sqrt(
            2 * np.pi * var_good
        )
        with np.errstate(invalid="ignore", divide="ignore"):
            att = np.where(
                ib > 0,
                (ndtr((ib - intensity) / sigma) - ndtr(-intensity / sigma))
                / np.maximum(ib, 1e-300),
                0.0,
            )
        num = prior * att
        den = num + (1.0 - prior) * good
        new_p = np.where(den > 0, num / np.maximum(den, 1e-300), prior)
        w = (1.0 - new_p) / var_good
        new_ibar = np.bincount(
            inverse, weights=w * intensity, minlength=groups
        ) / np.maximum(np.bincount(inverse, weights=w, minlength=groups), 1e-300)
        done = np.allclose(new_ibar, ibar, rtol=1e-6, atol=1e-9) and np.allclose(
            new_p, p, atol=1e-6
        )
        ibar, p = new_ibar, new_p
        if done:
            break
    sigma_ibar = 1.0 / np.sqrt(
        np.maximum(
            np.bincount(inverse, weights=(1.0 - p) / variance, minlength=groups), 1e-300
        )
    )
    # Judged by the equivalents where the reflection is strong enough and has
    # them; otherwise -- weak, perhaps weak because attenuated, or observed
    # fewer than three times -- by its background alone, the prior.
    judged = (ibar / sigma_ibar >= min_i_sigma)[inverse] & ~few
    probability = np.where(judged, p, np.where(informed, prior, 0.0))
    flagged = probability > threshold
    return dict(
        ratio=ratio,
        prior=prior,
        probability=probability,
        ibar=ibar[inverse],
        judged=judged,
        flagged=flagged,
        iterations=it + 1,
    )


def load_and_judge(expt_path: str, refl_path: str, **options) -> Judged:
    e = expt_module.load(expt_path)
    if not e.experiments or e[0].crystal is None or not e[0].crystal.hall:
        raise ValueError(f"{expt_path} has no crystal with a space group")
    t = refl.load(refl_path)
    c = t.columns
    for name in (
        "intensity.scale.value",
        "intensity.scale.variance",
        "inverse_scale_factor",
        "background.mean",
        "flags",
    ):
        if name not in c:
            raise ValueError(
                f"{refl_path} has no {name}: a scaled table of mxi_integrate's is wanted"
            )
    flags = np.asarray(c["flags"]).astype(np.int64)
    # Every observation scaling considered -- those its own outlier test
    # rejected among them, as unharmed observations can be.
    considered = (flags & (eq.SCALED | eq.OUTLIER_IN_SCALING)) != 0
    isf = np.asarray(c["inverse_scale_factor"], float)
    var = np.asarray(c["intensity.scale.variance"], float)
    keep = considered & (isf > 0) & (var > 0)
    rows = np.nonzero(keep)[0]
    intensity = np.asarray(c["intensity.scale.value"], float)[keep] / isf[keep]
    variance = var[keep] / isf[keep] ** 2
    background = np.asarray(c["background.mean"], float)[keep]
    expected = (
        np.asarray(c["background.expected"], float)[keep]
        if "background.expected" in c
        else np.full(len(rows), np.nan)
    )
    d = np.asarray(c["d"], float)[keep] if "d" in c else np.full(len(rows), np.nan)
    keys = eq.asu_keys(
        np.asarray(c["miller_index"]).reshape(-1, 3)[keep], e[0].crystal.hall
    )
    r = judge(keys, intensity, variance, background, expected, d, **options)
    return Judged(
        rows,
        keys,
        intensity,
        np.sqrt(variance),
        background,
        r["ratio"],
        r["prior"],
        r["probability"],
        r["ibar"],
        r["judged"],
        r["flagged"],
        d,
        r["iterations"],
    )


def report(j: Judged, shells: int = 10) -> str:
    groups = len(np.unique(j.keys))
    judged_groups = len(np.unique(j.keys[j.judged]))
    hit_groups = len(np.unique(j.keys[j.flagged]))
    lines = [
        f"{len(j.rows)} observations of {groups} reflections; {judged_groups} reflections strong enough to "
        f"judge; {j.iterations} iterations",
        f"{int(j.flagged.sum())} observations more likely attenuated than good, in {hit_groups} reflections",
        "",
        f"  {'d (A)':>15} {'obs':>9} {'flagged':>8} {'%':>6} {'ratio, flagged':>15} {'ratio, kept':>12}",
    ]
    ok = np.isfinite(j.d) & (j.d > 0)
    s = 1.0 / np.where(ok, j.d, np.nan)
    edges = np.nanquantile(s, np.linspace(0, 1, shells + 1))
    for k in range(shells):
        sel = (
            ok
            & (s >= edges[k])
            & ((s <= edges[k + 1]) if k == shells - 1 else (s < edges[k + 1]))
        )
        if not sel.any():
            continue
        f = sel & j.flagged
        kept = sel & ~j.flagged
        rf = f"{np.median(j.ratio[f]):15.2f}" if f.any() else f"{'':>15}"
        lines.append(
            f"  {j.d[sel].max():6.2f} -{j.d[sel].min():7.2f} {sel.sum():9d} {f.sum():8d} "
            f"{100 * f.sum() / sel.sum():6.3f} {rf} {np.median(j.ratio[kept]):12.2f}"
        )
    lines += [
        "",
        "  ratio is each observation's background over its reflection's 75th percentile: an",
        "  attenuated observation's is low, signal and background dimmed together.",
    ]
    return "\n".join(lines)


def show(j: Judged, table: refl.ReflectionTable, hall: str, index) -> str:
    k = eq.asu_keys(np.array([index]), hall)[0]
    sel = np.nonzero(j.keys == k)[0]
    hkl = np.asarray(table.columns["miller_index"]).reshape(-1, 3)
    lines = [
        (
            f"== {list(index)}: {len(sel)} observations; Ibar {j.ibar[sel[0]]:.2f}"
            if len(sel)
            else f"== {list(index)}: no observations scaling considered"
        ),
        f"  {'hkl':>13} {'I':>10} {'sigma':>8} {'bg':>7} {'ratio':>6} {'prior':>6} {'P(att)':>7}  decision",
    ]
    for i in sel[np.argsort(-j.intensity[sel])]:
        lines.append(
            f"  {str(hkl[j.rows[i]].tolist()):>13} {j.intensity[i]:10.2f} {j.sigma[i]:8.2f} "
            f"{j.background[i]:7.2f} {j.ratio[i]:6.2f} {j.prior[i]:6.2f} {j.probability[i]:7.3f}  "
            f"{'attenuated' if j.flagged[i] else 'kept'}"
            f"{'' if j.judged[i] else ', by its background alone'}"
        )
    return "\n".join(lines)


def annotate(table: refl.ReflectionTable, j: Judged) -> refl.ReflectionTable:
    n = table.nrows
    prob = np.zeros(n)
    prob[j.rows] = j.probability
    ratio = np.ones(n)
    ratio[j.rows] = j.ratio
    table.columns["attenuation.probability"] = prob
    table.types["attenuation.probability"] = "double"
    table.columns["attenuation.background_ratio"] = ratio
    table.types["attenuation.background_ratio"] = "double"
    return table


def apply_to(
    target_path: str, scaled: refl.ReflectionTable, j: Judged
) -> tuple[refl.ReflectionTable, int]:
    """mxi_scale's input with the flagged observations' integrated flags cleared,
    so the next round leaves them out. The rows are the scaled table's: mxi_scale
    writes every row of its input, in order -- checked by the Miller indices."""
    t = refl.load(target_path)
    if t.nrows != scaled.nrows or not np.array_equal(
        np.asarray(t.columns["miller_index"]).reshape(-1, 3),
        np.asarray(scaled.columns["miller_index"]).reshape(-1, 3),
    ):
        raise ValueError(
            f"{target_path}'s rows are not the scaled table's: give the table mxi_scale scaled"
        )
    flags = np.asarray(t.columns["flags"]).astype(np.int64).copy()
    rows = j.rows[j.flagged]
    flags[rows] &= ~(INTEGRATED_SUM | INTEGRATED_PRF)
    flags[rows] |= eq.EXCLUDED_FOR_SCALING
    t.columns["flags"] = flags.astype(np.asarray(t.columns["flags"]).dtype)
    return t, len(rows)
