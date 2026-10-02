"""Where profile fitting fails, from an integrated table's profile.failure.

mxi_integrate writes each reflection's reason for not being profile fitted as
profile.failure -- 0 fitted, else a code -- since scaling takes profile-fitted
observations only and a data set that loses a third of them says nothing, by
the count alone, of why. This gives the count of each reason and then, for the
failures, the fraction failing binned by each thing that could be the cause:
position in the scan, |zeta|, the box's depth in images and width in pixels,
place on the detector, resolution and summed I/sigma. A cause shows as the one
variable whose bins differ sharply where the others are flat.

    mxeq failures integrated.refl
    mxeq failures integrated.refl --code 4     # one reason only
"""

from __future__ import annotations

import numpy as np

from . import refl

REASONS = {
    1: "the box rejected before fitting",
    2: "no reference profile for its place and scan block",
    3: "the profile not carried onto its pixels",
    4: "the profile zero over the whole foreground: off its box",
    5: "none of the foreground measured: wholly in a gap, say",
    6: "under --least-measured of the profile measured",
    7: "the least squares had no solution",
    8: "the device's fit failed (without --gpu says why)",
}


def _table(name, values, failing, edges=None, labels=None, bins=10):
    lines = [f"  {name}", f"    {'bin':<22} {'n':>9} {'failing':>9} {'%':>7}"]
    v = np.asarray(values, float)
    ok = np.isfinite(v)
    if edges is None:
        inside = v[ok]
        if inside.size == 0:
            return lines + ["    (no values)"]
        edges = np.unique(np.quantile(inside, np.linspace(0, 1, bins + 1)))
        if len(edges) < 2:
            edges = np.array([inside.min(), inside.max() + 1])
        edges[-1] = np.nextafter(edges[-1], np.inf)
    for k in range(len(edges) - 1):
        sel = ok & (v >= edges[k]) & (v < edges[k + 1])
        n = int(sel.sum())
        f = int((sel & failing).sum())
        label = labels[k] if labels else f"[{edges[k]:.4g}, {edges[k + 1]:.4g})"
        pct = f"{100.0 * f / n:7.1f}" if n else f"{'':>7}"
        lines.append(f"    {label:<22} {n:>9} {f:>9} {pct}")
    return lines


def analyse(path: str, code: int | None = None) -> str:
    t = refl.load(path)
    c = t.columns
    if "profile.failure" not in c:
        raise ValueError(
            f"{path} has no profile.failure: integrate with a newer mxi_integrate"
        )
    why = np.asarray(c["profile.failure"]).ravel().astype(int)
    lines = [f"{path}: {t.nrows} reflections"]
    codes, counts = np.unique(why, return_counts=True)
    for k, n in zip(codes.tolist(), counts.tolist()):
        reason = "profile fitted" if k == 0 else REASONS.get(k, f"code {k}")
        lines.append(f"  {n:>10}  {k}  {reason}")
    failing = why == code if code is not None else why != 0
    what = (
        REASONS.get(code, f"code {code}")
        if code is not None
        else "not profile fitted, any reason"
    )
    lines += ["", f"Failing ({what}), by:"]

    px = np.asarray(c["xyzcal.px"]).reshape(-1, 3) if "xyzcal.px" in c else None
    if px is not None:
        lines += _table("position in the scan (image)", px[:, 2], failing, bins=12)
    if "zeta" in c:
        lines += _table("|zeta|", np.abs(np.asarray(c["zeta"]).ravel()), failing)
    if "bbox" in c:
        bbox = np.asarray(c["bbox"]).reshape(-1, 6)
        depth = bbox[:, 5] - bbox[:, 4]
        top = int(min(depth.max(), 30)) if depth.size else 1
        edges = np.arange(int(depth.min()), top + 2) if depth.size else None
        lines += _table(
            "the box's depth in images",
            depth,
            failing,
            edges=edges,
            labels=[str(k) for k in edges[:-1]] if edges is not None else None,
        )
        lines += _table("the box's width in pixels", bbox[:, 1] - bbox[:, 0], failing)
    if px is not None:
        x, y = px[:, 0], px[:, 1]
        cols = np.digitize(x, np.quantile(x, [1 / 3, 2 / 3]))
        rows = np.digitize(y, np.quantile(y, [1 / 3, 2 / 3]))
        cell = rows * 3 + cols
        names = [
            f"{r} {s}"
            for r in ("top", "middle", "bottom")
            for s in ("left", "centre", "right")
        ]
        lines += _table(
            "place on the detector",
            cell,
            failing,
            edges=np.arange(10) - 0.5,
            labels=names,
        )
    if "d" in c:
        lines += _table("resolution d (A)", np.asarray(c["d"]).ravel(), failing)
    if "intensity.sum.value" in c and "intensity.sum.variance" in c:
        v = np.asarray(c["intensity.sum.variance"]).ravel()
        isig = np.where(
            v > 0,
            np.asarray(c["intensity.sum.value"]).ravel()
            / np.sqrt(np.where(v > 0, v, 1)),
            np.nan,
        )
        lines += _table(
            "summed I/sigma",
            isig,
            failing,
            edges=[-np.inf, 0, 1, 3, 10, 30, np.inf],
            labels=["below 0", "0-1", "1-3", "3-10", "10-30", "30 and over"],
        )
    return "\n".join(lines)
