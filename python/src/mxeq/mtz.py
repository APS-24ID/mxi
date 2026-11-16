"""Two unmerged, scaled MTZ files compared: DIALS's and mxi's, say.

Graeme, 15 November: the main result of either pipeline is its scaled.mtz,
so compare those, side by side. Four parts:

1. The files: space group, cell, wavelength, rows, batches, resolution, and
   the columns each has -- differences here are the export's to explain.
2. The indexing: where the lattice has more symmetry than the space group
   (insulin's I 2 3 in a cubic lattice), the two may be indexed differently
   but equally rightly. Each observation's original indices (from M/ISYM)
   and batch fix the relation exactly where merged data cannot: the lattice
   operation under which most observations match is the one used.
3. The observations, matched by original indices and batch (each counted
   from the file's first batch): intensity and sigma ratios, overall and
   against batch -- scaling's differences show there, against the image.
4. The merged data, both merged here the same way (inverse-variance means
   over the Laue group's equivalents; Friedel mates apart with --anomalous):
   per resolution shell, counts, CC, ratio, <I/sigma>, multiplicity.

Only the intensities and the file structure are compared; what the batch
headers hold (orientations, goniometer) is not decoded.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np

_OFF = 512  # index offset for packing: |h|, |k|, |l| below 512
_BATCH_BITS = 20


@dataclass
class Unmerged:
    label: str
    path: str
    spacegroup: str
    cell: tuple
    wavelengths: list
    columns: list  # (label, type)
    nrows: int
    batches: list  # the batch headers' numbers
    d_range: tuple  # (low, high)
    hkl: np.ndarray  # original indices, N x 3 (from M/ISYM)
    batch: np.ndarray  # counted from the file's first batch
    i: np.ndarray
    sigi: np.ndarray
    rotations: np.ndarray  # the space group's rotations, acting h -> h R
    unit_cell: object  # gemmi.UnitCell, for d


def load(path: str, label: str, columns=("I", "SIGI")) -> Unmerged:
    """An unmerged MTZ: its header, and its observations with their original
    indices restored from M/ISYM."""
    import gemmi

    mtz = gemmi.read_mtz_file(path)
    labels = mtz.column_labels()
    for need in ("M/ISYM", "BATCH", *columns):
        if need not in labels:
            raise ValueError(
                f"{path} has no {need} column"
                + (": is it unmerged?" if need in ("M/ISYM", "BATCH") else "")
                + f" (it has {', '.join(labels)})"
            )
    header_columns = [(c.label, c.type) for c in mtz.columns]
    wavelengths = sorted(
        {round(d.wavelength, 6) for d in mtz.datasets if d.wavelength > 0}
    )
    batches = [b.number for b in mtz.batches]
    d_range = (mtz.resolution_low(), mtz.resolution_high())
    if not mtz.switch_to_original_hkl():
        raise ValueError(f"{path}: its indices could not be put back as measured")
    data = np.array(mtz.array, copy=True)
    col = {name: k for k, name in enumerate(labels)}
    hkl = data[:, [col["H"], col["K"], col["L"]]].round().astype(np.int64)
    batch = data[:, col["BATCH"]].round().astype(np.int64)
    batch = batch - (batch.min() if len(batch) else 0)
    rotations = np.array(
        [
            np.array(op.rot, dtype=np.int64) // op.DEN
            for op in mtz.spacegroup.operations().sym_ops
        ]
    )
    c = mtz.cell
    return Unmerged(
        label=label,
        path=path,
        spacegroup=mtz.spacegroup.hm,
        cell=(c.a, c.b, c.c, c.alpha, c.beta, c.gamma),
        wavelengths=wavelengths,
        columns=header_columns,
        nrows=mtz.nreflections,
        batches=batches,
        d_range=d_range,
        hkl=hkl,
        batch=batch,
        i=data[:, col[columns[0]]].astype(np.float64),
        sigi=data[:, col[columns[1]]].astype(np.float64),
        rotations=rotations,
        unit_cell=c,
    )


def _pack(hkl: np.ndarray) -> np.ndarray:
    h = hkl + _OFF
    return (h[..., 0] << 20) | (h[..., 1] << 10) | h[..., 2]


def _unpack(code: np.ndarray) -> np.ndarray:
    return (
        np.stack([(code >> 20) & 1023, (code >> 10) & 1023, code & 1023], axis=-1)
        - _OFF
    )


def canonical(hkl: np.ndarray, rotations: np.ndarray, anomalous: bool) -> np.ndarray:
    """One code per set of equivalents: the largest packed index among each
    reflection's equivalents under the rotations -- and their inverses,
    Friedel's, unless anomalous. Any consistent choice serves to match; it
    need not be CCP4's asymmetric unit."""
    ops = rotations if anomalous else np.concatenate([rotations, -rotations])
    eq = np.einsum("ni,kij->nkj", hkl, ops)
    return _pack(eq).max(axis=1)


def lattice_operations(cell, tolerance: float = 0.01) -> np.ndarray:
    """The lattice's own symmetry, as the integer matrices M (entries -1, 0,
    1, determinant 1) acting h -> h M that keep the reciprocal metric: the
    candidates for how two indexings of one lattice can differ."""
    basis = np.eye(3, dtype=np.int64)
    g = np.empty((3, 3))
    for i in range(3):
        g[i, i] = cell.calculate_1_d2(gemmi_hkl(basis[i]))
    for i, j in itertools.combinations(range(3), 2):
        g[i, j] = g[j, i] = 0.5 * (
            cell.calculate_1_d2(gemmi_hkl(basis[i] + basis[j])) - g[i, i] - g[j, j]
        )
    m = np.array(list(itertools.product((-1, 0, 1), repeat=9)), dtype=np.int64)
    m = m.reshape(-1, 3, 3)
    m = m[np.round(np.linalg.det(m)).astype(np.int64) == 1]
    kept = np.einsum("nij,jk,nlk->nil", m, g, m)
    scale = np.abs(g).max()
    ok = np.all(np.abs(kept - g) <= tolerance * scale, axis=(1, 2))
    return m[ok]


def gemmi_hkl(v) -> list:
    return [int(x) for x in v]


def _observation_codes(o: Unmerged, m: np.ndarray | None = None) -> np.ndarray:
    hkl = o.hkl if m is None else o.hkl @ m
    return (_pack(hkl) << _BATCH_BITS) | (o.batch & ((1 << _BATCH_BITS) - 1))


def _match(codes_a: np.ndarray, codes_b: np.ndarray):
    """Rows of A and of B with the same code, each code once."""
    ua, ia = np.unique(codes_a, return_index=True)
    ub, ib = np.unique(codes_b, return_index=True)
    _, ka, kb = np.intersect1d(ua, ub, assume_unique=True, return_indices=True)
    return ia[ka], ib[kb]


def choose_indexing(a: Unmerged, b: Unmerged, sample: int = 100000, seed: int = 0):
    """The lattice operation M taking B's indices to A's under which most of
    B's observations (a sample of them) match one of A's by original indices
    and batch; with the identity's count beside it."""
    rng = np.random.default_rng(seed)
    pick = (
        rng.choice(len(b.hkl), size=sample, replace=False)
        if len(b.hkl) > sample
        else np.arange(len(b.hkl))
    )
    sub = Unmerged(**{**b.__dict__, "hkl": b.hkl[pick], "batch": b.batch[pick]})
    codes_a = _observation_codes(a)
    best, best_n, identity_n = np.eye(3, dtype=np.int64), -1, 0
    for m in lattice_operations(a.unit_cell):
        n = len(_match(codes_a, _observation_codes(sub, m))[0])
        if np.array_equal(m, np.eye(3, dtype=np.int64)):
            identity_n = n
        if n > best_n:
            best, best_n = m, n
    return best, best_n, identity_n, len(pick)


@dataclass
class Merged:
    code: np.ndarray
    i: np.ndarray
    sigi: np.ndarray
    n: np.ndarray


def merge(o: Unmerged, m: np.ndarray, anomalous: bool) -> Merged:
    """Inverse-variance means over each set of equivalents, observations with
    a positive sigma only; indices taken through M first."""
    ok = o.sigi > 0
    keys = canonical(o.hkl[ok] @ m, o.rotations, anomalous)
    code, inverse = np.unique(keys, return_inverse=True)
    w = 1.0 / o.sigi[ok] ** 2
    sw = np.bincount(inverse, w)
    i = np.bincount(inverse, w * o.i[ok]) / sw
    return Merged(code, i, 1.0 / np.sqrt(sw), np.bincount(inverse))


def _cc(x, y):
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _median(x):
    return float(np.median(x)) if len(x) else float("nan")


def compare(
    a: Unmerged,
    b: Unmerged,
    shells: int = 10,
    anomalous: bool = False,
    csv: str | None = None,
) -> str:
    out = []
    w = max(len(a.label), len(b.label), 8)

    # 1. The files.
    out.append("The files")
    for o in (a, b):
        out.append(f"  {o.label:<{w}} {o.path}")
    rows = [
        ("space group", a.spacegroup, b.spacegroup),
        (
            "cell",
            " ".join(f"{x:.3f}" for x in a.cell[:3])
            + " "
            + " ".join(f"{x:.2f}" for x in a.cell[3:]),
            " ".join(f"{x:.3f}" for x in b.cell[:3])
            + " "
            + " ".join(f"{x:.2f}" for x in b.cell[3:]),
        ),
        (
            "wavelength",
            ", ".join(f"{x:.6f}" for x in a.wavelengths) or "none",
            ", ".join(f"{x:.6f}" for x in b.wavelengths) or "none",
        ),
        ("rows", str(a.nrows), str(b.nrows)),
        (
            "batches",
            f"{len(a.batches)}"
            + (f" ({min(a.batches)} to {max(a.batches)})" if a.batches else ""),
            f"{len(b.batches)}"
            + (f" ({min(b.batches)} to {max(b.batches)})" if b.batches else ""),
        ),
        (
            "resolution",
            f"{a.d_range[0]:.2f} to {a.d_range[1]:.2f} A",
            f"{b.d_range[0]:.2f} to {b.d_range[1]:.2f} A",
        ),
    ]
    out.append(f"  {'':12s} {a.label:>40s} {b.label:>40s}")
    for name, x, y in rows:
        mark = "" if x == y else "  differs"
        out.append(f"  {name:12s} {x:>40s} {y:>40s}{mark}")
    ca, cb = dict(a.columns), dict(b.columns)
    only_a = [f"{k}({ca[k]})" for k in ca if k not in cb]
    only_b = [f"{k}({cb[k]})" for k in cb if k not in ca]
    typed = [f"{k}: {ca[k]} against {cb[k]}" for k in ca if k in cb and ca[k] != cb[k]]
    out.append(f"  columns in both: {sum(1 for k in ca if k in cb)}")
    out.append(f"  only in {a.label}: {', '.join(only_a) or 'none'}")
    out.append(f"  only in {b.label}: {', '.join(only_b) or 'none'}")
    if typed:
        out.append("  types differ: " + "; ".join(typed))

    # 2. The indexing.
    m, n_best, n_identity, n_sample = choose_indexing(a, b)
    identity = np.array_equal(m, np.eye(3, dtype=np.int64))
    out.append("")
    out.append("The indexing")
    ops = len(lattice_operations(a.unit_cell))
    out.append(
        f"  {ops} lattice operations tried; of {n_sample} of {b.label}'s observations, "
        f"{n_best} match {a.label}'s by original indices and batch"
        + (
            " as indexed"
            if identity
            else f" under h' = h M, M = {m.tolist()} ({n_identity} as indexed)"
        )
    )
    if n_best < 0.5 * n_sample:
        out.append(
            "  fewer than half match: the batches may be numbered or ordered "
            "differently, or the data are not the same sweep"
        )

    # 3. The observations.
    ia, ib = _match(_observation_codes(a), _observation_codes(b, m))
    ok = (a.sigi[ia] > 0) & (b.sigi[ib] > 0)
    ia, ib = ia[ok], ib[ok]
    out.append("")
    out.append("The observations, matched by original indices and batch")
    out.append(
        f"  {len(ia)} matched: {len(ia) / max(len(a.i), 1):.1%} of {a.label}'s "
        f"{len(a.i)}, {len(ia) / max(len(b.i), 1):.1%} of {b.label}'s {len(b.i)}"
    )
    ra, rb = a.i[ia], b.i[ib]
    strong = (ra / a.sigi[ia] > 5) & (rb / b.sigi[ib] > 5)
    k = _median(rb[strong] / ra[strong])
    out.append(
        f"  intensity: CC {_cc(ra, rb):.4f}; {b.label}/{a.label}, median over "
        f"{int(strong.sum())} with I/sigma above 5 in both, {k:.4f}"
    )
    out.append(
        f"  sigma: {b.label}/{a.label}, median {_median(b.sigi[ib] / a.sigi[ia]):.4f}"
    )
    if len(ia):
        batch = a.batch[ia]
        edges = np.quantile(batch, np.linspace(0, 1, 11))
        out.append(
            f"  against batch (from each file's first), ten bins: "
            f"{b.label}/{a.label} for the strong, then sigma's"
        )
        out.append(f"  {'batches':>15s} {'n':>8s} {'I ratio':>9s} {'sigma ratio':>12s}")
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = (batch >= lo) & (batch <= hi if hi == edges[-1] else batch < hi)
            s = sel & strong
            out.append(
                f"  {int(lo):>7d}-{int(hi):<7d} {int(sel.sum()):8d} "
                f"{_median(rb[s] / ra[s]):9.4f} "
                f"{_median(b.sigi[ib][sel] / a.sigi[ia][sel]):12.4f}"
            )

    # 4. The merged data.
    ma = merge(a, np.eye(3, dtype=np.int64), anomalous)
    mb = merge(b, m, anomalous)
    _, ka, kb = np.intersect1d(
        ma.code, mb.code, assume_unique=True, return_indices=True
    )
    hkl = _unpack(ma.code[ka])
    d = np.asarray(
        a.unit_cell.calculate_d_array(hkl.astype(np.int32)), dtype=np.float64
    )
    out.append("")
    out.append(
        "The merged data, merged here the same way"
        + (" (Friedel mates apart)" if anomalous else "")
    )
    out.append(
        f"  unique: {len(ma.code)} in {a.label}, {len(mb.code)} in {b.label}, "
        f"{len(ka)} in both"
    )
    order = np.argsort(-d)
    bins = np.array_split(order, max(1, min(shells, len(order))))
    out.append(
        f"  {'d max':>7s} {'d min':>7s} {'n':>7s} {'CC':>7s} {'ratio':>7s} "
        f"{'I/s ' + a.label[:6]:>12s} {'I/s ' + b.label[:6]:>12s} "
        f"{'mult ' + a.label[:5]:>11s} {'mult ' + b.label[:5]:>11s}"
    )

    def row(sel, first, last):
        ia_, ib_ = ka[sel], kb[sel]
        xa, xb = ma.i[ia_], mb.i[ib_]
        s = (xa / ma.sigi[ia_] > 3) & (xb / mb.sigi[ib_] > 3)
        return (
            f"  {first:7.2f} {last:7.2f} {len(sel):7d} {_cc(xa, xb):7.3f} "
            f"{_median(xb[s] / xa[s]):7.3f} {np.mean(xa / ma.sigi[ia_]):12.2f} "
            f"{np.mean(xb / mb.sigi[ib_]):12.2f} {np.mean(ma.n[ia_]):11.2f} "
            f"{np.mean(mb.n[ib_]):11.2f}"
        )

    for sel in bins:
        if len(sel):
            out.append(row(sel, d[sel].max(), d[sel].min()))
    if len(order):
        out.append(row(order, d.max(), d.min()))

    if csv:
        with open(csv, "w") as f:
            f.write(
                f"h,k,l,d,I_{a.label},sigI_{a.label},n_{a.label},I_{b.label},sigI_{b.label},n_{b.label}\n"
            )
            for j in range(len(ka)):
                h_ = hkl[j]
                f.write(
                    f"{h_[0]},{h_[1]},{h_[2]},{d[j]:.4f},{ma.i[ka[j]]:.6g},{ma.sigi[ka[j]]:.6g},"
                    f"{ma.n[ka[j]]},{mb.i[kb[j]]:.6g},{mb.sigi[kb[j]]:.6g},{mb.n[kb[j]]}\n"
                )
        out.append(f"\nWrote {len(ka)} matched unique reflections to {csv}")
    return "\n".join(out)
