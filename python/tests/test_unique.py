"""mxeq unique splits a difference in I/sigma into its parts, and finds the
setting the second data set was indexed in.

Planted: the same observations with every sigma 1.2 times larger, which must
read as -16.7 per cent in 1/sigma and nowhere else; every other observation
dropped, which must read as the root of two in multiplicity; and the second
data set indexed in P3's twin setting, (k, h, -l), which must be found.
"""

import gemmi
import numpy as np

from mxeq import equivalents as eq
from mxeq import unique

HALL = "P 3"


def observations(label="a", seed=1, uniques=600, per=6):
    rng = np.random.default_rng(seed)
    ops = [
        np.array(op.rot, dtype=np.int64) // op.DEN
        for op in gemmi.symops_from_hall(HALL).sym_ops
    ]
    reps = []
    seen = set()
    while len(reps) < uniques:
        h = rng.integers(-15, 16, 3)
        if not h.any():
            continue
        k = int(eq.asu_keys(h[None, :], HALL)[0])
        if k not in seen:
            seen.add(k)
            reps.append(h)
    reps = np.array(reps)
    truth = rng.exponential(2000.0, uniques) + 50.0
    hkl, intensity, variance = [], [], []
    for r, t in zip(reps, truth):
        for _ in range(per):
            h = r @ ops[rng.integers(len(ops))] * rng.choice([-1, 1])
            v = t + 20.0
            hkl.append(h)
            intensity.append(t + np.sqrt(v) * rng.normal())
            variance.append(v)
    hkl = np.array(hkl)
    d = 40.0 / np.linalg.norm(hkl, axis=1)
    return unique.Observations(
        label, HALL, hkl, np.array(intensity), np.array(variance), d
    )


def overall(text):
    lines = text.splitlines()
    start = max(i for i, line in enumerate(lines) if "medians, in per cent" in line)
    row = next(
        line for line in lines[start:] if line.strip().startswith("overall")
    ).split()
    cc, total, mean, mult, inv_sigma, scatter = (float(x) for x in row[1:7])
    return cc, total, mean, mult, inv_sigma, scatter


def relative_scatter(text):
    lines = text.splitlines()
    start = max(i for i, line in enumerate(lines) if "medians, in per cent" in line)
    return float(
        next(
            line for line in lines[start:] if line.strip().startswith("overall")
        ).split()[7]
    )


def test_sigmas_larger_by_a_fifth_read_as_sigma_and_nothing_else():
    a = observations()
    b = unique.Observations("b", HALL, a.hkl, a.intensity, a.variance * 1.44, a.d)
    cc, total, mean, mult, inv_sigma, scatter = overall(unique.compare(a, b, shells=3))
    assert cc > 0.9999 and abs(mean) < 0.01 and abs(mult) < 0.01 and abs(scatter) < 0.01
    assert abs(inv_sigma - 20.0) < 0.05  # a's 1/sigma is 1.2 times b's
    assert abs(total - 20.0) < 0.05


def test_half_the_observations_read_as_multiplicity():
    a = observations()
    keep = np.arange(len(a.intensity)) % 2 == 0
    b = unique.Observations(
        "b", HALL, a.hkl[keep], a.intensity[keep], a.variance[keep], a.d[keep]
    )
    _, _, _, mult, inv_sigma, _ = overall(unique.compare(a, b, shells=3))
    assert abs(mult - 100 * (np.sqrt(2) - 1)) < 1.0, mult
    assert abs(inv_sigma) < 0.5


def test_the_twin_setting_is_found():
    a = observations()
    twin = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
    b = unique.Observations("b", HALL, a.hkl @ twin, a.intensity, a.variance, a.d)
    text = unique.compare(a, b, shells=3)
    assert "reindexed by" in text, text.splitlines()[:3]
    cc = overall(text)[0]
    assert cc > 0.9999


def test_overloaded_observations_are_counted_where_a_table_flags_them():
    # DIALS flags a reflection with a pixel over the trusted range; mxi does
    # not. Flagged observations scaling kept are counted, per shell and overall.
    a = observations()
    over = np.zeros(len(a.intensity), bool)
    over[:30] = True
    b = unique.Observations("b", HALL, a.hkl, a.intensity, a.variance, a.d, over, 45)
    text = unique.compare(a, b, shells=3)
    assert "b: 45 rows flagged overloaded, 30 of them used in scaling" in text
    assert "a: 0 rows flagged overloaded" in text and "never sets the flag" in text
    lines = text.splitlines()
    start = lines.index("  b, over the reflections both have")
    row = next(
        line for line in lines[start:] if line.strip().startswith("overall")
    ).split()
    assert int(row[-1]) == 30


def test_a_difference_of_overall_scale_cancels_in_the_relative_scatter():
    # One data set the other times 0.8, sigmas and all: the mean intensity and
    # the scatter differ by a quarter, the relative scatter, I/sigma and
    # sigma's part not at all.
    a = observations()
    b = unique.Observations("b", HALL, a.hkl, 0.8 * a.intensity, 0.64 * a.variance, a.d)
    text = unique.compare(a, b, shells=3)
    cc, total, mean, mult, inv_sigma, scatter = overall(text)
    assert abs(mean - 25.0) < 0.05 and abs(scatter - 25.0) < 0.05
    assert abs(relative_scatter(text)) < 0.05
    assert abs(total) < 0.05 and abs(inv_sigma + 20.0) < 0.05


def test_two_absurd_observations_do_not_choose_a_reindexing():
    # Graeme's mxi table held two observations of one reflection at -25 sigma,
    # from a defective pixel in their backgrounds. Merged, that reflection was
    # so far out that the Pearson correlation of the two data sets fell to
    # 0.44, and a wrong operator -- not even a symmetry of the lattice -- was
    # taken for the indexing; every statistic after it was nonsense. The
    # operator is chosen by a rank correlation, which a few such values
    # cannot move.
    a = observations()
    intensity = a.intensity.copy()
    target = np.flatnonzero(
        np.all(a.hkl == a.hkl[0], axis=1) | np.all(a.hkl == -a.hkl[0], axis=1)
    )
    rows = target[:2] if len(target) >= 2 else np.array([0, 1])
    intensity[rows] = -150.0 * np.max(a.intensity)
    b = unique.Observations("b", HALL, a.hkl, intensity, a.variance, a.d)
    text = unique.compare(a, b, shells=3)
    assert "b's indices as they are" in text, text.splitlines()[:3]
