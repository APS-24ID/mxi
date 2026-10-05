"""mxeq attenuation: each observation good or attenuated, from its equivalents
and its background together, attenuation only lowering.

Planted: a reflection like ferritin's (3,1,1) -- 13 unharmed near 950 on normal
backgrounds and 14 attenuated, signal and background dimmed together by T from
0.01 to 0.6 -- whose 14 must go and 13 stay, Ibar near 950; thousands of
ordinary reflections, under half a per cent flagged; the flare, backgrounds three
times with intensities right, none flagged; and a weak reflection left alone.
"""

import numpy as np

from mxeq import attenuation as at


def counts_variance(intensity, background):
    return np.maximum(intensity, 0.0) + 100.0 * background + 10.0


def group(rng, key, truth, good, attenuated, bg=25.0, flare=0):
    i, v, b, k = [], [], [], []
    for _ in range(good):
        bb = bg * rng.uniform(0.85, 1.15)
        var = counts_variance(truth, bb)
        i.append(truth + np.sqrt(var) * rng.normal())
        v.append(var)
        b.append(bb)
    for _ in range(attenuated):
        t = rng.uniform(0.01, 0.6)
        bb = t * bg * rng.uniform(0.85, 1.15)
        var = counts_variance(t * truth, bb)
        i.append(t * truth + np.sqrt(var) * rng.normal())
        v.append(var)
        b.append(bb)
    for _ in range(flare):
        bb = 3.0 * bg
        var = counts_variance(truth, bb)
        i.append(truth + np.sqrt(var) * rng.normal())
        v.append(var)
        b.append(bb)
    n = good + attenuated + flare
    return [key] * n, i, v, b


def assemble(groups):
    keys, i, v, b = (
        np.concatenate([np.asarray(g[x], float) for g in groups]) for x in range(4)
    )
    return keys.astype(np.int64), i, v, b


def test_the_upper_group_is_the_truth_and_the_attenuated_go():
    rng = np.random.default_rng(1)
    keys, i, v, b = assemble([group(rng, 0, 950.0, 13, 14)])
    r = at.judge(keys, i, v, b, np.full(len(i), np.nan), np.full(len(i), 50.0))
    assert not r["flagged"][:13].any()
    assert r["flagged"][13:].all()
    assert abs(r["ibar"][0] - 950.0) < 30.0


def test_ordinary_reflections_and_the_flare_are_left_alone():
    rng = np.random.default_rng(2)
    gs = [
        group(rng, k, rng.exponential(2000.0) + 200.0, 12, 0, bg=rng.uniform(0.5, 5.0))
        for k in range(3000)
    ]
    gs += [group(rng, 5000 + k, 950.0, 10, 0, flare=4) for k in range(50)]
    keys, i, v, b = assemble(gs)
    r = at.judge(keys, i, v, b, np.full(len(i), np.nan), np.full(len(i), 3.0))
    ordinary = keys < 5000
    assert r["flagged"][ordinary].mean() < 0.005
    assert not r["flagged"][~ordinary].any()


def test_a_weak_reflection_is_judged_by_its_background_alone():
    # Too weak for its equivalents to judge -- perhaps weak because attenuated
    # -- so the background decides: the dimmed go, the normal stay; and a weak
    # reflection on normal backgrounds keeps every observation.
    rng = np.random.default_rng(3)
    keys, i, v, b = assemble([group(rng, 0, 5.0, 6, 3, bg=2.0)])
    r = at.judge(keys, i, v, b, np.full(len(i), np.nan), np.full(len(i), 3.0))
    assert not r["judged"].any()
    assert not r["flagged"][:6].any()
    assert r["flagged"][6:][b[6:] < 0.5 * np.percentile(b, 75)].all()
    keys, i, v, b = assemble([group(rng, 0, 5.0, 9, 0, bg=2.0)])
    r = at.judge(keys, i, v, b, np.full(len(i), np.nan), np.full(len(i), 3.0))
    assert not r["flagged"].any()


def test_a_reflection_observed_once_is_judged_against_the_background_model():
    # No group to compare with: the background model's expectation is the
    # reference where the table has it; without it, nothing can be said.
    keys = np.array([0, 1, 2], np.int64)
    i = np.array([0.5, 900.0, 0.5])
    v = np.array([10.0, 1000.0, 10.0])
    b = np.array([0.4, 1.0, 0.4])
    expected = np.array([1.0, 1.0, np.nan])
    r = at.judge(keys, i, v, b, expected, np.full(3, 50.0))
    assert r["flagged"][0] and not r["flagged"][1] and not r["flagged"][2]


def test_attenuation_only_lowers_so_the_upper_group_wins_when_backgrounds_say_nothing():
    # 6 at 950 and 12 at 475, every background normal: a model letting
    # attenuation go either way takes the larger group; attenuation only
    # lowers, so the upper is the truth.
    rng = np.random.default_rng(5)
    i, v, b = [], [], []
    for truth, n in ((950.0, 6), (475.0, 12)):
        for _ in range(n):
            var = counts_variance(truth, 25.0)
            i.append(truth + np.sqrt(var) * rng.normal())
            v.append(var)
            b.append(25.0 * rng.uniform(0.9, 1.1))
    i, v, b = np.array(i), np.array(v), np.array(b)
    r = at.judge(
        np.zeros(len(i), np.int64),
        i,
        v,
        b,
        np.full(len(i), np.nan),
        np.full(len(i), 50.0),
    )
    assert abs(r["ibar"][0] - 950.0) < 40.0
    assert r["flagged"][6:].all() and not r["flagged"][:6].any()
