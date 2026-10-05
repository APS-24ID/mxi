"""mxeq background-model fits R(s) . G(phi) to the backgrounds, flags those
that do not fit, and filters the table by them.

Planted: backgrounds made from a known radial curve with a water-like bump, a
rotation wave of 15 per cent over a full turn, the polarisation, solid angle and
efficiency the experiment's geometry gives, counting noise and a 3 per cent
intrinsic spread; then a cluster in a shadow, a tenth, and a ring in a flare,
three times. Every planted shadow and flare must be flagged and few of the rest;
G's shape recovered; the filtered table's flags right.
"""

import json

import numpy as np

from mxeq import background_model as bm
from mxeq import refl

PX, DIST, CENTRE = 0.075, 200.0, (1000.0, 1000.0)


def write(tmp_path, n=40000, seed=8):
    rng = np.random.default_rng(seed)
    e = {
        "__id__": "ExperimentList",
        "experiment": [{"__id__": "Experiment", "beam": 0, "detector": 0, "scan": 0}],
        "beam": [
            {
                "direction": [0.0, 0.0, 1.0],
                "wavelength": 1.0,
                "polarization_normal": [0.0, 1.0, 0.0],
                "polarization_fraction": 0.99,
            }
        ],
        "detector": [
            {
                "panels": [
                    {
                        "origin": [-CENTRE[0] * PX, CENTRE[1] * PX, -DIST],
                        "fast_axis": [1.0, 0, 0],
                        "slow_axis": [0, -1.0, 0],
                        "pixel_size": [PX, PX],
                        "image_size": [2000, 2000],
                        "mu": 3.0,
                        "thickness": 0.45,
                    }
                ]
            }
        ],
        "scan": [{"image_range": [1, 3600], "properties": {"oscillation": [0.0, 0.1]}}],
    }
    (tmp_path / "a.expt").write_text(json.dumps(e))
    g = bm.geometry(str(tmp_path / "a.expt"))
    radius = np.sqrt(rng.uniform(20.0**2, 950.0**2, n))
    az = rng.uniform(-np.pi, np.pi, n)
    x, y = CENTRE[0] + radius * np.cos(az), CENTRE[1] + radius * np.sin(az)
    two_theta = np.arctan(radius * PX / DIST)
    d = 1.0 / (2.0 * np.sin(two_theta / 2.0))
    s = 1.0 / d
    r_true = 1.5 * np.exp(-s / 0.4) + 0.6 * np.exp(-(((s - 0.3) / 0.05) ** 2)) + 0.05
    phi = rng.uniform(0.0, 360.0, n)
    g_true = 1.0 + 0.15 * np.cos(np.radians(2 * phi))
    p, omega, q = bm.known_factors(g, x, y)
    n_pix = rng.integers(500, 3000, n).astype(float)
    mean = r_true * g_true * p * omega * q * rng.lognormal(0.0, 0.03, n)
    bg = rng.poisson(mean * n_pix) / n_pix
    shadow = (radius < 70.0) & (np.abs(az - 2.5) < 0.5)
    flare = (radius > 90.0) & (radius < 110.0) & (np.abs(az) < 1.0)
    bg[shadow] *= 0.1
    bg[flare] *= 3.0
    t = refl.ReflectionTable(nrows=n)
    flags = np.full(n, (1 << 8) | (1 << 9) | 1, dtype=np.int64)
    for name, values, kind in (
        ("background.mean", bg, "double"),
        ("num_pixels.background", n_pix.astype(np.int32), "int"),
        ("d", d, "double"),
        ("xyzcal.px", np.column_stack([x, y, phi / 0.1]), "vec3<double>"),
        (
            "xyzcal.mm",
            np.column_stack([x * PX, y * PX, np.radians(phi)]),
            "vec3<double>",
        ),
        ("flags", flags, "std::size_t"),
    ):
        t.columns[name] = values
        t.types[name] = kind
    refl.write(str(tmp_path / "a.refl"), t)
    return shadow, flare, phi, g_true


def test_planted_shadows_and_flares_are_flagged_and_few_else(tmp_path):
    shadow, flare, _, _ = write(tmp_path)
    r = bm.run(str(tmp_path / "a.expt"), str(tmp_path / "a.refl"), z_max=5.0)
    assert shadow.sum() > 20 and flare.sum() > 20
    assert r.fit.flagged_low[shadow].mean() > 0.95
    assert r.fit.flagged_high[flare].mean() > 0.95
    clean = ~(shadow | flare)
    assert (r.fit.flagged_low | r.fit.flagged_high)[clean].mean() < 0.005
    assert 0.02 < r.fit.tau < 0.045  # the planted 3 per cent


def test_the_rotation_is_recovered(tmp_path):
    write(tmp_path)
    r = bm.run(str(tmp_path / "a.expt"), str(tmp_path / "a.refl"))
    grid, fitted = r.fit.phi_curve
    truth = 1.0 + 0.15 * np.cos(np.radians(2 * grid))
    assert np.corrcoef(np.log(fitted), np.log(truth))[0, 1] > 0.98
    amplitude = (fitted.max() - fitted.min()) / (truth.max() - truth.min())
    assert 0.85 < amplitude < 1.15


def test_the_filtered_table_has_its_flags_right(tmp_path):
    # By default only the too-low go; with reject_high, the too-high too.
    write(tmp_path)
    base = (1 << 8) | (1 << 9) | 1
    for reject_high in (False, True):
        r = bm.run(
            str(tmp_path / "a.expt"), str(tmp_path / "a.refl"), reject_high=reject_high
        )
        flags = np.asarray(r.table.columns["flags"]).astype(np.int64)
        gone = r.fit.flagged_low | (r.fit.flagged_high if reject_high else False)
        assert r.fit.flagged_high.any()
        assert np.all((flags[gone] & ((1 << 8) | (1 << 9))) == 0)
        assert np.all(flags[gone] & (1 << 24))
        assert np.all(flags[~gone] == base)
    assert (
        "background.expected" in r.table.columns and "background.z" in r.table.columns
    )


def test_the_polarisation_is_one_along_the_beam_and_cos_squared_in_the_plane():
    g = bm.Geometry(
        np.array([-10.0, 10.0, -100.0]),
        np.array([1.0, 0, 0]),
        np.array([0, -1.0, 0]),
        (0.1, 0.1),
        np.array([0, 0, 1.0]),
        np.array([0, 1.0, 0]),
        1.0,
        0.0,
        0.0,
        0.0,
        0.1,
    )
    # Pixel (100, 100) is on the beam. Pixel (1100, 100) is 100 mm along +x at
    # 100 mm: 2theta 45 degrees, in the plane the beam is polarised in, where a
    # fully polarised beam scatters cos^2(2theta) = 0.5 of what it does forward.
    p, _, _ = bm.known_factors(g, np.array([100.0, 1100.0]), np.array([100.0, 100.0]))
    assert abs(p[0] - 1.0) < 1e-9
    assert abs(p[1] - np.cos(np.radians(45.0)) ** 2) < 1e-9


def test_the_figure_maps_the_removed_and_the_kept(tmp_path):
    pytest = __import__("pytest")
    pytest.importorskip("matplotlib")
    from mxeq.plots import background as plot

    shadow, flare, _, _ = write(tmp_path)
    r = bm.run(str(tmp_path / "a.expt"), str(tmp_path / "a.refl"))
    px = np.asarray(r.table.columns["xyzcal.px"]).reshape(-1, 3)
    bg = np.asarray(r.table.columns["background.mean"], float).ravel()
    plot.draw_model(r.fit, bg, px[:, 0], px[:, 1], CENTRE, str(tmp_path / "m.png"))
    assert (tmp_path / "m.png").stat().st_size > 20000
    assert len(r.fit.flagged_shell) == len(bg) and not r.fit.flagged_shell.any()


def test_annotate_only_adds_the_columns_and_changes_no_flag(tmp_path):
    write(tmp_path)
    r = bm.run(str(tmp_path / "a.expt"), str(tmp_path / "a.refl"), annotate_only=True)
    flags = np.asarray(r.table.columns["flags"]).astype(np.int64)
    assert np.all(flags == ((1 << 8) | (1 << 9) | 1))
    assert (r.fit.flagged_low | r.fit.flagged_high).any()  # found, but left in
    assert (
        "background.z" in r.table.columns and "background.expected" in r.table.columns
    )
    assert "annotate-only" in r.report


def test_reject_inner_leaves_out_the_innermost_whole_and_stops_at_its_edge(tmp_path):
    # Beside a backstop, half the backgrounds a tenth, half three times: no
    # background there is normal, z cannot tell good from bad, and one
    # attenuated to the model's level would pass. --reject-inner takes the
    # region whole, outward until a shell's z spread is normal.
    write(tmp_path)
    t = refl.load(str(tmp_path / "a.refl"))
    px = np.asarray(t.columns["xyzcal.px"]).reshape(-1, 3)
    radius = np.hypot(px[:, 0] - CENTRE[0], px[:, 1] - CENTRE[1])
    rng = np.random.default_rng(3)
    inside = radius < 60.0
    bg = np.asarray(t.columns["background.mean"], float).copy()
    bg[inside] *= np.where(rng.uniform(size=inside.sum()) < 0.5, 0.1, 3.0)
    t.columns["background.mean"] = bg
    refl.write(str(tmp_path / "a.refl"), t)
    r = bm.run(str(tmp_path / "a.expt"), str(tmp_path / "a.refl"), reject_inner=3.0)
    gone = r.fit.flagged_shell | r.fit.flagged_low | r.fit.flagged_high
    assert gone[inside].all()
    beyond = radius > 90.0
    assert r.fit.flagged_shell[beyond].mean() < 0.01
    assert "the innermost left out to d" in r.report
