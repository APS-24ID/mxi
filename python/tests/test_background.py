"""mxeq background: each reflection's background against the median at its
resolution, by shell, and the figure.

Planted: a background falling smoothly with resolution, a polarisation wave of
ten per cent with the azimuth, and a cluster beside the beam in a shadow, a
tenth of what it should be. The table must show the cluster in the innermost
shell and nowhere else, and the ratio's median must be 1.
"""

import json

import numpy as np
import pytest

from mxeq import background, refl


def planted(tmp_path, n=6000):
    rng = np.random.default_rng(4)
    px, distance, centre = 0.075, 200.0, (1000.0, 1000.0)
    e = {
        "__id__": "ExperimentList",
        "experiment": [{"__id__": "Experiment", "beam": 0, "detector": 0}],
        "beam": [{"direction": [0.0, 0.0, 1.0], "wavelength": 1.0}],
        "detector": [
            {
                "panels": [
                    {
                        "origin": [-centre[0] * px, centre[1] * px, -distance],
                        "fast_axis": [1.0, 0.0, 0.0],
                        "slow_axis": [0.0, -1.0, 0.0],
                        "pixel_size": [px, px],
                        "image_size": [2000, 2000],
                    }
                ]
            }
        ],
    }
    (tmp_path / "a.expt").write_text(json.dumps(e))
    radius = np.sqrt(rng.uniform(30.0**2, 900.0**2, n))
    azimuth = rng.uniform(-np.pi, np.pi, n)
    x, y = centre[0] + radius * np.cos(azimuth), centre[1] + radius * np.sin(azimuth)
    two_theta = np.arctan(radius * px / distance)
    d = 1.0 / (2.0 * np.sin(two_theta / 2.0))
    bg = (2.0 * np.exp(-((1.0 / d) ** 2) / 0.2) + 0.1) * (
        1.0 + 0.1 * np.cos(2 * azimuth)
    )
    bg *= rng.normal(1.0, 0.03, n)
    shadow = (radius < 60.0) & (np.abs(azimuth - 0.5) < 0.6)
    bg[shadow] *= 0.1
    t = refl.ReflectionTable(nrows=n)
    for name, values, kind in (
        ("background.mean", bg, "double"),
        ("d", d, "double"),
        ("xyzcal.px", np.column_stack([x, y, rng.uniform(0, 100, n)]), "vec3<double>"),
        ("num_pixels.background", np.full(n, 200, dtype=np.int32), "int"),
    ):
        t.columns[name] = values
        t.types[name] = kind
    refl.write(str(tmp_path / "a.refl"), t)
    return int(shadow.sum())


def test_a_shadow_shows_in_the_innermost_shell_and_nowhere_else(tmp_path):
    n_shadow = planted(tmp_path)
    assert n_shadow > 5
    b = background.load(str(tmp_path / "a.expt"), str(tmp_path / "a.refl"))
    assert abs(b.centre[0] - 1000.0) < 1e-6 and abs(b.centre[1] - 1000.0) < 1e-6
    r = background.ratios(b, background.curve(b))
    assert abs(np.median(r) - 1.0) < 0.02
    text = background.table(b, shells=10)
    block = text.split("by resolution")[1].split("by azimuth")[0]
    rows = [line.split() for line in block.splitlines() if "%" in line]
    below = [float(row[-2].rstrip("%")) for row in rows]
    assert below[0] > 0.0 and all(v == 0.0 for v in below[1:]), below


def test_the_figure_is_drawn(tmp_path):
    pytest.importorskip("matplotlib")
    from mxeq.plots import background as plot

    planted(tmp_path)
    b = background.load(str(tmp_path / "a.expt"), str(tmp_path / "a.refl"))
    plot.draw(b, str(tmp_path / "b.png"), (50.0, 3.0), 200.0)
    assert (tmp_path / "b.png").stat().st_size > 10000


def test_a_range_limits_the_table_and_gives_azimuth_and_image(tmp_path):
    planted(tmp_path)
    b = background.load(str(tmp_path / "a.expt"), str(tmp_path / "a.refl"))
    text = background.table(b, shells=4, d_range=(100.0, 3.0))
    assert "in 100 to 3 A" in text and "by azimuth around the beam, 100 to 3 A" in text
    assert "by image, 100 to 3 A" in text
    block = text.split("by resolution")[1].split("by azimuth")[0]
    ds = [float(line.split()[0]) for line in block.splitlines() if "%" in line]
    assert max(ds) <= 100.0
