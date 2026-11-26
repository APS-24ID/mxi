"""mxeq lattice-map: offsets of indexed spots from their lattice points, in
the crystal's axes and in each reflection's frame (26 November). Planted
offsets, each found where it belongs."""

import numpy as np

from mxeq import expt, latticemap, refl


def rotation(axis, angle):
    axis = np.asarray(axis, dtype=float) / np.linalg.norm(axis)
    k = np.array(
        [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]]
    )
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


def setup(n=4000, seed=1, scan_varying=False):
    rng = np.random.default_rng(seed)
    real = rotation([0.3, -0.5, 0.8], 0.7) @ np.diag([50.0, 56.0, 61.0])
    real = real.T  # a, b, c as rows
    crystal = expt.Crystal(a=real)
    A = crystal.setting
    if scan_varying:
        # The crystal turning 0.5 degrees across the scan, about y.
        points = np.array(
            [rotation([0, 1, 0], np.radians(0.5) * i / 100) @ A for i in range(101)]
        )
        crystal = expt.Crystal(a=real, a_at_scan_points=points)
    hkl = rng.integers(-20, 21, size=(n, 3))
    hkl = hkl[np.any(hkl != 0, axis=1)]
    z = rng.uniform(0, 100, len(hkl))
    e = expt.Experiment(
        crystal=crystal,
        scan=expt.Scan(image_range=(1, 100), oscillation=(0.0, 0.1)),
        goniometer=expt.Goniometer(rotation_axis=np.array([1.0, 0.0, 0.0])),
    )
    settings = latticemap._settings(crystal, e.scan, z)
    r = np.einsum("nij,nj->ni", settings, hkl.astype(float))
    return e, hkl, z, r


def table(hkl, z, r):
    t = refl.ReflectionTable(nrows=len(hkl))
    t.columns["miller_index"] = hkl
    t.columns["rlp"] = r
    t.columns["id"] = np.zeros(len(hkl), dtype=int)
    t.columns["xyzobs.px.value"] = np.stack([np.zeros_like(z), np.zeros_like(z), z], 1)
    return t


def gather(e, hkl, z, r):
    return latticemap.gather(table(hkl, z, r), expt.ExperimentList(experiments=[e]))


def test_a_perfect_lattice_has_no_offsets():
    o = gather(*setup())
    assert np.abs(o.dhkl).max() < 1e-9
    assert np.nanmax(np.abs(np.concatenate([o.radial, o.along, o.across]))) < 1e-9


def test_a_cell_too_large_shows_radially_and_in_proportion_to_the_index():
    e, hkl, z, r = setup()
    o = gather(e, hkl, z, r * 1.001)
    assert np.allclose(o.radial, 1e-3, atol=1e-9)
    assert np.nanmax(np.abs(o.along)) < 1e-6 and np.nanmax(np.abs(o.across)) < 1e-6
    assert np.allclose(o.dhkl, 1e-3 * hkl, atol=1e-9)


def test_a_turn_about_the_rotation_axis_shows_along_the_rotation():
    e, hkl, z, r = setup()
    turn = np.radians(0.05)
    o = gather(e, hkl, z, r @ rotation([1, 0, 0], turn).T)
    pred = r / np.linalg.norm(r, axis=1)[:, None]
    reach = np.linalg.norm(np.cross([1.0, 0.0, 0.0], pred), axis=1)
    far = reach > 0.05
    assert np.allclose(o.along[far] / reach[far], np.degrees(turn), rtol=1e-3)
    assert np.nanmax(np.abs(o.radial)) < 1e-6
    assert np.nanmax(np.abs(o.across)) < 1e-4 * np.degrees(turn) * 100


def test_a_scan_varying_crystal_is_taken_at_each_spots_image():
    o = gather(*setup(scan_varying=True))
    assert np.abs(o.dhkl).max() < 1e-9


def test_random_spots_show_as_a_tail():
    e, hkl, z, r = setup(n=8000)
    rng = np.random.default_rng(2)
    r = r + rng.normal(0, 2e-5, r.shape)  # a core
    n = len(r) // 10
    r[:n] = r[:n] + rng.uniform(-3e-3, 3e-3, (n, 3))  # a tenth, anywhere near
    o = gather(e, hkl, z, r)
    tail = latticemap.summary(o.dhkl[:, 0])["beyond_3"]
    assert 0.07 < tail < 0.12, tail
    text = latticemap.report(o, shells=4)
    assert "beyond 3 widths" in text and "histogram" in text
