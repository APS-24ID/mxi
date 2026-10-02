"""mxeq failures shows where profile fitting fails.

A table with failures planted past image 50 only, of one code: the counts by
code must be right, and the table by position in the scan must read none failing
before that image and all after -- what a cause looks like, against a variable
that has nothing to do with it, flat.
"""

import numpy as np

from mxeq import failures, refl


def test_a_planted_cause_shows_against_the_scan_and_not_elsewhere(tmp_path):
    n = 4000
    rng = np.random.default_rng(2)
    z = rng.uniform(0, 100, n)
    t = refl.ReflectionTable(nrows=n)
    columns = {
        "xyzcal.px": (
            np.column_stack([rng.uniform(0, 2000, n), rng.uniform(0, 2000, n), z]),
            "vec3<double>",
        ),
        "zeta": (rng.uniform(0.1, 1.0, n), "double"),
        "d": (rng.uniform(1.2, 10.0, n), "double"),
        "intensity.sum.value": (rng.uniform(0, 1000, n), "double"),
        "intensity.sum.variance": (np.full(n, 100.0), "double"),
        "profile.failure": (np.where(z >= 50, 4, 0).astype(np.int32), "int"),
    }
    bbox = np.zeros((n, 6), dtype=np.int32)
    bbox[:, 1] = 9
    bbox[:, 3] = 9
    bbox[:, 4] = np.floor(z).astype(np.int32) - 3
    bbox[:, 5] = bbox[:, 4] + 7
    columns["bbox"] = (bbox, "int6")
    for name, (values, kind) in columns.items():
        t.columns[name] = values
        t.types[name] = kind
    path = tmp_path / "planted.refl"
    refl.write(str(path), t)

    out = failures.analyse(str(path))
    lines = out.splitlines()
    failing = int(np.sum(z >= 50))
    assert any(line.split()[:2] == [str(failing), "4"] for line in lines), out
    assert any(line.split()[:2] == [str(n - failing), "0"] for line in lines), out
    scan = lines[lines.index("  position in the scan (image)") + 2 :][:12]
    percents = [float(line.split()[-1]) for line in scan]
    starts = [float(line.split()[0].strip("[,")) for line in scan]
    assert all(p == 0.0 for p, s in zip(percents, starts) if s < 49), scan
    assert all(p == 100.0 for p, s in zip(percents, starts) if s >= 51), scan
    zeta = lines[lines.index("  |zeta|") + 2 :][:10]
    assert all(
        30.0 < float(line.split()[-1]) < 70.0 for line in zeta
    ), zeta  # flat: no cause
