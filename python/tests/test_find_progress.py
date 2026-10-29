"""mxi find's progress, as integration's: from --progress-from images (100), a
bar on a terminal and a line each tenth of the way otherwise -- here, the
output a pipe; none below the threshold; the spots the same either way."""

import json
import os
import subprocess
import sys

import numpy as np
import pytest

from mxeq import rflx

BUILD = os.path.dirname(os.environ.get("MXI_FIND", "")) or os.environ.get(
    "MXI_BUILD", ""
)
FIND = os.path.join(BUILD, "mxi_find")

pytestmark = pytest.mark.skipif(not os.path.exists(FIND), reason="mxi is not built")


def test_a_long_series_reports_its_progress(tmp_path):
    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(tmp_path), "120"],
        check=True,
        capture_output=True,
    )
    from mxeq.fixtures import nxmx

    nxmx.write_expt(str(tmp_path / "one.expt"), 120)
    doc = json.load(open(tmp_path / "one.expt"))
    doc["experiment"][0]["imageset"] = 0
    doc["imageset"] = [
        {"__id__": "ImageSequence", "template": str(tmp_path / "series.nxs")}
    ]
    json.dump(doc, open(tmp_path / "l.expt", "w"))
    runs = {}
    for name, extra in (("with", []), ("without", ["--progress-from", "1000"])):
        run = subprocess.run(
            [
                FIND,
                str(tmp_path / "l.expt"),
                "-o",
                str(tmp_path / f"{name}.rflx"),
                *extra,
            ],
            capture_output=True,
            text=True,
            cwd=tmp_path,
        )
        assert run.returncode == 0, run.stderr
        runs[name] = run.stdout
    lines = [l for l in runs["with"].splitlines() if "%:" in l]
    assert len(lines) == 9 and "of 120 images" in lines[0]  # 10% .. 90%
    assert not any("%:" in l for l in runs["without"].splitlines())
    a, b = rflx.load_table(str(tmp_path / "with.rflx")), rflx.load_table(
        str(tmp_path / "without.rflx")
    )
    assert a.nrows == b.nrows and np.array_equal(
        np.asarray(a["xyzobs.px.value"]), np.asarray(b["xyzobs.px.value"])
    )
