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
    # One thread and four: progress is counted as each chunk of frames
    # finishes, and a chunk is 4 x threads + 4 frames, so the tenths reached
    # depend on the threads -- every tenth on one, fewer on four. Expecting all
    # nine, this passed on one core and failed on GitHub's four (1.1.3-rc1).
    for name, extra in (
        ("one", ["-j", "1"]),
        ("four", ["-j", "4"]),
        ("without", ["--progress-from", "1000"]),
    ):
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
    for name in ("one", "four"):
        lines = [l for l in runs[name].splitlines() if "%:" in l]
        percent = [int(l.split("%")[0]) for l in lines]
        assert lines and all("of 120 images" in l for l in lines), name
        assert percent == sorted(set(percent)) and all(p < 100 for p in percent), name
    one = [l for l in runs["one"].splitlines() if "%:" in l]
    assert len(one) == 9, "on one thread, chunks of 8 frames: every tenth"
    assert not any("%:" in l for l in runs["without"].splitlines())
    a = rflx.load_table(str(tmp_path / "one.rflx"))
    for other in ("four", "without"):
        b = rflx.load_table(str(tmp_path / f"{other}.rflx"))
        assert a.nrows == b.nrows and np.array_equal(
            np.asarray(a["xyzobs.px.value"]), np.asarray(b["xyzobs.px.value"])
        ), other
