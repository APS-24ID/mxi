"""The processing's summary a .rflx carries (docs/output.md): each program's
step, after the steps its input carried. mxi_find twice -- on a synthetic
series' experiment list, then on its own output -- must leave two steps, the
first carried. Needs MXI_FIND."""

import os
import subprocess
import sys

import pytest

from mxeq import rflx

FIND = os.environ.get("MXI_FIND")
needs_find = pytest.mark.skipif(not FIND, reason="set MXI_FIND")


@needs_find
def test_a_step_is_added_after_those_the_input_carried(tmp_path):
    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(tmp_path / "s"), "12"],
        check=True,
        capture_output=True,
    )
    expt = str(tmp_path / "s" / "series.expt")
    master = str(tmp_path / "s" / "series.nxs")
    first = subprocess.run(
        [FIND, "-e", expt, master, "-o", "first.rflx"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert first.returncode == 0, first.stderr + first.stdout
    one = rflx.summary_dict(str(tmp_path / "first.rflx"))
    assert [s["program"] for s in one["steps"]] == ["mxi_find"]
    # Its figures are this run's: the spots written.
    table = rflx.load_table(str(tmp_path / "first.rflx"))
    assert one["steps"][0]["spots"] == len(table.columns["flags"])
    # Spots per image a list a sweep: Graeme's four sweeps had kept only one.
    per_image = one["steps"][0]["spots_per_image"]
    assert isinstance(per_image, list) and isinstance(per_image[0], list)
    # By resolution only where the geometry reads in full, as here it does not;
    # where it does, fifty bins equal in 1/d^2.
    by_d = one["steps"][0].get("spots_by_resolution")
    assert by_d is None or len(by_d["spots"]) == 50

    second = subprocess.run(
        [FIND, "-e", "first.rflx", master, "-o", "second.rflx"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert second.returncode == 0, second.stderr + second.stdout
    two = rflx.summary_dict(str(tmp_path / "second.rflx"))
    assert [s["program"] for s in two["steps"]] == ["mxi_find", "mxi_find"]
    assert two["steps"][0]["time"] == one["steps"][0]["time"]


REPORT = os.environ.get("MXI_REPORT")


@pytest.mark.skipif(not (FIND and REPORT), reason="set MXI_FIND and MXI_REPORT")
def test_the_report_draws_each_step_it_finds(tmp_path):
    # mxi_report on the spot finder's file: its step's section, its plots,
    # each plot's data and layout JSON a browser's plotly can take.
    import json
    import re

    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(tmp_path / "s"), "12"],
        check=True,
        capture_output=True,
    )
    s = tmp_path / "s"
    found = subprocess.run(
        [
            FIND,
            "-e",
            str(s / "series.expt"),
            str(s / "series.nxs"),
            "-o",
            "strong.rflx",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert found.returncode == 0, found.stderr + found.stdout
    r = subprocess.run(
        [REPORT, "strong.rflx", "-o", "report.html"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr + r.stdout
    page = (tmp_path / "report.html").read_text(encoding="utf-8")
    assert "<h2>Spot finding" in page
    calls = re.findall(
        r"Plotly\.newPlot\('plot\d+', (\[.*?\]), (\{.*?\}), \{responsive: true\}\)",
        page,
    )
    assert calls, "no plots"
    titles = [json.loads(layout)["title"] for _, layout in calls]
    assert "Spots per image" in titles
    for traces, _ in calls:
        assert all(len(t["x"]) == len(t["y"]) for t in json.loads(traces))
