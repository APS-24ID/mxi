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
