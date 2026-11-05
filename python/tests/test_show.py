"""mxi show: everything a file holds, or only what is asked for --
--models, --unit-cell (a crystal's cell a line, nothing else), --reflections,
--summary, --step NAME -- in the usual order."""

import os
import subprocess

import pytest

BUILD = os.path.dirname(os.environ.get("MXI_FIND", "")) or os.environ.get(
    "MXI_BUILD", ""
)
MXI = os.path.join(BUILD, "mxi")

pytestmark = pytest.mark.skipif(not os.path.exists(MXI), reason="mxi is not built")


@pytest.fixture
def indexed(tmp_path):
    from test_verify import indexed as make

    make(tmp_path)  # the planted series, a lattice: indexed.expt, indexed.refl
    run = subprocess.run(
        [MXI, "convert", "indexed.expt", "indexed.refl", "-o", "indexed.rflx"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    return tmp_path


def show(where, *words, code=0):
    run = subprocess.run(
        [MXI, "show", "indexed.rflx", *words], capture_output=True, text=True, cwd=where
    )
    assert run.returncode == code, run.stderr
    return run.stdout, run.stderr


def test_the_unit_cells_a_line_each_and_nothing_else(indexed):
    out, _ = show(indexed, "--unit-cell")
    lines = out.splitlines()
    assert lines[0] == "# experiment sweep a b c alpha beta gamma space_group"
    rows = [l.split() for l in lines[1:]]
    assert rows and all(len(r) >= 9 for r in rows), "nine columns, the group last"
    models, _ = show(indexed, "--models")
    for r in rows:
        assert (
            f"experiment {r[0]}: {r[2]} {r[3]} {r[4]} A, {r[5]} {r[6]} {r[7]} deg"
            in models
        )
    assert "indexed.rflx" not in out, "no name: the cells as they are"


def test_each_part_alone_and_together_in_order(indexed):
    everything, _ = show(indexed)
    for header in ("Experiments:", "Reflections:", "Processing:"):
        assert header in everything
    models, _ = show(indexed, "--models")
    assert (
        "Experiments:" in models
        and "Reflections:" not in models
        and "Processing:" not in models
    )
    refl, _ = show(indexed, "--reflections")
    assert "Reflections:" in refl and "Experiments:" not in refl
    summary, _ = show(indexed, "--summary")
    assert "Processing:" in summary and "Reflections:" not in summary
    both, _ = show(indexed, "--summary", "--models")
    assert both.index("Experiments:") < both.index("Processing:"), "the usual order"


def test_one_steps_record(indexed):
    # strong.rflx, as spot finding wrote it, its summary with it.
    def step(*words, code=0):
        run = subprocess.run(
            [MXI, "show", "strong.rflx", *words],
            capture_output=True,
            text=True,
            cwd=indexed,
        )
        assert run.returncode == code, run.stderr
        return run.stdout, run.stderr

    out, _ = step("--step", "find")
    assert "mxi_find" in out and "Experiments:" not in out
    assert out == step("--step", "mxi_find")[0]
    _, err = step("--step", "scale", code=1)
    assert "records no step mxi_scale" in err
    # A file with no summary -- a DIALS pair, converted -- says so.
    out, _ = show(indexed, "--step", "find", code=1)
    assert "no summary in the file" in out
