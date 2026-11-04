"""mxi process: the eight steps in turn, each option to the steps it
belongs to; a failing step stops it, named, with its log and its command."""

import os
import subprocess
import sys

import pytest

from test_ui import takes

BUILD = os.path.dirname(os.environ.get("MXI_FIND", "")) or os.environ.get(
    "MXI_BUILD", ""
)
MXI = os.path.join(BUILD, "mxi")

pytestmark = pytest.mark.skipif(not os.path.exists(MXI), reason="mxi is not built")

ALL = [
    "--gpu",
    "-j",
    "4",
    "--cell",
    "78,78,78,90,90,90",
    "--space-group",
    "I213",
    "--d-min-auto",
    "--no-joint",
]


def dry(*words):
    run = subprocess.run(
        [MXI, "process", *words, "--dry-run"], capture_output=True, text=True
    )
    assert run.returncode == 0, run.stderr
    return run.stdout.splitlines()


def test_each_option_goes_to_its_steps():
    assert dry("a.nxs", "b.nxs", *ALL) == [
        "mxi import a.nxs b.nxs -j 4",
        "mxi find imported.rflx --gpu -j 4",
        "mxi index strong.rflx --known-cell 78,78,78,90,90,90 --known-space-group I213 --no-joint -j 4",
        "mxi refine indexed.rflx -j 4",
        "mxi integrate refined.rflx --gpu -j 4",
        "mxi symmetry integrated.rflx -j 4",
        "mxi scale symmetrized.rflx --space-group I213 --d-min-auto -j 4",
        "mxi export scaled.rflx -j 4",
    ]
    assert dry("a.nxs") == [
        "mxi import a.nxs",
        "mxi find imported.rflx",
        "mxi index strong.rflx",
        "mxi refine indexed.rflx",
        "mxi integrate refined.rflx",
        "mxi symmetry integrated.rflx",
        "mxi scale symmetrized.rflx",
        "mxi export scaled.rflx",
    ]


def test_a_space_group_without_a_cell_goes_to_scaling_alone():
    lines = dry("a.nxs", "--space-group", "I213")
    assert lines[2] == "mxi index strong.rflx"
    assert lines[6] == "mxi scale symmetrized.rflx --space-group I213"


def test_every_option_it_passes_is_one_its_step_takes():
    # Cannot drift from the steps: each flag on each command, in its --help.
    for line in dry("a.nxs", *ALL):
        words = line.split()
        step = words[1]
        helped = subprocess.run([MXI, step, "--help"], capture_output=True, text=True)
        for flag in [w for w in words[2:] if w.startswith("-")]:
            assert takes(helped.stdout + helped.stderr, flag), f"mxi {step} {flag}"


def test_a_failing_step_stops_it_named_with_its_log_and_command(tmp_path):
    # The planted series is no NXmx mxi import reads: step 1 fails.
    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(tmp_path), "6"],
        check=True,
        capture_output=True,
    )
    run = subprocess.run(
        [MXI, "process", "series.nxs", "-j", "2"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 1
    assert "stopped at step 1 of 8, mxi import, which failed (exit 1)" in run.stderr
    assert (
        "mxi_import.log" in run.stderr
        and "  mxi import series.nxs -j 2\n" in run.stderr
    )
    assert "step 1 of 8: mxi import series.nxs -j 2" in run.stdout
    assert "step 2 of 8" not in run.stdout and not (tmp_path / "strong.rflx").exists()
    assert (tmp_path / "mxi_process.log").exists()


def test_it_needs_the_images_and_refuses_what_it_does_not_take():
    none = subprocess.run([MXI, "process"], capture_output=True, text=True)
    assert none.returncode == 2 and "master files" in none.stderr
    odd = subprocess.run(
        [MXI, "process", "a.nxs", "--max-lattices", "2"], capture_output=True, text=True
    )
    assert odd.returncode == 2 and "--max-lattices" in odd.stderr


@pytest.fixture
def strong(tmp_path):
    from test_verify import indexed

    indexed(tmp_path)  # the planted series, its spots found: strong.rflx
    return tmp_path


def dry_in(where, *words):
    run = subprocess.run(
        [MXI, "process", *words, "--dry-run"], capture_output=True, text=True, cwd=where
    )
    return run.returncode, run.stdout.splitlines(), run.stderr


def test_a_file_picks_up_after_the_step_that_wrote_it(strong):
    code, lines, notes = dry_in(strong, "strong.rflx", "-j", "2")
    assert code == 0 and lines == [
        "mxi index strong.rflx -j 2",
        "mxi refine indexed.rflx -j 2",
        "mxi integrate refined.rflx -j 2",
        "mxi symmetry integrated.rflx -j 2",
        "mxi scale symmetrized.rflx -j 2",
        "mxi export scaled.rflx -j 2",
    ]
    assert "strong.rflx, written by mxi_find: from mxi index on" in notes


def test_a_cluster_is_scaled_afresh_and_exported(strong):
    # A cluster carries the summary it was cut from, scaling's last; its
    # writer, mxi correlation, says it is to be scaled again.
    import shutil

    import h5py

    shutil.copy(strong / "strong.rflx", strong / "cluster_0.rflx")
    with h5py.File(strong / "cluster_0.rflx", "a") as f:
        f.attrs["creator"] = "mxi_correlation"
    code, lines, notes = dry_in(
        strong,
        "cluster_0.rflx",
        "--d-min-auto",
        "--gpu",
        "--cell",
        "78,78,78,90,90,90",
        "--space-group",
        "I213",
    )
    assert code == 0 and lines == [
        "mxi scale cluster_0.rflx --space-group I213 --d-min-auto",
        "mxi export scaled.rflx",
    ]
    assert "--cell not used" in notes and "--gpu not used" in notes
    assert "--space-group not used" not in notes, "scaling takes it"


def test_without_its_writer_the_summarys_last_step_and_with_neither_refused(strong):
    import h5py

    with h5py.File(strong / "strong.rflx", "a") as f:
        del f.attrs["creator"]
    code, lines, notes = dry_in(strong, "strong.rflx")
    assert (
        code == 0 and lines[0] == "mxi index strong.rflx"
    ), "the summary: mxi_find last"
    with h5py.File(strong / "strong.rflx", "a") as f:
        del f["summary"]
    code, lines, notes = dry_in(strong, "strong.rflx")
    assert code == 1 and "says neither the step that wrote it" in notes


def test_a_failure_picking_up_says_how_to_go_on(strong):
    run = subprocess.run(
        [MXI, "process", "strong.rflx", "-j", "2"],
        capture_output=True,
        text=True,
        cwd=strong,
    )
    assert run.returncode != 0
    failed = run.stderr.split("stopped at step ")[1].split(",")[0]  # "2 of 6"
    assert failed.endswith("of 6")
    assert "or all that remains:\n  mxi process " in run.stderr and "-j 2" in run.stderr
