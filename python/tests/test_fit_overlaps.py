"""mxi integrate fitting overlaps together (the default; --no-fit-overlaps
each alone), the plumbing: reflections whose foregrounds
overlap found, their boxes held, fitted together (fit_together, tested on its
own in tests/test_fit_together.cc), the results put over the fits alone -- on
the CPU and through the device's batches (emulated) alike -- and nothing else
touched. On the planted series, its profile model widened (sigma_b 0.6) so that
its foregrounds overlap: its predictions are not where its spots are, so the
intensities mean nothing here; what is tested is what changes, and what does
not."""

import os
import subprocess

import numpy as np
import pytest

from mxeq import refl, rflx

BUILD = os.path.dirname(os.environ.get("MXI_FIND", "")) or os.environ.get(
    "MXI_BUILD", ""
)
MXI = os.path.join(BUILD, "mxi")

pytestmark = pytest.mark.skipif(not os.path.exists(MXI), reason="mxi is not built")

WIDE = ["--sigma-b", "0.6", "--sigma-m", "0.1", "-j", "2"]


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    from test_verify import indexed

    where = tmp_path_factory.mktemp("overlaps")
    indexed(where)
    return where


def integrate(where, out, *words):
    run = subprocess.run(
        [MXI, "integrate", "indexed.expt", "indexed.refl", "-o", out, *WIDE, *words],
        capture_output=True,
        text=True,
        cwd=where,
    )
    assert run.returncode == 0, run.stderr[-2000:]
    return run.stdout, refl.load(str(where / out))


def outside_identical(a, b, together):
    for k in a.columns:
        x, y = np.asarray(a.columns[k])[~together], np.asarray(b.columns[k])[~together]
        assert np.array_equal(x, y, equal_nan=x.dtype.kind == "f"), k


def test_overlaps_fitted_together_and_nothing_else_touched(planted):
    _, alone = integrate(planted, "alone.rflx", "--no-fit-overlaps")
    out, both = integrate(planted, "together.rflx")
    assert "Overlapping foregrounds:" in out and "Fitted together:" in out
    assert "profile.fitted_together" not in alone.columns, "without it, no new column"
    together = np.asarray(both.columns["profile.fitted_together"]).astype(bool)
    assert together.sum() >= 4, "the widened foregrounds overlap"
    outside_identical(alone, both, together)
    # The members: every one flagged overlapped as before; profile fitting
    # changed, and the summed intensities without their neighbours' counts,
    # each variance at least its variance with no signal.
    flags = np.asarray(both.columns["flags"]).astype(np.int64)
    assert np.all((flags[together] >> 12) & 1)
    for k in ("intensity.prf.value", "intensity.sum.value"):
        assert not np.array_equal(
            np.asarray(alone.columns[k])[together],
            np.asarray(both.columns[k])[together],
        ), k
    v = np.asarray(both.columns["intensity.sum.variance"])[together]
    floor = np.asarray(both.columns["intensity.sum.variance_background"])[together]
    assert np.all(v >= floor) and np.all(v > 0)
    steps = rflx.summary_dict(str(planted / "together.rflx"))["steps"]
    assert steps[-1]["fitted_together"] == int(together.sum())
    assert (
        "fitted_together"
        not in rflx.summary_dict(str(planted / "alone.rflx"))["steps"][-1]
    )


def test_through_the_devices_batches_the_same_fits_together(planted):
    # The device's fits arrive when a batch is collected: the fits together
    # must still be what is written, exactly the CPU's (both fitted on the
    # CPU, in double, from the same boxes).
    _, cpu = integrate(planted, "cpu.rflx")
    _, alone = integrate(
        planted, "emu_alone.rflx", "--gpu-emulate", "--no-fit-overlaps"
    )
    _, emu = integrate(planted, "emu.rflx", "--gpu-emulate")
    together = np.asarray(emu.columns["profile.fitted_together"]).astype(bool)
    assert np.array_equal(
        together, np.asarray(cpu.columns["profile.fitted_together"]).astype(bool)
    )
    outside_identical(alone, emu, together)
    for k in (
        "intensity.prf.value",
        "intensity.prf.variance",
        "intensity.sum.value",
        "intensity.sum.variance",
    ):
        assert np.array_equal(
            np.asarray(cpu.columns[k])[together], np.asarray(emu.columns[k])[together]
        ), k


@pytest.mark.parametrize("words", [["--no-mask-neighbours"], ["--two-pass"]])
def test_where_it_cannot_run_each_alone_said_only_if_asked(planted, words):
    out, t = integrate(planted, "refused.rflx", "--fit-overlaps", *words)
    assert "--fit-overlaps" in out and "profile.fitted_together" not in t.columns
    out, t = integrate(planted, "quiet.rflx", *words)
    assert "--fit-overlaps" not in out and "profile.fitted_together" not in t.columns
