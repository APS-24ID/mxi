"""mxi integrate chooses one pass or two by memory (22 November).

A single pass holds each closed shoebox until the scan block after its own is
complete; on crowded data that is most of its memory -- a 900-image set with c
200 A, 21 GB of footprint in one pass and 7.2 in two, in the same time, too
much for a small laptop. So the single pass's holding is predicted from the
boxes before an image is read, and over the budget -- --max-memory, else half
the machine's memory -- two passes are taken, and said. --one-pass insists.
The two passes so chosen are --two-pass's own, column for column.
"""

import os
import subprocess

import numpy as np
import pytest

from mxeq import refl

BUILD = os.environ.get("MXI_BUILD")
MXI = os.path.join(BUILD, "mxi") if BUILD else None
needs_build = pytest.mark.skipif(not MXI, reason="set MXI_BUILD")
# The planted series' indexed spots have no shoeboxes to estimate a profile
# model from: given, as the other tests on it give it.
MODEL = ["--sigma-b", "0.6", "--sigma-m", "0.1", "-j", "2"]


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    from test_verify import indexed

    where = tmp_path_factory.mktemp("budget")
    indexed(where)
    return where


def integrate(where, out, *words):
    return subprocess.run(
        [MXI, "integrate", "indexed.expt", "indexed.refl", "-o", out]
        + MODEL
        + list(words),
        capture_output=True,
        text=True,
        cwd=where,
    )


def same_tables(a, b):
    for k in b.columns:
        u, v = np.asarray(a[k]), np.asarray(b[k])
        if u.shape != v.shape:
            return False
        if not np.array_equal(u, v, equal_nan=u.dtype.kind == "f"):
            return False
    return a.nrows == b.nrows and set(a.columns) == set(b.columns)


@needs_build
def test_over_the_budget_two_passes_are_taken_and_said(planted):
    run = integrate(planted, "auto.rflx", "--max-memory", "1e-6")
    assert run.returncode == 0, run.stderr[-2000:]
    assert "Two passes over the images" in run.stdout
    assert "--max-memory" in run.stdout and "--one-pass to insist" in run.stdout
    two = integrate(planted, "two.rflx", "--two-pass")
    assert two.returncode == 0, two.stderr[-2000:]
    a = refl.load(str(planted / "auto.rflx"))
    b = refl.load(str(planted / "two.rflx"))
    assert same_tables(a, b), "the two passes chosen are --two-pass's own"


@needs_build
def test_under_the_budget_one_pass_says_nothing(planted):
    run = integrate(planted, "one.rflx", "--max-memory", "1000")
    assert run.returncode == 0, run.stderr[-2000:]
    assert "Two passes" not in run.stdout and "One pass, as asked" not in run.stdout


@needs_build
def test_one_pass_insists_and_says_so(planted):
    run = integrate(planted, "insist.rflx", "--max-memory", "1e-6", "--one-pass")
    assert run.returncode == 0, run.stderr[-2000:]
    assert "One pass, as asked (--one-pass)" in run.stdout
    assert "Two passes" not in run.stdout


@needs_build
def test_one_pass_and_two_pass_together_are_refused(planted):
    run = integrate(planted, "both.rflx", "--one-pass", "--two-pass")
    assert run.returncode == 2
    assert "--one-pass or --two-pass, not both" in run.stderr


@needs_build
def test_timing_says_the_prediction(planted):
    run = integrate(planted, "timed.rflx", "--timing")
    assert run.returncode == 0, run.stderr[-2000:]
    assert "memory, predicted from the boxes: one pass" in run.stdout
