"""Every command takes -j N and --threads N alike, and says so in its help;
a count that is not a number is refused (APS-24ID/mxi#1)."""

import os
import re
import subprocess

import pytest

BUILD = os.path.dirname(os.environ.get("MXI_FIND", "")) or os.environ.get(
    "MXI_BUILD", ""
)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
COMMANDS = [
    c
    for c in re.findall(
        r'^\s*\{"([a-z]+)", mxi_', open(os.path.join(ROOT, "apps/mxi.cc")).read(), re.M
    )
]

pytestmark = pytest.mark.skipif(
    not os.path.exists(os.path.join(BUILD, "mxi")), reason="mxi is not built"
)


def test_every_command_is_found():
    assert len(COMMANDS) == 23 and "find" in COMMANDS and "refine" in COMMANDS


@pytest.mark.parametrize("command", COMMANDS)
def test_every_command_takes_j_and_threads(command):
    mxi = os.path.join(BUILD, "mxi")
    for flag in ("-j", "--threads"):
        run = subprocess.run(
            [mxi, command, flag, "2", "--help"], capture_output=True, text=True
        )
        assert run.returncode == 0, (flag, run.stderr)
        assert "-j, --threads N" in run.stdout + run.stderr, "said in its help"
    bad = subprocess.run([mxi, command, "-j", "many"], capture_output=True, text=True)
    assert bad.returncode != 0 and "needs a number of threads" in bad.stderr


@pytest.mark.parametrize(
    "old", ["--fft-threads", "--jacobian-threads", "--normal-threads"]
)
def test_the_one_part_counts_are_gone(old):
    for command in ("index", "refine"):
        run = subprocess.run(
            [os.path.join(BUILD, "mxi"), command, old, "2"],
            capture_output=True,
            text=True,
        )
        assert run.returncode != 0 and "unknown option" in run.stderr
