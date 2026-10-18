"""Every program answers --help: its usage on standard output, exit 0.

Seven programs asked the shared parser has("--help"), which it never set --
it reports a help request as args.help -- so their --help was an error, exit
2 and nothing on standard output (October 2026). Needs MXI_FIND, for the
build directory every program is in.
"""

import os
import pathlib
import subprocess

import pytest

FIND = os.environ.get("MXI_FIND")
BUILD = pathlib.Path(FIND).parent if FIND else None
PROGRAMS = (
    sorted(
        p.name
        for p in BUILD.glob("mxi_*")
        if p.is_file()
        and os.access(p, os.X_OK)
        and p.suffix == ""
        and p.name != "mxi_tests"  # the test runner, not a program
    )
    if BUILD
    else []
)


@pytest.mark.skipif(not FIND, reason="set MXI_FIND")
@pytest.mark.parametrize("program", PROGRAMS)
def test_help_is_an_answer_not_an_error(program, tmp_path):
    out = subprocess.run(
        [str(BUILD / program), "--help"], capture_output=True, text=True, cwd=tmp_path
    )
    assert out.returncode == 0, f"{program} --help: exit {out.returncode}"
    assert "usage" in out.stdout.lower(), f"{program} --help printed no usage"
