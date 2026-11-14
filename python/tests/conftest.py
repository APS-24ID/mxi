import sys
from pathlib import Path

# fixtures.py sits beside the tests and is imported by name rather than as part
# of the package: it generates test data and has no business shipping.
sys.path.insert(0, str(Path(__file__).parent))

# The programs' locations as absolute paths, before any test changes
# directory: tests run the programs from temporary directories, where a
# relative MXI_FIND -- as the release workflow gave, ../build/mxi_find -- names
# nothing, and 41 of them failed on GitHub (26 October). Read when pytest
# loads this file, before the test modules read the variables themselves.
import os  # noqa: E402

# The programs from the build, where they are not given one by one: MXI_BUILD
# alone had run only the tests that read it, and every test asking MXI_FIND,
# MXI_INTEGRATE and the rest was skipped -- 129 of them on Graeme's Mac (13
# November). Each is the build directory's link to mxi, mxi_<command>; one
# set already stands. MXI_BUILD itself from MXI_FIND, where only that is set.
_PROGRAMS = (
    "import find index refine integrate symmetry scale export max report convert"
).split()
if not os.environ.get("MXI_BUILD") and os.environ.get("MXI_FIND"):
    os.environ["MXI_BUILD"] = os.path.dirname(os.environ["MXI_FIND"])
_build = os.environ.get("MXI_BUILD")
if _build:
    for _command in _PROGRAMS:
        _variable = "MXI_" + _command.upper()
        _program = os.path.join(_build, "mxi_" + _command)
        if not os.environ.get(_variable) and os.path.exists(_program):
            os.environ[_variable] = _program

# Every MXI_ variable naming a file or directory that exists -- the programs,
# the build, the data sets -- not a list of them, which would miss the next.
for _name, _value in list(os.environ.items()):
    if _name.startswith("MXI_") and _value and os.path.exists(_value):
        os.environ[_name] = os.path.abspath(_value)
