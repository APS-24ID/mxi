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

# Every MXI_ variable naming a file or directory that exists -- the programs,
# the build, the data sets -- not a list of them, which would miss the next.
for _name, _value in list(os.environ.items()):
    if _name.startswith("MXI_") and _value and os.path.exists(_value):
        os.environ[_name] = os.path.abspath(_value)
