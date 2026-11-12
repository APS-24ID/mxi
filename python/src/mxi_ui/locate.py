"""Where mxi is, to run it: the MXI environment variable; this Python's own
programs' directory, where the wheel installs mxi; beside this package, as a
build lays it out (build/mxi beside build/mxi_ui); else the PATH. No tkinter
here."""

import os
import shutil
import sys
import sysconfig


def find_mxi():
    given = os.environ.get("MXI")
    if given:
        return given
    name = "mxi.exe" if os.name == "nt" else "mxi"  # docs/windows.md
    # Beside this package, as a build lays it out (build/mxi, build/mxi_ui):
    # first, so that a build's mxi-ui runs that build's mxi.
    here = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    for where in (
        here,
        # Installed from the wheel: mxi among this Python's own programs --
        # the environment's bin, Scripts on Windows -- whether or not the
        # environment is the PATH's (Graeme, 8 November: everything in the
        # wheel).
        sysconfig.get_path("scripts"),
        os.path.dirname(sys.executable),
    ):
        candidate = os.path.join(where or "", name)
        if where and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return shutil.which("mxi")
