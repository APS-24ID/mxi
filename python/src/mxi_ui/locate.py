"""Where mxi is, to run it: the MXI environment variable; beside this
package, as in a release archive (bin/mxi beside bin/mxi_ui); else the
PATH. No tkinter here."""

import os
import shutil


def find_mxi():
    given = os.environ.get("MXI")
    if given:
        return given
    here = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    # mxi.exe on Windows (docs/windows.md); shutil.which tries it itself.
    beside = os.path.join(here, "mxi.exe" if os.name == "nt" else "mxi")
    if os.path.isfile(beside) and os.access(beside, os.X_OK):
        return beside
    return shutil.which("mxi")
