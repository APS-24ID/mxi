#!/usr/bin/env python3
"""mxi-view: the images of an experiment list, with the spots and
predictions mxi found on them, and the detector's untrusted regions drawn and
saved -- mxi-ui's viewer in a window of its own (mxi_ui/viewer.py,
docs/ui.md, docs/mask.md). Save writes masked.rflx in the directory it is run
in (-o to name another file), through mxi untrusted, and mxi_mask.json.

It had been napari's, until mxi-ui's viewer drew the regions and summed the
images too: Tk now, with numpy, h5py and hdf5plugin, as mxi-ui is; no napari.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

VERSION = "@PROJECT_VERSION@"
# From the wheel (mxi-python), no build to stamp it in: the wheel's own.
if VERSION.startswith("@"):
    try:
        from importlib.metadata import version as _version

        VERSION = _version("mxi-python")
    except Exception:  # noqa: BLE001 -- run from a checkout, not installed
        VERSION = "unknown"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="mxi-view",
        description="The images of an experiment list, with the spots and "
        "predictions mxi found on them; and the detector's untrusted regions "
        "drawn and saved to masked.rflx through mxi untrusted.",
    )
    parser.add_argument(
        "experiments", help="a .rflx or .expt: imported, strong, integrated ..."
    )
    parser.add_argument(
        "-o", default="masked.rflx", help="where Save writes (masked.rflx)"
    )
    parser.add_argument(
        "--sum", type=int, default=1, help="images summed for each shown (1)"
    )
    parser.add_argument(
        "--sweep", type=int, default=0, help="the sweep shown first (0)"
    )
    parser.add_argument("--version", action="version", version=f"mxi-view {VERSION}")
    args = parser.parse_args(argv)
    if not os.path.exists(args.experiments):
        print(f"mxi-view: no such file: {args.experiments}", file=sys.stderr)
        return 2
    from mxi_ui import deps

    problems = [p for p in deps.problems() if "matplotlib" not in p]
    if problems:
        print("mxi-view cannot run here:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        return 2
    import tkinter as tk

    from mxi_ui.viewer import Viewer

    root = tk.Tk()
    root.withdraw()
    try:
        viewer = Viewer(
            root,
            args.experiments,
            f"mxi-view {os.path.basename(args.experiments)}",
            output=args.o,
            where=os.getcwd(),
            summing=args.sum,
            sweep=args.sweep,
        )
    except Exception as e:  # noqa: BLE001 -- said, whatever it is
        print(f"mxi-view: {e}", file=sys.stderr)
        root.destroy()
        return 1
    # The viewer's window is the program's: closed, it ends.
    viewer.top.bind(
        "<Destroy>", lambda e: root.quit() if e.widget is viewer.top else None
    )
    root.mainloop()
    try:
        root.destroy()
    except tk.TclError:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
