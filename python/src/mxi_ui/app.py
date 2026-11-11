"""mxi-ui's entry: mxi-ui [DIRECTORY], or mxi-ui from the wheel."""

import argparse
import sys

from mxi_ui import deps
from mxi_ui.version import VERSION


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mxi-ui",
        description="An interface to mxi: the workflow, its images and its results.",
    )
    parser.add_argument(
        "directory", nargs="?", help="the working directory (the current one)"
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="say whether it can run here, and what is missing",
    )
    parser.add_argument("--version", action="version", version=f"mxi-ui {VERSION}")
    args = parser.parse_args(argv)

    found = deps.problems(need_display=not args.check)
    if args.check:
        if found:
            for problem in found:
                print(f"mxi-ui: {problem}")
            return 1
        from mxi_ui.locate import find_mxi

        print(
            f"mxi-ui {VERSION}: ready (Python {sys.version.split()[0]}; mxi: {find_mxi()})"
        )
        return 0
    if found:
        for problem in found:
            print(f"mxi-ui: {problem}", file=sys.stderr)
        return 1

    import tkinter as tk

    from mxi_ui.window import MainWindow

    root = tk.Tk()
    MainWindow(root, args.directory)
    root.mainloop()
    return 0
