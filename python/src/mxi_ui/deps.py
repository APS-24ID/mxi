"""Whether mxi-ui can run here, and what to install where it cannot. Checked
before anything imports tkinter, so that a missing piece is a sentence naming
what to install, not a traceback. tkinter comes with Python, but some systems
package it apart: Debian and Ubuntu's python3-tk, Homebrew's python-tk."""

import importlib.util
import os
import platform
import sys

PACKAGES = ("numpy", "h5py", "hdf5plugin", "matplotlib")


def _linux_family(os_release="/etc/os-release"):
    try:
        text = open(os_release).read()
    except OSError:
        return ""
    fields = dict(
        line.split("=", 1)
        for line in text.splitlines()
        if "=" in line and not line.startswith("#")
    )
    ids = (
        " ".join(fields.get(k, "").strip('"') for k in ("ID", "ID_LIKE"))
        .lower()
        .split()
    )
    for family in (
        "debian",
        "ubuntu",
        "fedora",
        "rhel",
        "centos",
        "arch",
        "suse",
        "opensuse",
    ):
        if family in ids:
            return family
    return ""


def tkinter_hint(system=None, family=None, version=None):
    """What to install for tkinter, for this system."""
    system = system or platform.system()
    major, minor = (version or sys.version_info)[:2]
    if system == "Darwin":
        return (
            f"With Homebrew's Python: brew install python-tk@{major}.{minor}. "
            "The python.org installers include it; with conda: conda install tk."
        )
    if system == "Linux":
        family = _linux_family() if family is None else family
        if family in ("debian", "ubuntu"):
            return "sudo apt install python3-tk"
        if family in ("fedora", "rhel", "centos"):
            return "sudo dnf install python3-tkinter"
        if family == "arch":
            return "sudo pacman -S tk"
        if family in ("suse", "opensuse"):
            return "sudo zypper install python3-tk"
        return "install your distribution's Tk package for Python (often python3-tk)"
    if system == "Windows":
        return "run the python.org installer again, with 'tcl/tk and IDLE' ticked"
    return "install Tk for this Python"


def problems(need_display=True):
    """What stops mxi-ui running here, each with what to do: [] if nothing."""
    found = []
    if importlib.util.find_spec("tkinter") is None:
        found.append(
            f"tkinter is missing from this Python ({sys.executable}): {tkinter_hint()}"
        )
    else:
        try:
            import tkinter

            if tkinter.TkVersion < 8.6:
                found.append(
                    f"Tk {tkinter.TkVersion} is too old: 8.6 or later is needed"
                )
        except ImportError as e:
            found.append(f"tkinter cannot be loaded ({e}): {tkinter_hint()}")
    missing = [p for p in PACKAGES if importlib.util.find_spec(p) is None]
    if missing:
        found.append(
            f"{', '.join(missing)} missing: pip install mxi-python (the wheel from the release) "
            f"or pip install {' '.join(missing)}"
        )
    if (
        need_display
        and platform.system() == "Linux"
        and not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))
    ):
        found.append("no display: run it on a desktop, or over ssh with -X")
    return found
