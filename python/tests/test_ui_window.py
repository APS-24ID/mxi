"""mxi ui's window, where there is a display and tkinter: it opens on the
first step not done, runs a command into its log, cancels one, and marks a
step done when its output appears."""

import os
import subprocess
import sys
import time

import pytest

tkinter = pytest.importorskip("tkinter")


def _can_open_a_window():
    probe = "import tkinter; tkinter.Tk().destroy()"
    try:
        return (
            subprocess.run(
                [sys.executable, "-c", probe], capture_output=True, timeout=20
            ).returncode
            == 0
        )
    except subprocess.TimeoutExpired:
        return False


pytestmark = pytest.mark.skipif(
    not _can_open_a_window(), reason="no display for a window"
)


@pytest.fixture
def window(tmp_path):
    from mxi_ui.window import MainWindow

    for name in ("imported.rflx", "strong.rflx"):
        (tmp_path / name).write_text("")
    root = tkinter.Tk()
    root.withdraw()
    w = MainWindow(root, str(tmp_path))
    root.update()
    yield w
    w.runner.cancel()
    root.destroy()


def drive(w, until, seconds=20):
    end = time.time() + seconds
    while time.time() < end:
        w.root.update()
        if until():
            return
        time.sleep(0.02)
    raise AssertionError("timed out")


def log_text(w):
    return w.log.get("1.0", "end")


def test_it_opens_on_the_first_step_not_done(window):
    assert window.selected.key == "index"
    assert (
        window.tree.set("find", "state") == "done"
        and window.tree.set("index", "state") == ""
    )
    tabs = [window.tabs.tab(t, "text") for t in window.tabs.tabs()]
    assert tabs == ["Images", "Log", "Summary", "Plots"]
    assert "mxi index strong.rflx" in window.summary.get("1.0", "end")


def test_a_command_runs_into_the_log_and_marks_its_step_done(window, tmp_path):
    script = "open('indexed.rflx', 'w').close(); print('indexed', flush=True)"
    assert window.run([sys.executable, "-c", script], title="a stand-in for index")
    drive(window, lambda: "[exit 0]" in log_text(window))
    assert "indexed" in log_text(window) and window.status.get() == "Done"
    assert window.tree.set("index", "state") == "done"
    assert str(window.run_button["state"]) == "normal"


def test_a_running_command_can_be_cancelled(window):
    window.run(
        [
            sys.executable,
            "-c",
            "import time; print('waiting', flush=True); time.sleep(60)",
        ]
    )
    drive(window, lambda: "waiting" in log_text(window))
    window.cancel()
    drive(window, lambda: "[exit" in log_text(window))
    assert window.status.get().startswith("Failed")
