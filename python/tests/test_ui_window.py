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
    assert tabs == ["Setup", "Log", "Summary", "Images", "Plots"]
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


@pytest.fixture
def stand_in(tmp_path, monkeypatch):
    """An mxi that prints what it was given, so that Run's command is seen."""
    fake = tmp_path / "fake-mxi"
    fake.write_text(
        f"#!{sys.executable}\nimport sys\nprint('given:', ' '.join(sys.argv[1:]))\n"
    )
    fake.chmod(0o755)
    monkeypatch.setenv("MXI", str(fake))
    return fake


def choose(w, key):
    w.tree.selection_set(key)
    w.root.update()


def test_a_form_builds_the_command_as_it_is_filled(window):
    choose(window, "scale")
    form = window.forms["scale"]
    form["options"]["--anomalous"].set(True)
    form["options"]["--d-min"].set("1.8")
    form["threads"].set("4")
    window.root.update()
    assert (
        window.command.get()
        == "mxi scale symmetrized.rflx --anomalous --d-min 1.8 -j 4"
    )
    assert str(window.run_button["state"]) == "normal" and window.trouble.get() == ""
    form["options"]["--d-min"].set("two")
    window.root.update()
    assert "Resolution limit" in window.trouble.get()
    assert str(window.run_button["state"]) == "disabled"
    window.copy_command()
    assert window.root.clipboard_get() == window.command.get()


def test_a_forms_entries_are_kept_from_step_to_step(window):
    choose(window, "index")
    window.forms["index"]["options"]["--max-lattices"].set("2")
    choose(window, "refine")
    choose(window, "index")
    assert window.command.get() == "mxi index strong.rflx --max-lattices 2"


def test_run_runs_the_forms_command(stand_in, tmp_path):
    from mxi_ui.window import MainWindow

    (tmp_path / "imported.rflx").write_text("")
    root = tkinter.Tk()
    root.withdraw()
    w = MainWindow(root, str(tmp_path))
    try:
        choose(w, "find")
        w.forms["find"]["options"]["--min-spot-size"].set("5")
        w.root.update()
        w.run_selected()
        drive(w, lambda: "[exit 0]" in log_text(w))
        assert "$ mxi find imported.rflx --min-spot-size 5" in log_text(w)
        assert "given: find imported.rflx --min-spot-size 5" in log_text(w)
    finally:
        w.runner.cancel()
        root.destroy()


def test_the_inputs_browse_button_is_shown(window):
    # Packed after its entry, it had been dropped where the form's column was
    # narrower than the row asked -- seen in a screenshot, not by a test.
    window.root.deiconify()
    window.root.geometry("1100x720")
    choose(window, "scale")
    for _ in range(10):
        window.root.update()
    buttons = [
        b
        for frame in window.fields.winfo_children()
        for b in frame.winfo_children()
        if str(b.winfo_class()) == "TButton" and b.cget("text") == "Browse..."
    ]
    assert len(buttons) == 1 and buttons[0].winfo_ismapped()
    assert buttons[0].winfo_width() > 20


def test_a_fields_explanation_is_a_tooltip_not_a_line_beside_it(window):
    from mxi_ui.tooltip import DELAY_MS

    choose(window, "scale")
    beside = [
        w for w in window.fields.grid_slaves() if int(w.grid_info()["column"]) == 2
    ]
    assert not beside, "no explanations beside the fields"
    tooltip = next(t for t in window.tooltips if t.text.startswith("--d-min:"))
    assert tooltip.text == "--d-min: nothing beyond this d"
    tooltip.widget.event_generate("<Enter>")
    end = time.time() + (DELAY_MS + 2000) / 1000
    while time.time() < end and not tooltip.showing():
        window.root.update()
        time.sleep(0.02)
    assert tooltip.showing(), "shown once the pointer rests"
    tooltip.widget.event_generate("<Leave>")
    window.root.update()
    assert not tooltip.showing(), "gone when it leaves"
    texts = [t.text for t in window.tooltips]
    assert any(t.startswith("-j:") for t in texts) and any(
        t.startswith("-o:") for t in texts
    )
    assert any(t.startswith("--anomalous:") for t in texts), "a tick box's too"


def test_run_stays_shown_and_named_however_narrow_the_window(window):
    window.root.deiconify()
    window.root.minsize(1, 1)
    choose(window, "scale")
    for width in (1100, 700, 560):
        window.root.geometry(f"{width}x600")
        for _ in range(10):
            window.root.update()
        assert window.run_button.winfo_ismapped(), width
        assert window.run_button.winfo_width() > 30, width
    assert window.run_button.cget("text") == "Run scale"
    choose(window, "index")
    assert window.run_button.cget("text") == "Run index"


def test_the_shortcut_runs_the_step(stand_in, tmp_path):
    from mxi_ui.window import RUN_KEY, MainWindow

    (tmp_path / "imported.rflx").write_text("")
    root = tkinter.Tk()
    w = MainWindow(root, str(tmp_path))
    try:
        choose(w, "find")
        root.update()
        root.focus_force()
        root.event_generate(f"<{RUN_KEY}-Return>", when="tail")
        drive(w, lambda: "[exit 0]" in log_text(w))
        assert "given: find imported.rflx" in log_text(w)
    finally:
        w.runner.cancel()
        root.destroy()


def test_the_sidebars_rows_fit_the_font(window):
    from tkinter import font as tkfont

    linespace = tkfont.nametofont("TkDefaultFont").metrics("linespace")
    height = int(tkinter.ttk.Style(window.root).lookup("Treeview", "rowheight"))
    assert height >= linespace


def test_run_is_in_view_on_every_tab_and_a_step_keeps_the_tab(window):
    window.root.deiconify()
    window.root.geometry("900x480")
    choose(window, "index")
    for tab in window.tabs.tabs():
        window.tabs.select(tab)
        for _ in range(5):
            window.root.update()
        bottom = window.run_button.winfo_rooty() + window.run_button.winfo_height()
        assert (
            window.run_button.winfo_ismapped() and window.run_button.winfo_height() > 15
        )
        assert bottom <= window.root.winfo_rooty() + window.root.winfo_height()
    window.tabs.select(window.summary.master)
    choose(window, "refine")
    assert window.tabs.select() == str(window.summary.master), "the tab kept"
    assert "mxi refine indexed.rflx" in window.summary.get("1.0", "end")


def test_a_steps_summary_and_plots_come_from_its_file(tmp_path, monkeypatch):
    from ui_records import install

    from mxi_ui.window import MainWindow

    install(monkeypatch, tmp_path)
    root = tkinter.Tk()
    root.withdraw()
    w = MainWindow(root, str(tmp_path))
    try:
        choose(w, "scale")
        text = w.summary.get("1.0", "end")
        assert "Merging statistics" in text and "High resolution limit" in text
        assert "Resolution at CC½ 0.3 (Å)" in text
        assert w.canvas is not None and len(w.figure.axes) >= 4
        titles = [a.get_title() for a in w.figure.axes]
        assert "R factors by resolution" in titles
        (tmp_path / "scaled.rflx").unlink()
        choose(w, "refine")
        choose(w, "scale")
        labels = [
            c.cget("text")
            for c in w.plots.winfo_children()
            if c.winfo_class() == "TLabel"
        ]
        assert labels == ["Plots appear once scale has run."] and w.canvas is None
    finally:
        root.destroy()


def test_the_plots_are_drawn_to_the_tabs_size(tmp_path, monkeypatch):
    # Sized from its own wish, the figure had been drawn wider than the tab,
    # its right-hand plots cut off; it is the canvas's size now, which fills
    # the tab, when the tab is shown and when the window changes size.
    from ui_records import install

    from mxi_ui.window import MainWindow

    install(monkeypatch, tmp_path)
    root = tkinter.Tk()
    w = MainWindow(root, str(tmp_path))
    try:
        choose(w, "scale")
        w.tabs.select(w.plots)
        for geometry in ("1100x720", "900x600"):
            root.geometry(geometry)
            for _ in range(15):
                root.update()
            canvas = w.canvas.get_tk_widget()
            drawn = w.figure.get_size_inches() * w.figure.dpi
            assert abs(drawn[0] - canvas.winfo_width()) <= 2, geometry
            assert canvas.winfo_width() <= w.plots.winfo_width(), geometry
    finally:
        root.destroy()


def test_a_file_without_the_steps_record_says_which_it_has(tmp_path, monkeypatch):
    from ui_records import chain, install

    from mxi_ui.window import MainWindow

    records = chain()
    del records["index"]["program"]  # as if written some other way
    records["index"]["program"] = "mxi_convert"
    install(monkeypatch, tmp_path, records)
    root = tkinter.Tk()
    root.withdraw()
    w = MainWindow(root, str(tmp_path))
    try:
        choose(w, "index")
        text = w.summary.get("1.0", "end")
        assert "indexed.rflx records no step of mxi index: it records" in text
        assert "mxi_find" in text and "mxi_convert" in text
    finally:
        root.destroy()
