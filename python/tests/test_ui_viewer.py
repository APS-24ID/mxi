"""mxi ui's image viewer on a planted series, where there is a display: the
markers drawn are those mxi_tools.Markers puts on each image, each on its
spot at any zoom and pan; the images stepped through; the readout; an indexed
file's indexed spots apart; the main window's View images."""

import os
import subprocess
import sys
import time

import pytest

tkinter = pytest.importorskip("tkinter")

from test_ui_window import _can_open_a_window  # noqa: E402

BUILD = os.path.dirname(os.environ.get("MXI_FIND", "")) or os.environ.get(
    "MXI_BUILD", ""
)
MXI = os.path.join(BUILD, "mxi")

pytestmark = [
    pytest.mark.skipif(not _can_open_a_window(), reason="no display for a window"),
    pytest.mark.skipif(not os.path.exists(MXI), reason="mxi is not built"),
]


def settle(root, n=10):
    for _ in range(n):
        root.update()
        time.sleep(0.01)


@pytest.fixture
def strong(tmp_path):
    from test_verify import indexed

    indexed(tmp_path)  # the planted series, a whole geometry, its spots found
    return tmp_path


def open_viewer(path):
    from mxi_ui.viewer import Viewer

    root = tkinter.Tk()
    root.withdraw()
    v = Viewer(root, str(path))
    v.top.geometry("800x700")
    settle(root)
    v.fit()
    settle(root)
    return root, v


def centres(v, tag):
    out = []
    for item in v.canvas.find_withtag(tag):
        x0, y0, x1, y1 = v.canvas.coords(item)
        out.append(((x0 + x1) / 2, (y0 + y1) / 2))
    return sorted(out)


def expected(v, key):
    yx, _ = v.layers[0][key].on(v._image() - 1)
    w, h = v._canvas_size()
    pts = [v.view.to_canvas(x + 0.5, y + 0.5) for y, x in yx]
    return sorted(p for p in pts if -10 <= p[0] <= w + 10 and -10 <= p[1] <= h + 10)


def test_each_marker_is_on_its_spot_at_any_zoom(strong):
    root, v = open_viewer(strong / "strong.rflx")
    try:
        v.image.set(3)
        v.show()
        settle(root)
        for action in (None, "in", "in", "pan", "out"):
            if action == "in":
                v.view.zoom(1.25, 400, 300)
            elif action == "out":
                v.view.zoom(0.8, 400, 300)
            elif action == "pan":
                v.view.pan(37, -21)
            v.draw()
            settle(root, 3)
            drawn, wanted = centres(v, "unindexed"), expected(v, "unindexed")
            assert len(drawn) == len(wanted) > 0, action
            assert all(
                abs(a - c) < 1e-6 and abs(b - d) < 1e-6
                for (a, b), (c, d) in zip(drawn, wanted)
            ), action
            # And on the spot in the picture: the darkest pixel near a marker.
            x, y = v.view.to_detector(*drawn[len(drawn) // 2])
            patch = v.frame[
                max(0, int(y) - 3) : int(y) + 4, max(0, int(x) - 3) : int(x) + 4
            ]
            assert (
                patch.max()
                == v.frame[
                    max(0, int(y) - 6) : int(y) + 7, max(0, int(x) - 6) : int(x) + 7
                ].max()
                > 0
            ), action
    finally:
        root.destroy()


def test_the_images_step_within_the_sweep(strong):
    root, v = open_viewer(strong / "strong.rflx")
    try:
        assert (v.first, v.last) == (1, 6) and v._image() == 1
        v.step(-1)
        assert v._image() == 1
        for _ in range(10):
            v.step(1)
        assert v._image() == 6
    finally:
        root.destroy()


def test_the_readout_under_the_cursor(strong):
    root, v = open_viewer(strong / "strong.rflx")
    try:
        cx, cy = v.view.to_canvas(100.5, 80.5)

        class Event:
            x, y = cx, cy

        v._hover(Event)
        text = v.readout.get()
        assert (
            "x    100." in text
            and "y     80." in text
            and "counts" in text
            and " Å" in text
        )
    finally:
        root.destroy()


def test_an_indexed_files_spots_are_their_own_layer(strong):
    run = subprocess.run(
        [
            MXI,
            "convert",
            str(strong / "indexed.expt"),
            str(strong / "indexed.refl"),
            "-o",
            str(strong / "indexed.rflx"),
        ],
        capture_output=True,
        text=True,
        cwd=strong,
    )
    assert run.returncode == 0, run.stderr
    root, v = open_viewer(strong / "indexed.rflx")
    try:
        layers = v.layers[0]
        assert layers.get("indexed") is not None and layers.get("unindexed") is None
        boxes = [c.cget("text") for c in v.toggles.winfo_children()]
        assert boxes == [f"Indexed spots ({len(layers['indexed'])})"]
        v.image.set(3)
        v.show()
        settle(root)
        assert v.canvas.find_withtag("indexed") and not v.canvas.find_withtag(
            "unindexed"
        )
        v.showing["indexed"].set(False)
        v.draw()
        assert not v.canvas.find_withtag("indexed"), "a layer hidden"
    finally:
        root.destroy()


def test_view_images_opens_the_steps_file_or_what_it_reads(strong):
    from mxi_ui.window import MainWindow

    root = tkinter.Tk()
    root.withdraw()
    w = MainWindow(root, str(strong))
    try:
        w.tree.selection_set("find")
        root.update()
        assert w.images_file().endswith("strong.rflx")
        viewer = w.view_images()
        assert viewer is not None and "Find spots" in viewer.top.title()
        w.tree.selection_set("index")
        root.update()
        assert w.images_file().endswith("strong.rflx"), "not run: what it reads"
    finally:
        root.destroy()


def test_view_images_with_nothing_to_view_says_so(tmp_path, monkeypatch):
    from tkinter import messagebox

    from mxi_ui.window import MainWindow

    said = []
    monkeypatch.setattr(messagebox, "showinfo", lambda *a: said.append(a))
    root = tkinter.Tk()
    root.withdraw()
    w = MainWindow(root, str(tmp_path))
    try:
        w.tree.selection_set("find")
        root.update()
        assert w.images_file() is None and w.view_images() is None
        assert said and "run Import first" in said[0][1]
    finally:
        root.destroy()


def test_contrast_and_colours_recolour_without_reducing_again(strong, monkeypatch):
    from mxi_ui import render

    root, v = open_viewer(strong / "strong.rflx")
    try:
        calls = []
        real = render.reduce
        monkeypatch.setattr(render, "reduce", lambda *a: calls.append(1) or real(*a))
        v.ceiling.set(v.ceiling.get() / 2)
        v.draw()
        v.scheme.set("heat")
        v.draw()
        assert calls == [], "the reduction kept"
        v.view.zoom(1.25, 100, 100)
        v.draw()
        assert calls == [1], "a new view reduced afresh"
    finally:
        root.destroy()


def test_zoom_without_a_wheel(strong):
    # Graeme has no wheel: buttons and keys zoom, a level each.
    root, v = open_viewer(strong / "strong.rflx")
    try:
        v.top.deiconify()
        settle(root)
        start = v.view.scale
        w, h = v._canvas_size()
        middle = v.view.to_detector(w / 2, h / 2)
        buttons = {
            b.cget("text"): b
            for b in v.top.winfo_children()[0].winfo_children()
            if b.winfo_class() == "TButton"
        }
        buttons["+"].invoke()
        assert v.view.scale < start
        assert all(
            abs(a - b) < 1e-6 for a, b in zip(v.view.to_canvas(*middle), (w / 2, h / 2))
        ), "about the middle"
        buttons["−"].invoke()
        assert v.view.scale == start
        v.top.focus_force()
        settle(root)
        v.top.event_generate("<KeyPress-plus>")
        v.top.event_generate("<KeyPress-equal>")
        settle(root)
        assert v.view.scale < start
        v.top.event_generate("<KeyPress-0>")
        settle(root)
        assert v.view.scale == start, "0 fits"
    finally:
        root.destroy()


def test_a_trackpads_burst_of_scrolls_zooms_a_level_a_tenth_of_a_second(strong):
    from mxi_ui.viewer import WHEEL_INTERVAL

    root, v = open_viewer(strong / "strong.rflx")
    try:
        start = v.view.scale

        class E:
            x, y = 100, 100

        taken = [
            v._wheel(E, 1, now=t)
            for t in (10.0, 10.01, 10.02, 10.05, 10.0 + WHEEL_INTERVAL + 0.01)
        ]
        assert taken == [True, False, False, False, True]
        levels = 0
        s = start
        from mxi_ui import render

        while s > v.view.scale + 1e-9:
            s = render.next_level(s, inwards=True)
            levels += 1
        assert levels == 2, "two levels, not five"
    finally:
        root.destroy()
