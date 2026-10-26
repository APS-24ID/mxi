"""mxi_mask, the window driven by synthetic mouse events, with no display.

Two sweeps of the planted series: a rectangle dragged and a polygon clicked
and closed on the first; Undo; Same as last on the second, which takes the
first's; Clear and Undo; Save, which writes both image sets' regions into the
file and the last mask beside it, mxi_mask.json.
"""

import importlib.machinery
import importlib.util
import json
import os
import subprocess
import sys
import types

import pytest

from mxeq import rflx

BUILD = os.path.dirname(os.environ.get("MXI_FIND", "")) or os.environ.get(
    "MXI_BUILD", ""
)
MASK = os.path.join(BUILD, "mxi_mask")

pytestmark = pytest.mark.skipif(
    not (os.path.exists(MASK) and os.path.exists(os.path.join(BUILD, "mxi_untrusted"))),
    reason="mxi_mask and mxi_untrusted are not built",
)


def load_editor():
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    loader = importlib.machinery.SourceFileLoader("mxi_mask_module", MASK)
    spec = importlib.util.spec_from_loader("mxi_mask_module", loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def two_sweeps(tmp_path):
    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(tmp_path), "6"],
        check=True,
        capture_output=True,
    )
    from mxeq.fixtures import nxmx

    nxmx.write_expt(str(tmp_path / "one.expt"), 6)
    doc = json.load(open(tmp_path / "one.expt"))
    master = str(tmp_path / "series.nxs")
    first = dict(doc["experiment"][0], imageset=0)
    second = dict(first, imageset=1, identifier="1c0ffee0-0000-4000-8000-000000000002")
    doc["experiment"] = [first, second]
    doc["imageset"] = [
        {"__id__": "ImageSequence", "template": master},
        {"__id__": "ImageSequence", "template": master},
    ]
    json.dump(doc, open(tmp_path / "two.expt", "w"))
    subprocess.run(
        [
            os.path.join(BUILD, "mxi_convert"),
            str(tmp_path / "two.expt"),
            "-o",
            str(tmp_path / "two.rflx"),
        ],
        check=True,
        capture_output=True,
        cwd=tmp_path,
    )
    return str(tmp_path / "two.rflx")


def event(editor, x, y, button=1, dblclick=False):
    return types.SimpleNamespace(
        inaxes=editor.ax, xdata=x, ydata=y, button=button, dblclick=dblclick, key=None
    )


def test_drawing_undo_same_as_last_and_save(tmp_path, monkeypatch):
    path = two_sweeps(tmp_path)  # before the move: PYTHONPATH may be relative
    monkeypatch.chdir(tmp_path)
    module = load_editor()
    editor = module.Editor(path, module.read_experiments(path), 6)
    assert len(editor.sweeps) == 2

    # A rectangle, dragged corner to corner.
    editor.on_press(event(editor, 100, 50))
    editor.on_motion(event(editor, 180, 90))
    editor.on_release(event(editor, 260, 170))
    assert editor.sweep["shapes"] == [
        {"panel": 0, "rectangle": [100.0, 50.0, 260.0, 170.0]}
    ]

    # A polygon: three clicks, closed with the right button.
    editor.choose("polygon")
    for x, y in ((20, 200), (80, 240), (20, 250)):
        editor.on_press(event(editor, x, y))
    editor.on_press(event(editor, 0, 0, button=3))
    assert editor.sweep["shapes"][-1]["polygon"] == [[20, 200], [80, 240], [20, 250]]

    # Undo takes the polygon away, and only it.
    editor.undo()
    assert len(editor.sweep["shapes"]) == 1

    # The second sweep: Same as last takes the first's.
    editor.go(+1)
    assert editor.current == 1 and editor.sweep["shapes"] == []
    editor.same_as_last()
    assert editor.sweep["shapes"] == editor.sweeps[0]["shapes"]
    editor.clear()
    assert editor.sweep["shapes"] == []
    editor.undo()
    assert len(editor.sweep["shapes"]) == 1

    # A circle on it too, then Save.
    editor.choose("circle")
    editor.on_press(event(editor, 330, 200))
    editor.on_release(event(editor, 355, 200))
    editor.save()
    # Into masked.rflx, the input as it was (Graeme, 22 October).
    assert "untrusted" not in rflx.experiments_dict(path)["imageset"][0]
    written = rflx.experiments_dict(str(tmp_path / "masked.rflx"))["imageset"]
    assert written[0]["untrusted"] == [
        {"panel": 0, "rectangle": [100.0, 50.0, 260.0, 170.0]}
    ]
    assert written[1]["untrusted"][-1] == {"panel": 0, "circle": [330.0, 200.0, 25.0]}
    last = json.load(open(tmp_path / "mxi_mask.json"))
    assert last["last"] == written[1]["untrusted"]
    assert set(last["imagesets"]) == {"0", "1"}
    editor.plt.close(editor.fig)


def test_a_saved_mask_applies_to_another_file_with_no_window(tmp_path, monkeypatch):
    # The last mask, kept: applied by mxi_untrusted to every image set.
    path = two_sweeps(tmp_path)
    monkeypatch.chdir(tmp_path)
    (tmp_path / "last.json").write_text(json.dumps([{"circle": [10, 10, 5]}]))
    run = subprocess.run(
        [
            os.path.join(BUILD, "mxi_untrusted"),
            path,
            "--set",
            str(tmp_path / "last.json"),
        ],
        capture_output=True,
        text=True,
    )
    assert run.returncode == 0, run.stderr
    sets = rflx.experiments_dict(str(tmp_path / "masked.rflx"))["imageset"]
    assert all("untrusted" not in s for s in rflx.experiments_dict(path)["imageset"])
    assert all(s["untrusted"] == [{"circle": [10, 10, 5]}] for s in sets)


def test_a_path_relative_to_the_working_directory_is_read(tmp_path, monkeypatch):
    # As it is typed: "mxi_mask imported.rflx" in the directory holding it.
    # Read through mxi_convert in a scratch directory, it had not been found
    # there (Graeme, 22 October).
    two_sweeps(tmp_path)
    monkeypatch.chdir(tmp_path)
    module = load_editor()
    document = module.read_experiments("two.rflx")
    assert len(document["imageset"]) == 2
    editor = module.Editor("two.rflx", document, 2)
    editor.on_press(event(editor, 10, 10))
    editor.on_release(event(editor, 30, 30))
    editor.save()
    assert (tmp_path / "masked.rflx").exists()
    editor.plt.close(editor.fig)
