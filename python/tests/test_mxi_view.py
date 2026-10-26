"""mxi_view and mxi_tools.py, mxi's Python tools.

The module without napari: untrusted regions to napari's shapes and back, in
the two pixel conventions; spots as napari points; its reflection reader
against mxeq's; frames summed against the file read directly. The window,
where napari and a display are there (xvfb-run in a sandbox): two sweeps of
the planted series, each spot point on a bright pixel, a rectangle drawn,
Same as last, Save into masked.rflx.
"""

import importlib.machinery
import importlib.util
import json
import os
import subprocess
import sys

import numpy as np
import pytest

from mxeq import rflx

BUILD = os.path.dirname(os.environ.get("MXI_FIND", "")) or os.environ.get(
    "MXI_BUILD", ""
)
TOOLS = os.path.join(BUILD, "mxi_tools.py")
VIEW = os.path.join(BUILD, "mxi_view")

pytestmark = pytest.mark.skipif(
    not (os.path.exists(TOOLS) and os.path.exists(VIEW)),
    reason="mxi_view and mxi_tools.py are not built",
)


def load(name, path):
    loader = importlib.machinery.SourceFileLoader(name, path)
    spec = importlib.util.spec_from_loader(name, loader)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


def tools():
    if BUILD not in sys.path:
        sys.path.insert(0, BUILD)
    return load("mxi_tools", TOOLS)


def test_regions_and_shapes_round_trip_in_both_conventions():
    t = tools()
    for region in (
        {"panel": 0, "rectangle": [100.0, 50.0, 260.0, 170.0]},
        {"panel": 0, "circle": [330.0, 200.0, 25.0]},
        {"panel": 0, "polygon": [[20.0, 200.0], [80.0, 240.0], [20.0, 250.0]]},
    ):
        kind, vertices = t.untrusted_to_shape(region)
        assert t.shape_to_untrusted(kind, vertices) == region
    # napari's pixel i is centred on i, DIALS's spans i to i + 1: a rectangle
    # on napari's corners -0.5 and 9.5 is DIALS's 0 to 10.
    r = t.shape_to_untrusted(
        "rectangle", [[-0.5, -0.5], [-0.5, 9.5], [4.5, 9.5], [4.5, -0.5]]
    )
    assert r["rectangle"] == [0.0, 0.0, 10.0, 5.0]
    # An ellipse not round becomes a polygon; a slice coordinate is dropped.
    e = t.shape_to_untrusted(
        "ellipse", [[3, 0, 0], [3, 0, 40], [3, 20, 40], [3, 20, 0]]
    )
    assert "polygon" in e and len(e["polygon"]) == 48
    assert t.shape_to_untrusted("line", [[0, 0], [5, 5]]) is None


def test_spots_become_napari_points_half_a_pixel_back():
    t = tools()
    table = {
        "id": np.array([0, 1, 0]),
        "xyzobs.px.value": np.array([[10.5, 20.5, 3.5], [1, 1, 1], [0.0, 0.0, 0.0]]),
    }
    p = t.points_of(table, "xyzobs.px.value", [0])
    assert p.tolist() == [[3.0, 20.0, 10.0]]  # (z, y, x); the zero row none


def planted(tmp_path, sets=1):
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
    doc["experiment"] = [first]
    doc["imageset"] = [{"__id__": "ImageSequence", "template": master}]
    if sets == 2:
        doc["experiment"].append(
            dict(first, imageset=1, identifier="1c0ffee0-0000-4000-8000-000000000002")
        )
        doc["imageset"].append({"__id__": "ImageSequence", "template": master})
    json.dump(doc, open(tmp_path / "list.expt", "w"))
    strong = tmp_path / "strong.rflx"
    run = subprocess.run(
        [
            os.path.join(BUILD, "mxi_find"),
            str(tmp_path / "list.expt"),
            "-o",
            str(strong),
        ],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    return strong, master


def test_reflections_and_frames_are_read_as_mxeq_and_h5py_read_them(tmp_path):
    t = tools()
    strong, master = planted(tmp_path)
    mine = t.read_reflections(str(strong), {"id", "xyzobs.px.value"})
    theirs = rflx.load_table(str(strong))
    assert np.array_equal(
        mine["xyzobs.px.value"], np.asarray(theirs["xyzobs.px.value"])
    )
    import h5py

    images = t.Images(master, window=3)
    with h5py.File(master, "r") as f:
        names = sorted(
            n for n in f["/entry/data"] if n == "data" or n.startswith("data_")
        )
        raw = np.sum([f["/entry/data"][names[0]][k] for k in range(3)], axis=0)
    assert images.shape == (6,) + raw.shape
    assert np.allclose(images[0], raw)
    assert images[0:2].shape == (2,) + raw.shape


needs_display = pytest.mark.skipif(
    not os.environ.get("DISPLAY") or importlib.util.find_spec("napari") is None,
    reason="napari and a display (xvfb-run) are needed for the window",
)


@needs_display
def test_the_window_shows_the_spots_and_saves_the_regions(tmp_path, monkeypatch):
    strong, master = planted(tmp_path, sets=2)
    monkeypatch.chdir(tmp_path)
    tools()
    view_module = load("mxi_view_module", VIEW)
    view = view_module.View(str(strong), "masked.rflx", 1)
    try:
        layers = [layer.name for layer in view.viewer.layers]
        assert {"images", "spots", "untrusted"} <= set(layers)
        # Each spot point on its spot: the brightest pixel near it, on its own
        # image, within a pixel -- the two conventions met (napari's pixel i
        # centred on i, DIALS's spanning i to i + 1). Measured: on it.
        points = view.viewer.layers["spots"].data
        images = view.viewer.layers["images"].data
        for z, y, x in points[:20]:
            frame = images[int(round(z))]
            yi, xi = int(round(y)), int(round(x))
            w = frame[yi - 4 : yi + 5, xi - 4 : xi + 5]
            j, i = np.unravel_index(np.argmax(w), w.shape)
            assert abs(yi - 4 + j - y) <= 1 and abs(xi - 4 + i - x) <= 1
        # A rectangle drawn on the first sweep; Same as last on the second.
        view.viewer.layers["untrusted"].add_rectangles(
            np.array([[49.5, 99.5], [169.5, 259.5]])
        )
        view.show(1)
        view.same_as_last()
        view.save()
        written = rflx.experiments_dict(str(tmp_path / "masked.rflx"))["imageset"]
        assert written[0]["untrusted"] == [
            {"panel": 0, "rectangle": [100.0, 50.0, 260.0, 170.0]}
        ]
        assert written[1]["untrusted"] == written[0]["untrusted"]
        assert "untrusted" not in rflx.experiments_dict(str(strong))["imageset"][0]
    finally:
        view.viewer.close()


def test_the_grey_scale_tops_out_well_above_the_background():
    # A diffraction image is nearly all background: the 99.5th percentile was
    # the background itself, and the whole image black (23 October).
    t = tools()
    frame = np.full((100, 100), 2.0)
    frame[10:13, 10:13] = 2500.0  # one spot
    top = t.display_top(frame)
    assert top >= 10.0 and top > 4 * np.median(frame)
    assert t.display_top(np.zeros((4, 4))) == 10.0


def test_compressed_columns_are_read_in_a_fresh_process(tmp_path):
    # mxi writes a numeric column of 256 elements or more bitshuffle-LZ4
    # (filter 32008). Read in a process that has imported nothing else -- here,
    # pytest had imported hdf5plugin through mxeq, and a test file too small
    # to be compressed hid it; Graeme's indexed.rflx did not.
    import h5py
    import hdf5plugin

    path = tmp_path / "compressed.rflx"
    xyz = np.arange(3000, dtype=float).reshape(1000, 3)
    with h5py.File(path, "w") as f:
        g = f.create_group("reflections/0")
        g.create_dataset(
            "xyzobs.px.value",
            data=xyz,
            chunks=(1000, 3),
            **hdf5plugin.Bitshuffle(cname="lz4"),
        )
        g.create_dataset("id", data=np.zeros(1000, np.int32))
    code = (
        f"import sys; sys.path.insert(0, {BUILD!r}); import mxi_tools; "
        f"t = mxi_tools.read_reflections({str(path)!r}); "
        "print(t['xyzobs.px.value'].shape, float(t['xyzobs.px.value'][-1, -1]))"
    )
    run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert run.stdout.strip() == "(1000, 3) 2999.0"
