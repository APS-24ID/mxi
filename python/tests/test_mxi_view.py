"""mxi_view and mxi_tools.py, mxi's Python tools.

mxi_tools: its reflection reader against mxeq's; frames summed against the
file read directly; markers on exactly the images of their boxes; masters of
several blocks; mxi found as an installed release lays it out. mxi view,
mxi ui's viewer since 1.1.5 (its own tests, test_ui_viewer.py), run from the
PATH. napari's viewer and its tests are gone.
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
MXI = os.path.join(BUILD, "mxi")

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
    assert len(images) == 6 and images.frame_shape == raw.shape
    assert np.allclose(images.summed(0), raw)


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


def test_a_spot_is_marked_on_exactly_the_images_of_its_box():
    # Its box's images, z0 to z1 - 1: the spot finder's are its pixels'. Not
    # from the centroid's variance -- its uncertainty, not its spread -- with
    # an image either side at least: a spot on one image was marked on its
    # neighbours, "highlighted on the next image" (Graeme, 26 October).
    t = tools()
    table = {
        "id": np.array([0, 0]),
        "xyzobs.px.value": np.array([[10.5, 20.5, 10.4], [5.5, 5.5, 3.5]]),
        "xyzobs.px.variance": np.array([[0.1, 0.1, 1.0], [0.1, 0.1, 0.08]]),
        "bbox": np.array([[9, 12, 19, 22, 8, 13], [4, 7, 4, 7, 3, 4]]),
    }
    m = t.Markers(table, "xyzobs.px.value", [0])

    def images(row_y):
        return [f for f in range(0, 20) if any(m.on(f)[0][:, 0] == row_y)]

    assert images(20.0) == [8, 9, 10, 11, 12]  # its box's images
    assert images(5.0) == [3]  # one image: not its neighbours
    yx, own = m.on(10)
    assert yx.tolist() == [[20.0, 10.0]] and own.tolist() == [True]  # its own: red
    assert m.on(9)[1].tolist() == [False]
    # Three images summed for each shown: on 6, its window 6-8 reaches 8.
    assert any(m.on(6, window=3)[0][:, 0] == 20.0) and not any(
        m.on(5, window=3)[0][:, 0] == 20.0
    )
    # With no box, its own image alone.
    bare = {"id": np.array([0]), "xyzcal.px": np.array([[1.5, 1.5, 4.5]])}
    p = t.Markers(bare, "xyzcal.px", [0])
    assert [f for f in range(10) if len(p.on(f)[0])] == [4]


def test_a_master_with_data_and_its_blocks_counts_each_frame_once(tmp_path):
    # Graeme's: a virtual "data" over the blocks data_000001 ... -- both counted
    # had shown every image twice, a slider to 7199 for 3600.
    import h5py

    t = tools()
    frames = np.arange(3 * 4 * 5, dtype=np.uint16).reshape(3, 4, 5)
    with h5py.File(tmp_path / "m.nxs", "w") as f:
        f.create_dataset("entry/data/data_000001", data=frames)
        f.create_dataset("entry/data/data", data=frames)
    images = t.Images(str(tmp_path / "m.nxs"))
    assert images.frames == 3 and len(images) == 3 and images.frame_shape == (4, 5)


def test_an_installed_layout_runs_mxis_commands(tmp_path, monkeypatch):
    # Since 1.1 an installed bin has mxi, mxi_view and mxi_tools.py, nothing
    # else: the tools run mxi untrusted, not a program mxi_untrusted, which
    # is in a build directory alone (as a link).
    import shutil

    strong, master = planted(tmp_path)
    installed = tmp_path / "bin"
    installed.mkdir()
    shutil.copy(os.path.realpath(os.path.join(BUILD, "mxi")), installed / "mxi")
    shutil.copy(TOOLS, installed / "mxi_tools.py")
    monkeypatch.setenv("PATH", "/usr/bin:/bin")
    t = load("mxi_tools_installed", str(installed / "mxi_tools.py"))
    assert t.mxi_command("untrusted") == [str(installed / "mxi"), "untrusted"]
    monkeypatch.chdir(tmp_path)
    t.save_untrusted(
        str(strong), {0: [{"panel": 0, "rectangle": [10.0, 10.0, 50.0, 40.0]}]}
    )
    written = rflx.experiments_dict(str(tmp_path / "masked.rflx"))["imageset"][0]
    assert written["untrusted"] == [{"panel": 0, "rectangle": [10.0, 10.0, 50.0, 40.0]}]


def test_mxi_view_runs_from_the_path_as_plain_mxi(tmp_path):
    # mxi finds the viewer beside itself: where it is, the system says --
    # /proc/self/exe, _NSGetExecutablePath -- not argv[0], which for plain
    # "mxi" from the PATH has no directory; on macOS that had sent it to
    # "./mxi_view" ("cannot run ./mxi_view", Graeme, 3 November). From
    # another directory, and through a link to it elsewhere on the PATH.
    import shutil

    installed = tmp_path / "bin"
    installed.mkdir()
    shutil.copy(os.path.realpath(os.path.join(BUILD, "mxi")), installed / "mxi")
    shutil.copy(VIEW, installed / "mxi_view")
    shutil.copy(TOOLS, installed / "mxi_tools.py")
    os.chmod(installed / "mxi_view", 0o755)
    links = tmp_path / "links"
    links.mkdir()
    os.symlink(installed / "mxi", links / "mxi")
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    for first in (installed, links):
        env = dict(os.environ, PATH=f"{first}{os.pathsep}{os.environ['PATH']}")
        run = subprocess.run(
            ["mxi", "view", "--version"],
            capture_output=True,
            text=True,
            cwd=elsewhere,
            env=env,
        )
        assert run.returncode == 0, (first.name, run.stderr)
        assert run.stdout.startswith("mxi view "), run.stdout


def test_mxi_view_opens_mxi_uis_viewer(tmp_path):
    # mxi view is mxi ui's viewer since napari's was retired: a file named
    # that is not there said first, whatever the display; with one, a window
    # opened and kept open, the options -o, --sum, --sweep and --mask taken.
    run = subprocess.run(
        [MXI, "view", "nowhere.rflx"], capture_output=True, text=True, cwd=tmp_path
    )
    assert run.returncode == 2 and "no such file: nowhere.rflx" in run.stderr
    helped = subprocess.run([MXI, "view", "--help"], capture_output=True, text=True)
    assert all(o in helped.stdout for o in ("-o", "--sum", "--sweep", "--mask"))
    if not os.environ.get("DISPLAY") or os.environ.get("MXI_NO_WINDOW_TESTS"):
        pytest.skip("no display for the window (or the release builds')")
    strong, _master = planted(tmp_path)
    import time

    p = subprocess.Popen(
        [MXI, "view", str(strong), "--sum", "3", "-o", "mine.rflx", "--mask"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=tmp_path,
    )
    time.sleep(4)
    alive = p.poll() is None
    p.terminate()
    _out, err = p.communicate(timeout=10)
    assert alive, err.decode()[-500:]
