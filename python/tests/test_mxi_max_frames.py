"""mxi_max reads the frames the series has, by its own keys, as mxi_find
does -- not 0 to the data's declared extent. A converted master's virtual
dataset can declare more frames than its mappings cover; mxi_max counted to
the extent and stopped at the first frame in no file, where mxi_find, taking
the series' keys, read them all (Graeme, converted CBF data, 15 October)."""

import os
import subprocess
import sys

import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("hdf5plugin")
MAX = os.environ.get("MXI_MAX")
pytestmark = pytest.mark.skipif(not MAX, reason="set MXI_MAX")


def planted(tmp_path, n=6):
    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(tmp_path), str(n)],
        check=True,
        capture_output=True,
    )
    return tmp_path / "series.nxs"


def declare_more(master, extra):
    """The master's virtual dataset made `extra` frames longer than its
    mappings: those frames in no file."""
    with h5py.File(master, "r+") as f:
        old = f["/entry/data/data"]
        shape, dtype = old.shape, old.dtype
        sources = [
            (s.file_name, s.dset_name, s.vspace.get_select_bounds())
            for s in old.virtual_sources()
        ]
        del f["/entry/data/data"]
        layout = h5py.VirtualLayout(shape=(shape[0] + extra,) + shape[1:], dtype=dtype)
        for name, path, (lo, hi) in sources:
            n = hi[0] - lo[0] + 1
            layout[lo[0] : lo[0] + n] = h5py.VirtualSource(
                name, path, shape=(n,) + shape[1:]
            )
        f["/entry/data"].create_virtual_dataset("data", layout, fillvalue=0)


def run(master):
    return subprocess.run(
        [MAX, str(master)], capture_output=True, text=True, cwd=master.parent
    )


def test_the_frames_mapped_are_read(tmp_path):
    master = planted(tmp_path)
    plain = run(master)
    assert plain.returncode in (0, 1), plain.stderr
    declare_more(master, 4)
    with h5py.File(master, "r") as f:
        assert f["/entry/data/data"].shape[0] == 10
    longer = run(master)
    assert longer.returncode in (0, 1), longer.stderr
    # The same six frames, the same answer.
    assert "6 frames" in longer.stdout, longer.stdout
    assert plain.stdout.split("read in")[0] == longer.stdout.split("read in")[0]
    assert plain.stdout.split("the largest")[1] == longer.stdout.split("the largest")[1]
