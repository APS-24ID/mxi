"""mxi_find on a series it cannot read: said clearly, and nothing written.

Graeme's two-lattice data: the master's /entry/data/data a virtual dataset
whose sources are in the master itself, ".", at /entry/data/data_000001 -- an
external link to the data file's /data. The upload renamed the data file,
adding "_1" to a name already taken, and mxi_find printed HDF5's error stack
six hundred times, said "No reflections found", and wrote an empty table:
images without spots, to anyone reading it so. Now each failure says the link
and the file that is not there, three at most, and a series none of whose
frames can be read writes nothing and says why.

The planted series of mxeq.fixtures.nxmx, its master rewritten in that shape:
read through both stages, the same spots as through its own master; one data
file renamed, the run failing and naming it; both, nothing written. Needs
MXI_FIND, h5py and hdf5plugin.
"""

import os
import subprocess
import sys

import pytest

from mxeq import refl

h5py = pytest.importorskip("h5py")
pytest.importorskip("hdf5plugin")
FIND = os.environ.get("MXI_FIND")
pytestmark = pytest.mark.skipif(not FIND, reason="set MXI_FIND")


def planted(tmp_path):
    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(tmp_path), "6"],
        check=True,
        capture_output=True,
    )
    return tmp_path / "series.nxs"


def through_links(master):
    """The master's data as Graeme's: a virtual dataset over the master
    itself, at external links to each data file's /data."""
    with h5py.File(master, "r+") as f:
        old = f["/entry/data/data"]
        shape, dtype = old.shape, old.dtype
        sources = [(s.file_name, s.dset_name) for s in old.virtual_sources()]
        del f["/entry/data/data"]
        layout = h5py.VirtualLayout(shape=shape, dtype=dtype)
        at = 0
        for k, (name, path) in enumerate(sources, start=1):
            with h5py.File(master.parent / name, "r") as g:
                n = g[path].shape[0]
            link = f"/entry/data/data_{k:06d}"
            f[link] = h5py.ExternalLink(name, "/" + path.lstrip("/"))
            layout[at : at + n] = h5py.VirtualSource(".", link, shape=(n,) + shape[1:])
            at += n
        f["/entry/data"].create_virtual_dataset("data", layout)


def find(directory, master, out):
    return subprocess.run(
        [FIND, str(master), "-o", out], capture_output=True, text=True, cwd=directory
    )


def test_frames_are_read_through_a_link_from_the_master_itself(tmp_path):
    master = planted(tmp_path)
    first = find(tmp_path, master, "own.refl")
    assert first.returncode == 0, first.stderr
    through_links(master)
    second = find(tmp_path, master, "linked.refl")
    assert second.returncode == 0, second.stderr
    a, b = refl.load(str(tmp_path / "own.refl")), refl.load(
        str(tmp_path / "linked.refl")
    )
    assert a.nrows == b.nrows > 0


def test_a_missing_data_file_is_named(tmp_path):
    master = planted(tmp_path)
    through_links(master)
    os.rename(tmp_path / "series_0001.h5", tmp_path / "series_0001_1.h5")
    r = find(tmp_path, master, "strong.refl")
    assert r.returncode != 0
    assert "series_0001.h5" in r.stderr and "there is no" in r.stderr
    assert r.stderr.count("HDF5-DIAG") == 0  # HDF5's own stack not printed


def test_a_series_none_of_whose_frames_can_be_read_writes_nothing(tmp_path):
    master = planted(tmp_path)
    through_links(master)
    for name in ("series_0000.h5", "series_0001.h5"):
        os.rename(tmp_path / name, tmp_path / name.replace(".h5", "_1.h5"))
    r = find(tmp_path, master, "strong.refl")
    assert r.returncode != 0
    assert "none of the 6 images could be read" in r.stderr
    assert "series_0000.h5" in r.stderr
    assert not (tmp_path / "strong.refl").exists()
