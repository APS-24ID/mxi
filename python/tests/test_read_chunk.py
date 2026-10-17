"""read_chunk (src/nxmx_reader.cc) is H5Dread_chunk with the bytes read by
pread outside HDF5's lock; MXI_READ_CHUNK=hdf5 sends every read through
H5Dread_chunk itself. The two must find the same spots, to the bit."""

import os
import subprocess

import numpy as np
import pytest

from mxeq import refl

FIND = os.environ.get("MXI_FIND")
needs = pytest.mark.skipif(not FIND, reason="set MXI_FIND")


@needs
def test_the_direct_read_is_h5dread_chunks_bytes(tmp_path):
    pytest.importorskip("hdf5plugin")
    import test_pixel_mask

    master = test_pixel_mask.series(tmp_path, mask_hot=False)
    tables = {}
    for how, extra in (("direct", {}), ("hdf5", {"MXI_READ_CHUNK": "hdf5"})):
        out = tmp_path / f"{how}.refl"
        r = subprocess.run(
            [FIND, str(master), "-o", str(out)],
            capture_output=True,
            text=True,
            cwd=tmp_path,
            env=dict(os.environ, **extra),
        )
        assert r.returncode == 0, r.stdout[-400:] + r.stderr
        tables[how] = refl.load(str(out))
    a, b = tables["direct"], tables["hdf5"]
    assert set(a.columns) == set(b.columns)
    for name in a.columns:
        assert np.array_equal(
            np.asarray(a.columns[name]), np.asarray(b.columns[name])
        ), name
