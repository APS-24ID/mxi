"""Frames are read with pread on one descriptor a data file, shared by every
reader thread: a descriptor a reader took threads x files of them, 16 x 18 past
macOS's 256 a process, and HDF5 could open nothing more (Graeme). Integration
with more threads than the limit allows descriptors must still run."""

import os
import resource
import subprocess

import pytest

BINARY = os.environ.get("MXI_INTEGRATE")
EXPT = os.environ.get("MXI_POSTREFINE_EXPT")
REFL = os.environ.get("MXI_POSTREFINE_REFL")
needs_data = pytest.mark.skipif(
    not (BINARY and EXPT and REFL),
    reason="set MXI_INTEGRATE, MXI_POSTREFINE_EXPT, MXI_POSTREFINE_REFL",
)


@needs_data
def test_sixteen_threads_need_far_fewer_than_twenty_descriptors(tmp_path):
    def limit():
        resource.setrlimit(resource.RLIMIT_NOFILE, (20, 20))

    r = subprocess.run(
        [
            BINARY,
            EXPT,
            REFL,
            "--threads",
            "16",
            "--last-image",
            "30",
            "-o",
            str(tmp_path / "out.rflx"),
        ],
        capture_output=True,
        text=True,
        preexec_fn=limit,
    )
    assert r.returncode == 0, r.stdout[-600:] + r.stderr
