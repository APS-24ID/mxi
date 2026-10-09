"""The older DECTRIS file writer's "nearly NeXus" master, as dxtbx reads it.

No NXdetector_module, no depends_on, no /entry/data/data: the detector's
geometry in geometry/orientation and geometry/translation, the scan as omega's
range_average, the frames in /entry/data/data_000001 and on, each a link to a
data file. mxi_import makes its models as dxtbx's FormatHDF5EigerNearlyNexus
does -- the orientation's and translation's z negated, omega about (-1, 0, 0)
from 0 in steps of omega_range_average rounded to 0.01, the recorded start not
used -- and the frame reader follows the links in order. Graeme's case: a 2023
Eiger 16M master, firmware 1.6.6.

A planted series of two data files of three 64x64 frames, a spot on frame 5:
imported, every model as dxtbx's rules give it; found, the spot on the frame
it is on, across the two files; and a data file missing, refused by name.
Needs MXI_IMPORT, MXI_FIND, h5py and hdf5plugin.
"""

import json
import os
import subprocess

import numpy as np
import pytest

from mxeq import refl

h5py = pytest.importorskip("h5py")
hdf5plugin = pytest.importorskip("hdf5plugin")
IMPORT = os.environ.get("MXI_IMPORT")
FIND = os.environ.get("MXI_FIND")
SIZE, PER_FILE = 64, 3
# dxtbx's own example of these fields: x along -x, y along -y, the detector
# 0.3 m away with the beam centre at (20.5, 30.5) pixels of 75 microns.
ORIENTATION = [-1.0, 0.0, 0.0, 0.0, -1.0, 0.0]
TRANSLATION = [20.5 * 75e-6, 30.5 * 75e-6, -0.3]


def plant(directory, files=2, omega_step=0.1000000015):
    rng = np.random.default_rng(0)
    frame_number = 0
    for k in range(1, files + 1):
        with h5py.File(directory / f"run_data_{k:06d}.h5", "w") as f:
            d = f.create_dataset(
                "entry/data/data",
                shape=(PER_FILE, SIZE, SIZE),
                dtype="u2",
                chunks=(1, SIZE, SIZE),
                **hdf5plugin.Bitshuffle(nelems=0, cname="lz4"),
            )
            for i in range(PER_FILE):
                frame = rng.poisson(0.1, (SIZE, SIZE)).astype("u2")
                if frame_number == 4:  # frame 5, the second file's second
                    frame[30:33, 40:43] += 400
                d[i] = frame
                frame_number += 1
    master = directory / "run_master.h5"
    with h5py.File(master, "w") as f:
        entry = f.create_group("entry")
        entry.attrs["NX_class"] = "NXentry"
        data = entry.create_group("data")
        data.attrs["NX_class"] = "NXdata"
        for k in range(1, files + 1):
            data[f"data_{k:06d}"] = h5py.ExternalLink(
                f"run_data_{k:06d}.h5", "/entry/data/data"
            )
        instrument = entry.create_group("instrument")
        beam = instrument.create_group("beam")
        beam.create_dataset("incident_wavelength", data=np.float32(0.9785649)).attrs[
            "units"
        ] = "angstrom"
        det = instrument.create_group("detector")
        det.attrs["NX_class"] = "NXdetector"
        det.create_dataset("description", data=b"Dectris Eiger 16M")
        det.create_dataset("detector_number", data=b"E-32-0116")
        for name, value in (
            ("x_pixel_size", 75e-6),
            ("y_pixel_size", 75e-6),
            ("sensor_thickness", 450e-6),
        ):
            det.create_dataset(name, data=np.float32(value)).attrs["units"] = "m"
        det.create_dataset("sensor_material", data=b"Si")
        det.create_dataset("bit_depth_image", data=np.int32(16))
        det.create_dataset("count_time", data=np.float32(0.01)).attrs["units"] = "s"
        det.create_dataset("geometry/orientation/value", data=np.array(ORIENTATION))
        det.create_dataset("geometry/translation/distances", data=np.array(TRANSLATION))
        spec = det.create_group("detectorSpecific")
        spec.create_dataset("countrate_correction_count_cutoff", data=np.uint32(25000))
        spec.create_dataset("nimages", data=np.uint32(PER_FILE * files))
        spec.create_dataset("ntrigger", data=np.uint32(1))
        spec.create_dataset("x_pixels_in_detector", data=np.uint32(SIZE))
        spec.create_dataset("y_pixels_in_detector", data=np.uint32(SIZE))
        gon = entry.create_group("sample").create_group("goniometer")
        gon.create_dataset("omega_range_average", data=np.float32(omega_step))
        gon.create_dataset("omega_start", data=np.float32(-40.0))
    return master


def run(directory, *words):
    return subprocess.run(list(words), capture_output=True, text=True, cwd=directory)


@pytest.mark.skipif(not IMPORT, reason="set MXI_IMPORT")
def test_the_models_are_dxtbx_s(tmp_path):
    master = plant(tmp_path)
    r = run(tmp_path, IMPORT, str(master), "-o", "imported.expt")
    assert r.returncode == 0, r.stderr + r.stdout
    assert "nearly-NeXus" in r.stdout and "omega_start, -40" in r.stdout
    e = json.load(open(tmp_path / "imported.expt"))
    p = e["detector"][0]["panels"][0]
    # (x, y, -z) of the translation, in mm, into DIALS's frame: (-x, y, -z).
    assert p["origin"] == pytest.approx([-20.5 * 0.075, 30.5 * 0.075, -300.0], abs=1e-6)
    assert p["fast_axis"] == pytest.approx([1.0, 0.0, 0.0])
    assert p["slow_axis"] == pytest.approx([0.0, -1.0, 0.0])
    assert p["image_size"] == [SIZE, SIZE]
    assert p["pixel_size"] == pytest.approx([0.075, 0.075], rel=1e-6)
    assert p["trusted_range"][1] == 25000
    g = e["goniometer"][0]
    assert (
        g["axes"] == [[1.0, 0.0, 0.0]]
        and g["names"] == ["omega"]
        and g["scan_axis"] == 0
    )
    s = e["scan"][0]
    assert s["image_range"] == [1, 2 * PER_FILE]  # the linked files' frames
    assert s["properties"]["oscillation"] == pytest.approx(
        [0.1 * i for i in range(2 * PER_FILE)]
    )


@pytest.mark.skipif(not FIND, reason="set MXI_FIND")
def test_the_frames_are_read_across_the_linked_files(tmp_path):
    master = plant(tmp_path)
    r = run(tmp_path, FIND, str(master), "-o", "strong.refl")
    assert r.returncode == 0, r.stderr + r.stdout
    assert f"{2 * PER_FILE} images of {SIZE} x {SIZE}" in r.stdout
    t = refl.load(str(tmp_path / "strong.refl"))
    xyz = np.asarray(t.columns["xyzobs.px.value"], float).reshape(-1, 3)
    planted = (np.abs(xyz[:, 0] - 41.5) < 1.5) & (np.abs(xyz[:, 1] - 31.5) < 1.5)
    assert planted.sum() == 1
    assert 4.0 <= xyz[planted][0, 2] < 5.0  # frame 5: z from 4 to 5


@pytest.mark.skipif(not FIND, reason="set MXI_FIND")
def test_a_missing_data_file_is_named(tmp_path):
    master = plant(tmp_path)
    (tmp_path / "run_data_000002.h5").unlink()
    r = run(tmp_path, FIND, str(master), "-o", "strong.refl")
    assert r.returncode != 0
    assert "data_000002" in r.stderr and "not there" in r.stderr
