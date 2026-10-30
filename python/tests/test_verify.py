"""mxi verify end to end: a planted series, its spots indexed to a lattice;
the report a lattice at a time, the images' mask read for what cannot be seen
(or --no-mask, said so), and a file with no lattice refused. The counting
itself is tested in C++ (tests/test_derivatives.cc): every prediction seen,
half, a mask, a supercell's eighth."""

import json
import os
import subprocess
import sys

import numpy as np
import pytest

from mxeq import refl, rflx

BUILD = os.path.dirname(os.environ.get("MXI_FIND", "")) or os.environ.get(
    "MXI_BUILD", ""
)
MXI = os.path.join(BUILD, "mxi")

pytestmark = pytest.mark.skipif(not os.path.exists(MXI), reason="mxi is not built")


def indexed(tmp_path):
    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(tmp_path), "6"],
        check=True,
        capture_output=True,
    )
    from mxeq.fixtures import nxmx

    nxmx.write_expt(str(tmp_path / "one.expt"), 6)
    doc = json.load(open(tmp_path / "one.expt"))
    # The fixture's list is what spot finding needs; prediction needs the whole
    # geometry: a beam, a goniometer, the panel's axes and origin, the scan's
    # oscillation -- six images of 1 degree, 60 mm.
    first = dict(doc["experiment"][0], imageset=0, beam=0, goniometer=0)
    doc["experiment"] = [first]
    doc["imageset"] = [
        {"__id__": "ImageSequence", "template": str(tmp_path / "series.nxs")}
    ]
    doc["beam"] = [
        {
            "__id__": "monochromatic",
            "direction": [0.0, 0.0, 1.0],
            "wavelength": 0.97,
            "divergence": 0.0,
            "sigma_divergence": 0.0,
            "polarization_normal": [0.0, 1.0, 0.0],
            "polarization_fraction": 0.999,
            "flux": 0.0,
            "transmission": 1.0,
        }
    ]
    doc["detector"][0]["panels"][0].update(
        {
            "fast_axis": [1.0, 0.0, 0.0],
            "slow_axis": [0.0, -1.0, 0.0],
            "origin": [-nxmx.WIDTH * 0.075 / 2, nxmx.HEIGHT * 0.075 / 2, -60.0],
            "thickness": 0.0,
            "material": "",
            "mu": 0.0,
            "gain": 1.0,
            "pedestal": 0.0,
            "identifier": "",
            "raw_image_offset": [0, 0],
            "mask": [],
            "px_mm_strategy": {"type": "SimplePxMmStrategy"},
        }
    )
    doc["goniometer"] = [
        {
            "rotation_axis": [1.0, 0.0, 0.0],
            "fixed_rotation": [1, 0, 0, 0, 1, 0, 0, 0, 1],
            "setting_rotation": [1, 0, 0, 0, 1, 0, 0, 0, 1],
        }
    ]
    doc["scan"][0].update({"batch_offset": 0, "oscillation": [0.0, 1.0]})
    json.dump(doc, open(tmp_path / "list.expt", "w"))
    run = subprocess.run(
        [MXI, "find", str(tmp_path / "list.expt"), "-o", str(tmp_path / "strong.rflx")],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    lattice = dict(first, crystal=0, identifier="1c0ffee0-0000-4000-8000-0000000000a1")
    doc["experiment"].append(lattice)
    doc["crystal"] = [
        {
            "__id__": "crystal",
            "real_space_a": [40.0, 0.0, 0.0],
            "real_space_b": [0.0, 40.0, 0.0],
            "real_space_c": [0.0, 0.0, 40.0],
            "space_group_hall_symbol": " P 1",
        }
    ]
    json.dump(doc, open(tmp_path / "indexed.expt", "w"))
    spots = rflx.load_table(str(tmp_path / "strong.rflx"))
    n = spots.nrows
    spots.columns["id"] = np.ones(n, dtype=spots.columns["id"].dtype)
    spots.columns["miller_index"] = np.tile(
        np.array([[1, 2, 3]], dtype=np.int32), (n, 1)
    )
    spots.columns["flags"] = np.full(n, 4, dtype=np.uint64)  # indexed
    spots.types["miller_index"] = "cctbx::miller::index<>"
    spots.types["flags"] = "std::size_t"
    spots.identifiers = {1: lattice["identifier"]}
    refl.write(str(tmp_path / "indexed.refl"), spots)
    return tmp_path / "indexed.expt", tmp_path / "indexed.refl"


def test_a_lattice_reported_with_the_images_mask(tmp_path):
    expt, table = indexed(tmp_path)
    run = subprocess.run(
        [MXI, "verify", str(expt), str(table), "--d-min", "3"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    out = run.stdout
    assert "Experiment 1:" in out and "predicted to 3.00 A (--d-min)" in out
    assert "unobservable" in out and "\n  all " in out
    assert "not taken into account" not in out, "the images' mask was read"
    assert "Experiment 0" not in out, "the crystal-less sweep is not a lattice"
    log = (tmp_path / "mxi_verify.log").read_text()
    assert "Experiment 1:" in log, "its table kept in its log"
    no_mask = subprocess.run(
        [MXI, "verify", str(expt), str(table), "--d-min", "3", "--no-mask"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert "masked pixels not taken into account: --no-mask" in no_mask.stdout


def test_a_file_with_no_lattice_is_refused(tmp_path):
    indexed(tmp_path)
    run = subprocess.run(
        [MXI, "verify", str(tmp_path / "strong.rflx")],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 1 and "has no lattice: index it first" in run.stderr
