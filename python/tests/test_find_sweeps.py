"""mxi_find on several sweeps: each experiment's spots from its own images,
into one table, every spot with its experiment's index as its id and each
identifier in the table's map -- as dials.find_spots writes several sweeps.

A planted series of twelve images, split in two -- images 1 to 6 and 7 to 12,
each experiment naming its own frames -- found alone and together. The table of
both must be the two halves exactly, column by column, shoeboxes included, the
first half's rows with id 0 and the second's with 1. Needs MXI_FIND.
"""

import json
import os
import subprocess
import sys

import numpy as np
import pytest

FIND = os.environ.get("MXI_FIND")
needs = pytest.mark.skipif(not FIND, reason="set MXI_FIND")


def lists(tmp_path):
    pytest.importorskip("h5py")
    pytest.importorskip("hdf5plugin")
    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(tmp_path / "s"), "12"],
        check=True,
        capture_output=True,
    )
    master = str(tmp_path / "s" / "series.nxs")
    base = json.load(open(tmp_path / "s" / "series.expt"))

    def half(first, last, identifier):
        scan = {"__id__": "Scan", "image_range": [first, last]}
        imageset = {
            "__id__": "ImageSequence",
            "template": master,
            "single_file_indices": list(range(first - 1, last)),
        }
        return scan, imageset, identifier

    halves = [
        half(1, 6, "a0000000-0000-4000-8000-000000000001"),
        half(7, 12, "b0000000-0000-4000-8000-000000000002"),
    ]

    def document(parts):
        return {
            "__id__": "ExperimentList",
            "experiment": [
                {
                    "__id__": "Experiment",
                    "identifier": ident,
                    "detector": 0,
                    "scan": i,
                    "imageset": i,
                }
                for i, (_, _, ident) in enumerate(parts)
            ],
            "detector": base["detector"],
            "scan": [p[0] for p in parts],
            "imageset": [p[1] for p in parts],
        }

    for name, parts in (("a", halves[:1]), ("b", halves[1:]), ("two", halves)):
        json.dump(document(parts), open(tmp_path / f"{name}.expt", "w"))


def find(tmp_path, name):
    r = subprocess.run(
        [
            FIND,
            str(tmp_path / f"{name}.expt"),
            "-o",
            str(tmp_path / f"{name}.refl"),
            "-j",
            "2",
        ],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stdout + r.stderr
    msgpack = pytest.importorskip("msgpack")
    with open(tmp_path / f"{name}.refl", "rb") as f:
        return msgpack.unpackb(f.read(), raw=False, strict_map_key=False)[2], r.stdout


@needs
def test_two_sweeps_are_the_two_halves_exactly(tmp_path):
    lists(tmp_path)
    a, _ = find(tmp_path, "a")
    b, _ = find(tmp_path, "b")
    two, out = find(tmp_path, "two")
    na, nb = a["nrows"], b["nrows"]
    assert na > 0 and nb > 0
    assert two["nrows"] == na + nb
    assert two["identifiers"] == {
        0: "a0000000-0000-4000-8000-000000000001",
        1: "b0000000-0000-4000-8000-000000000002",
    }
    for name, (_, (_, blob)) in two["data"].items():
        if name == "id":
            ids = np.frombuffer(blob, np.int32)
            assert (ids[:na] == 0).all() and (ids[na:] == 1).all()
            continue
        assert blob == a["data"][name][1][1] + b["data"][name][1][1], name
    assert (
        f"experiment 0: {na} reflections" in out
        and f"experiment 1: {nb} reflections" in out
    )


@needs
def test_a_master_named_with_several_sweeps_is_refused(tmp_path):
    lists(tmp_path)
    r = subprocess.run(
        [
            FIND,
            str(tmp_path / "two.expt"),
            "-x",
            str(tmp_path / "s" / "series.nxs"),
            "-o",
            str(tmp_path / "x.refl"),
        ],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 2 and "cannot belong to them all" in r.stderr
