"""mxi_import reads an NXmx master's geometry right, and its overrides do
what they say.

A master is planted with h5py -- metadata only, since mxi_import reads no
pixels -- awkwardly: lengths in metres, a detector chain with a two-theta
rotation in it, a fixed goniometer axis at an angle, a relative depends_on.
What the experiment list must hold is worked out here, independently, with
numpy's own rotations, and the overrides are each held to what they ask for.
And, where the data is to hand, the insulin master's mu, goniometer and scan
against dials.import's own output. Needs MXI_IMPORT.
"""

import json
import os
import subprocess

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

IMPORT = os.environ.get("MXI_IMPORT")
needs_import = pytest.mark.skipif(not IMPORT, reason="set MXI_IMPORT")


def rotation(axis, degrees):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    t = np.radians(degrees)
    k = np.array(
        [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]]
    )
    return np.eye(3) + np.sin(t) * k + (1 - np.cos(t)) * (k @ k)


def imgcif(v):
    return np.array([-v[0], v[1], -v[2]])


TWO_THETA = 15.0
DISTANCE_M = 0.25
OFFSET_M = 0.12
OFFSET_VEC = np.array([0.6, 0.8, 0.0])
CHI = 20.0
PIVOT_M = np.array([0.0, 0.005, 0.002])  # the two-theta arm's offset: not at the sample


def plant(
    path,
    chi=CHI,
    module_offset_offset=None,
    saturation=20000,
    meta=None,
    material=b"Silicon",
    wavelength=1.0,
):
    with h5py.File(path, "w") as f:
        e = f.create_group("entry")
        e.attrs["NX_class"] = "NXentry"
        beam = e.create_group("instrument/beam")
        beam.create_dataset("incident_wavelength", data=wavelength).attrs["units"] = (
            "angstrom"
        )
        det = e.create_group("instrument/detector")
        det.attrs["NX_class"] = "NXdetector"
        det.create_dataset("sensor_material", data=material)
        det.create_dataset("sensor_thickness", data=0.00045).attrs["units"] = "m"
        if meta is not None:
            # As DECTRIS masters give them: links into the _meta.h5 beside it,
            # which may be there or not.
            name = path.with_name(path.stem + "_meta.h5")
            if meta.get("write", True):
                with h5py.File(name, "w") as m:
                    g = m.create_group("_dectris")
                    if "cutoff" in meta:
                        g.create_dataset(
                            "countrate_correction_count_cutoff", data=meta["cutoff"]
                        )
                    if "bits" in meta:
                        g.create_dataset("bit_depth_image", data=meta["bits"])
            if "cutoff" in meta:
                det["saturation_value"] = h5py.ExternalLink(
                    name.name, "/_dectris/countrate_correction_count_cutoff"
                )
            if "bits" in meta:
                det["bit_depth_readout"] = h5py.ExternalLink(
                    name.name, "/_dectris/bit_depth_image"
                )
        elif saturation is not None:
            det.create_dataset("saturation_value", data=saturation)
        det.create_dataset("count_time", data=0.01)
        t = e.create_group("instrument/transformations")
        tt = t.create_dataset("two_theta", data=TWO_THETA)
        tt.attrs.update(
            {
                "transformation_type": "rotation",
                "vector": [-1.0, 0, 0],
                "units": "deg",
                "depends_on": ".",
                "offset": PIVOT_M,
                "offset_units": "m",
            }
        )
        dz = t.create_dataset("det_z", data=DISTANCE_M)
        dz.attrs.update(
            {
                "transformation_type": "translation",
                "vector": [0, 0, 1.0],
                "units": "m",
                "depends_on": "two_theta",
            }
        )
        mod = det.create_group("module")
        mod.attrs["NX_class"] = "NXdetector_module"
        mod.create_dataset("data_origin", data=[0, 0])
        mod.create_dataset("data_size", data=[300, 200])
        mo = mod.create_dataset("module_offset", data=OFFSET_M)
        mo.attrs.update(
            {
                "transformation_type": "translation",
                "vector": OFFSET_VEC,
                "units": "m",
                "depends_on": "/entry/instrument/transformations/det_z",
            }
        )
        if module_offset_offset is not None:
            # In metres, with no offset_units: as Diamond's Eiger masters write it.
            mo.attrs["offset"] = np.asarray(module_offset_offset, float)
        for name, vec in (
            ("fast_pixel_direction", [-1.0, 0, 0]),
            ("slow_pixel_direction", [0, -1.0, 0]),
        ):
            d = mod.create_dataset(name, data=7.5e-5)
            d.attrs.update(
                {
                    "transformation_type": "translation",
                    "vector": vec,
                    "units": "m",
                    "depends_on": "module_offset",
                }
            )
        s = e.create_group("sample")
        s.create_dataset("depends_on", data=b"/entry/sample/transformations/phi")
        g = s.create_group("transformations")
        om = g.create_dataset("omega", data=np.arange(10) * 0.5 + 10.0)
        om.attrs.update(
            {
                "transformation_type": "rotation",
                "vector": [-1.0, 0, 0],
                "units": "deg",
                "depends_on": ".",
            }
        )
        chi = g.create_dataset("chi", data=[chi])
        chi.attrs.update(
            {
                "transformation_type": "rotation",
                "vector": [0, 0, 1.0],
                "units": "deg",
                "depends_on": "omega",
            }
        )
        phi = g.create_dataset("phi", data=[0.0])
        phi.attrs.update(
            {
                "transformation_type": "rotation",
                "vector": [-1.0, -0.002, 0.001],
                "units": "deg",
                "depends_on": "chi",
            }
        )


def expected_panel():
    # p -> two_theta(det_z(module_offset(p))), McStas, then imgCIF; a
    # transformation with an offset is R p + offset, the offset after.
    r = rotation([-1, 0, 0], TWO_THETA)
    inner = (
        np.array([0, 0, DISTANCE_M * 1000])
        + OFFSET_VEC / np.linalg.norm(OFFSET_VEC) * OFFSET_M * 1000
    )
    origin = r @ inner + PIVOT_M * 1000
    return imgcif(origin), imgcif(r @ [-1.0, 0, 0]), imgcif(r @ [0, -1.0, 0])


def run(tmp_path, *extra):
    plant(tmp_path / "master.nxs")
    result = subprocess.run(
        [
            IMPORT,
            str(tmp_path / "master.nxs"),
            "-o",
            str(tmp_path / "imported.expt"),
            *extra,
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    return json.load(open(tmp_path / "imported.expt")), result.stdout


@needs_import
def test_a_planted_masters_geometry_is_read_as_worked_out_here(tmp_path):
    e, _ = run(tmp_path)
    p = e["detector"][0]["panels"][0]
    origin, fast, slow = expected_panel()
    assert np.allclose(p["origin"], origin, atol=1e-9)
    assert np.allclose(p["fast_axis"], fast, atol=1e-12)
    assert np.allclose(p["slow_axis"], slow, atol=1e-12)
    assert p["image_size"] == [200, 300] and np.allclose(
        p["pixel_size"], [0.075, 0.075]
    )
    assert p["trusted_range"] == [0.0, 20000.0] and p["thickness"] == pytest.approx(
        0.45
    )
    assert p["material"] == "Si" and p["mu"] > 0
    g = e["goniometer"][0]
    assert g["names"] == ["phi", "chi", "omega"] and g["scan_axis"] == 2
    assert np.allclose(g["axes"][0], imgcif([-1.0, -0.002, 0.001]))
    assert np.allclose(g["axes"][1], imgcif([0, 0, 1.0]))
    assert g["angles"] == [0.0, CHI, 0.0]
    s = e["scan"][0]
    assert s["image_range"] == [1, 10]
    assert np.allclose(s["properties"]["oscillation"], np.arange(10) * 0.5 + 10.0)
    assert np.allclose(s["properties"]["exposure_time"], 0.01)
    assert e["beam"][0]["wavelength"] == pytest.approx(1.0)


def beam_meets(panel):
    o, f, s = (
        np.asarray(panel[k], float) for k in ("origin", "fast_axis", "slow_axis")
    )
    n = np.cross(f, s)
    t = o.dot(n) / np.array([0, 0, -1.0]).dot(n)
    hit = np.array([0, 0, -t])
    return (
        (hit - o).dot(f) / panel["pixel_size"][0],
        (hit - o).dot(s) / panel["pixel_size"][1],
        abs(o.dot(n / np.linalg.norm(n))),
    )


@needs_import
def test_the_overrides_do_what_they_say(tmp_path):
    e, out = run(
        tmp_path,
        "--beam-centre",
        "123.5,77.25",
        "--distance",
        "180",
        "--wavelength",
        "0.8",
        "--image-range",
        "3,7",
    )
    p = e["detector"][0]["panels"][0]
    bx, by, d = beam_meets(p)
    assert (bx, by) == pytest.approx((123.5, 77.25), abs=1e-6)
    assert d == pytest.approx(180.0, abs=1e-6)
    _, fast, slow = expected_panel()
    assert np.allclose(p["fast_axis"], fast) and np.allclose(
        p["slow_axis"], slow
    )  # moved, not turned
    assert e["beam"][0]["wavelength"] == pytest.approx(0.8)
    e0, _ = run(tmp_path)
    assert (
        p["mu"] < e0["detector"][0]["panels"][0]["mu"]
    )  # shorter wavelength, less absorption
    s = e["scan"][0]
    assert s["image_range"] == [3, 7] and len(s["properties"]["oscillation"]) == 5
    assert "the beam centre given" in out and "the distance given" in out


@needs_import
def test_the_insulin_master_against_dials_import(tmp_path):
    master = os.environ.get("MXI_TEST_IMAGES")
    reference = os.environ.get("MXI_DIALS_IMPORTED")
    if not (
        master and reference and os.path.exists(master) and os.path.exists(reference)
    ):
        pytest.skip("set MXI_TEST_IMAGES and MXI_DIALS_IMPORTED")
    subprocess.run(
        [IMPORT, master, "-o", str(tmp_path / "ours.expt")],
        check=True,
        capture_output=True,
    )
    ours = json.load(open(tmp_path / "ours.expt"))
    theirs = json.load(open(reference))
    assert ours["goniometer"] == theirs["goniometer"]
    p, q = ours["detector"][0]["panels"][0], theirs["detector"][0]["panels"][0]
    assert p["mu"] == pytest.approx(q["mu"], rel=1e-12)
    for k in (
        "fast_axis",
        "slow_axis",
        "pixel_size",
        "trusted_range",
        "thickness",
        "material",
    ):
        assert p[k] == pytest.approx(q[k]) if k != "material" else p[k] == q[k], k
    assert ours["scan"][0]["image_range"] == theirs["scan"][0]["image_range"]
    assert np.allclose(
        ours["scan"][0]["properties"]["oscillation"],
        theirs["scan"][0]["properties"]["oscillation"],
    )


@needs_import
def test_several_masters_make_one_experiment_each_with_models_of_their_own(tmp_path):
    # As dials.import writes several sweeps: each experiment its own beam,
    # detector, goniometer, scan and image set, numbered past the ones before.
    plant(tmp_path / "a.nxs", chi=20.0)
    plant(tmp_path / "b.nxs", chi=45.0)
    result = subprocess.run(
        [
            IMPORT,
            str(tmp_path / "a.nxs"),
            str(tmp_path / "b.nxs"),
            "-o",
            str(tmp_path / "two.expt"),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    e = json.load(open(tmp_path / "two.expt"))
    assert len(e["experiment"]) == 2
    for i, x in enumerate(e["experiment"]):
        assert all(
            x[k] == i for k in ("beam", "detector", "goniometer", "scan", "imageset")
        )
    assert all(
        len(e[k]) == 2 for k in ("beam", "detector", "goniometer", "scan", "imageset")
    )
    assert e["goniometer"][0]["angles"] == [0.0, 20.0, 0.0]
    assert e["goniometer"][1]["angles"] == [0.0, 45.0, 0.0]
    assert e["experiment"][0]["identifier"] != e["experiment"][1]["identifier"]
    assert "2 sweeps, one experiment each" in result.stdout


@needs_import
def test_an_offset_without_offset_units_is_in_its_transformations_units(tmp_path):
    # nxmx's reading: Diamond's Eiger masters give offsets in metres with no
    # offset_units, and taking them as millimetres put the detector a
    # thousandth of the way out.
    off = np.array([0.01, 0.02, 0.0])
    plant(tmp_path / "master.nxs", module_offset_offset=off)
    result = subprocess.run(
        [IMPORT, str(tmp_path / "master.nxs"), "-o", str(tmp_path / "imported.expt")],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    p = json.load(open(tmp_path / "imported.expt"))["detector"][0]["panels"][0]
    r = rotation([-1, 0, 0], TWO_THETA)
    inner = (
        np.array([0, 0, DISTANCE_M * 1000])
        + OFFSET_VEC / np.linalg.norm(OFFSET_VEC) * OFFSET_M * 1000
        + off * 1000
    )
    assert np.allclose(p["origin"], imgcif(r @ inner + PIVOT_M * 1000), atol=1e-9)
    assert "taken in its units, m" in result.stdout


@needs_import
def test_compare_expt_compares_experiment_by_experiment(tmp_path):
    # Each experiment against its counterpart, through its own model indices:
    # a change to the second sweep's goniometer alone is found, and found there.
    from mxeq import compare_expt

    plant(tmp_path / "a.nxs", chi=20.0)
    plant(tmp_path / "b.nxs", chi=45.0)
    subprocess.run(
        [
            IMPORT,
            str(tmp_path / "a.nxs"),
            str(tmp_path / "b.nxs"),
            "-o",
            str(tmp_path / "two.expt"),
        ],
        check=True,
        capture_output=True,
    )
    same = compare_expt.compare_files(
        str(tmp_path / "two.expt"), str(tmp_path / "two.expt")
    )
    assert same.differences == 0
    e = json.load(open(tmp_path / "two.expt"))
    e["goniometer"][1]["angles"] = [0.0, 50.0, 0.0]
    json.dump(e, open(tmp_path / "changed.expt", "w"))
    changed = compare_expt.compare_files(
        str(tmp_path / "two.expt"), str(tmp_path / "changed.expt")
    )
    assert changed.differences == 1
    text = "\n".join(changed.lines)
    second = text.index("experiment 1 against")
    differing = [
        k
        for k, line in enumerate(text.splitlines())
        if line.strip().startswith("angles") and not line.endswith("same")
    ]
    starts = [len(line) + 1 for line in text.splitlines()]
    offsets = np.cumsum([0] + starts)[:-1]
    assert len(differing) == 1 and offsets[differing[0]] > second


def trusted(tmp_path, **planted):
    tmp_path.mkdir(parents=True, exist_ok=True)
    plant(tmp_path / "master.nxs", **planted)
    result = subprocess.run(
        [IMPORT, str(tmp_path / "master.nxs"), "-o", str(tmp_path / "imported.expt")],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    p = json.load(open(tmp_path / "imported.expt"))["detector"][0]["panels"][0]
    return p["trusted_range"][1], result.stdout


@needs_import
def test_the_trusted_range_is_the_lower_of_the_detectors_limit_and_the_types(tmp_path):
    # The detector's count limit, set by its exposure time, and the data type's:
    # 2^bits - 3, the two largest values being markers.
    top, _ = trusted(tmp_path / "a", meta={"cutoff": 133201, "bits": 32})
    assert top == 133201.0
    top, out = trusted(tmp_path / "b", meta={"cutoff": 100000, "bits": 16})
    assert top == 65533.0 and "is above what 16-bit data can hold" in out
    top, out = trusted(tmp_path / "c", meta={"bits": 16})
    assert top == 65533.0 and "no count limit could be read" in out


@needs_import
def test_a_link_into_a_missing_meta_file_is_said_so(tmp_path):
    # Not taken for absent: where it leads, and a fallback that distrusts no
    # count for its size, as dxtbx takes a file without one.
    top, out = trusted(tmp_path, meta={"cutoff": 133201, "bits": 32, "write": False})
    assert top == 2147483647.0
    assert (
        "saturation_value links to /_dectris/countrate_correction_count_cutoff in master_meta.h5"
        in out
    )
    assert "which cannot be opened" in out


# mu in 1/mm, from cctbx's eltbx tables -- parsed from its
# attenuation_coefficient.cpp and looked up by its own rule (the first point
# above the energy, the interval before it, log-log), independently of
# mxi_import's code -- as dials.import computes it. Either side of the Cd K
# edge at 26.711 keV, and past the table's end.
MU_FROM_CCTBX = [
    ("Silicon", 1.0, 4.207985875566818),
    ("CdTe", 1.0, 48.11463840098675),
    ("GaAs", 1.0, 80.07998629661657),
    ("CdTe", 0.46463323212022556, 6.113517802096811),  # 0.1 per cent below the edge
    ("CdTe", 0.46370489399411124, 18.199624769356834),  # 0.1 per cent above it
    ("CdTe", 0.3, 10.97168695998142),
]


@needs_import
@pytest.mark.parametrize("material, wavelength, mu", MU_FROM_CCTBX)
def test_mu_is_cctbx_s_for_each_sensor_material(tmp_path, material, wavelength, mu):
    plant(tmp_path / "master.nxs", material=material.encode(), wavelength=wavelength)
    result = subprocess.run(
        [IMPORT, str(tmp_path / "master.nxs"), "-o", str(tmp_path / "imported.expt")],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    p = json.load(open(tmp_path / "imported.expt"))["detector"][0]["panels"][0]
    assert p["mu"] == pytest.approx(mu, rel=1e-9)
    assert "no attenuation coefficient tabulated" not in result.stdout
