"""mxeq observations lists every observation of a reflection, its
equivalents found by the experiment list's space group, and lists a table
without mxi's columns -- DIALS's -- as well as one with them."""

import json

import numpy as np

from mxeq import observations, refl


def write(tmp_path, with_mxi_columns):
    e = {
        "__id__": "ExperimentList",
        "experiment": [{"__id__": "Experiment", "crystal": 0}],
        "crystal": [
            {
                "__id__": "crystal",
                "real_space_a": [50.0, 0, 0],
                "real_space_b": [0, 50.0, 0],
                "real_space_c": [0, 0, 70.0],
                "space_group_hall_symbol": " P 4",
            }
        ],
    }
    (tmp_path / "a.expt").write_text(json.dumps(e))
    hkl = np.array(
        [[1, 2, 3], [-2, 1, 3], [-1, -2, 3], [2, -1, 3], [1, 2, 4]], dtype=np.int32
    )
    n = len(hkl)
    t = refl.ReflectionTable(nrows=n)
    cols = {
        "miller_index": (hkl, "cctbx::miller::index<>"),
        "flags": (np.full(n, (1 << 8) | (1 << 9), dtype=np.int64), "std::size_t"),
        "xyzcal.px": (np.arange(3 * n, dtype=float).reshape(n, 3), "vec3<double>"),
        "intensity.sum.value": (np.arange(n, dtype=float) * 10, "double"),
        "intensity.prf.value": (np.arange(n, dtype=float) * 11, "double"),
    }
    if with_mxi_columns:
        cols["profile.failure"] = (np.array([0, 0, 4, 0, 0], dtype=np.int32), "int")
    for name, (v, kind) in cols.items():
        t.columns[name] = v
        t.types[name] = kind
    refl.write(str(tmp_path / "a.refl"), t)


def test_every_equivalent_is_listed_by_the_space_group(tmp_path):
    write(tmp_path, True)
    text = observations.listing(
        str(tmp_path / "a.expt"), str(tmp_path / "a.refl"), [[1, 2, 3]]
    )
    assert "== [1, 2, 3]: 4 observations" in text  # the four of P4, not [1, 2, 4]
    assert " why" in text and "[1, 2, 4]" not in text.split("observations")[1]


def test_a_table_without_mxis_columns_lists_too(tmp_path):
    write(tmp_path, False)
    text = observations.listing(
        str(tmp_path / "a.expt"), str(tmp_path / "a.refl"), [[2, -1, 3]]
    )
    assert "4 observations" in text and " why" not in text
