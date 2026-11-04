"""mxi ui's results from a step's summary record: the record its file holds,
its facts, its tables, its plots -- without a window."""

import math

from mxi_ui import results
from ui_records import chain, install


def test_a_steps_record_is_its_own_programs_from_its_own_file(tmp_path, monkeypatch):
    install(monkeypatch, tmp_path)
    path, rec = results.record(str(tmp_path), "refine")
    assert path.endswith("refined.rflx") and rec["program"] == "mxi_refine"
    path, rec = results.record(str(tmp_path), "export")  # an .mtz: no summary
    assert rec is None
    (tmp_path / "scaled.rflx").unlink()
    assert results.record(str(tmp_path), "scale")[1] is None


def test_facts_known_first_then_any_other_figure():
    facts = dict(results.facts(chain()["index"]))
    assert facts["Indexed fraction"] == "98.2%"
    assert facts["Cell"] == "78.1  78.1  78.1  90  90  90"
    assert (
        facts["A figure added later"] == "7"
    ), "a figure this does not know is still shown"
    assert "Cycles" not in facts, "a table is not a fact"


def test_scalings_merging_summary_laid_out_as_dials_does():
    (merging,) = [
        t for t in results.tables(chain()["scale"]) if t.title == "Merging statistics"
    ]
    assert merging.columns == ["", "overall", "low", "high"]
    rows = {r[0]: r[1:] for r in merging.rows}
    assert rows["High resolution limit"] == ["1.60", "4.30", "1.60"]
    assert rows["CC half"] == ["0.999", "0.999", "0.740"]
    assert list(rows)[:3] == [
        "High resolution limit",
        "Low resolution limit",
        "Completeness",
    ]


def test_each_step_has_its_plots():
    c = chain()
    titles = {k: [p.title for p in results.plots(c[k])] for k in c}
    assert titles["find"] == ["Spots per image", "Spots by resolution"]
    assert titles["index"] == ["RMSDs, cycle by cycle"]
    assert titles["refine"][0] == "RMSDs, cycle by cycle" and len(titles["refine"]) == 3
    assert titles["integrate"] == [
        "<I/σ> by resolution",
        "Profile correlation by resolution",
        "Background by resolution",
    ]
    assert titles["scale"] == [
        "CC½ by resolution",
        "<I/σ> by resolution",
        "Completeness and multiplicity",
        "R factors by resolution",
    ]
    assert titles["symmetry"] == [] and titles["import"] == []
    scale = results.plots(c["scale"])
    assert scale[3].log and all(p.resolution for p in scale)
    assert scale[2].right and scale[2].right[0].name == "multiplicity"


def test_a_resolution_axis_is_in_one_over_d_squared_ticked_in_d():
    plot = results.plots(chain()["integrate"])[0]
    xs = plot.series[0].x
    assert all(b > a for a, b in zip(xs, xs[1:])), "1/d^2 rising with resolution"
    positions, labels = results.d_ticks(xs)
    ds = [float(x) for x in labels]
    assert all(b < a for a, b in zip(ds, ds[1:])) and math.isclose(
        1 / ds[-1] ** 2, positions[-1], rel_tol=1e-2
    )


def test_a_record_with_missing_or_odd_values_does_not_break_it():
    odd = {
        "program": "mxi_scale",
        "merging_by_resolution": {
            "columns": ["d_max", "d_min", "cc_half"],
            "d_max": [None, 3.0],
            "d_min": [3.0, 2.0],
            "cc_half": [float("nan"), "x"],
        },
    }
    plots = results.plots(odd)
    assert plots and plots[0].series[0].y == [None, None]
    assert (
        results.facts(None) == []
        and results.plots(None) == []
        and results.tables({}) == []
    )


def test_cycles_and_images_are_counted_in_whole_numbers():
    c = chain()
    assert results.plots(c["refine"])[0].integer and results.plots(c["find"])[0].integer
    assert not results.plots(c["scale"])[0].integer
