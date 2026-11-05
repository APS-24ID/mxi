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
    assert facts["Cell"] == "78.100  78.100  78.100  90.000  90.000  90.000"
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
    assert not scale[3].log, "R factors on a linear axis, as usual"
    assert all(p.resolution for p in scale)
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


def test_large_numbers_whole_and_small_to_four_figures():
    facts = dict(
        results.facts({"program": "mxi_index", "volume": 236118.4, "d_min": 1.99395})
    )
    assert facts["Volume (Å³)"] == "236118" and facts["Resolution limit (Å)"] == "1.994"


def test_every_h5py_is_opened_with_the_filters():
    # The summary's longer arrays, like the reflections', are compressed:
    # opened with h5py alone, Graeme's indexed.rflx read no summary (4
    # November), as its reflections had not (23 October). Every import of
    # h5py in mxi_tools and mxi_ui is h5py_with_filters', none on its own.
    import ast
    import pathlib

    src = pathlib.Path(results.__file__).resolve().parents[1]
    files = [src / "mxi_tools.py", *sorted((src / "mxi_ui").glob("*.py"))]
    stray = []
    for path in files:
        tree = ast.parse(path.read_text())
        for function in [
            n for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.Module))
        ]:
            body = (
                function.body
                if isinstance(function, ast.FunctionDef)
                else function.body
            )
            for node in (
                body if isinstance(function, ast.Module) else ast.walk(function)
            ):
                if isinstance(node, ast.Import) and any(
                    a.name == "h5py" for a in node.names
                ):
                    owner = (
                        function.name
                        if isinstance(function, ast.FunctionDef)
                        else "<module>"
                    )
                    if owner != "h5py_with_filters":
                        stray.append(f"{path.name}:{node.lineno} in {owner}")
    assert not stray, stray


def test_a_table_cell_that_is_a_list_is_compact():
    rec = {
        "program": "mxi_index",
        "candidates": {
            "columns": ["cell", "volume", "indexed", "chosen"],
            "cell": [
                [
                    67.28432387505138,
                    67.44543540551338,
                    67.49954671142541,
                    109.51,
                    109.65,
                    109.22,
                ]
            ],
            "volume": [235918.4],
            "indexed": [154385],
            "chosen": [1],
        },
    }
    (table,) = results.tables(rec)
    assert table.rows == [
        ["67.28 67.45 67.50 109.51 109.65 109.22", "235918", "154385", "1"]
    ]


def test_a_clusters_records_read_as_facts():
    find = {
        "program": "mxi_find",
        "spots": 40,
        "whole_run_left_out": ["spots_by_resolution"],
    }
    assert (
        dict(results.facts(find))["Left out (the whole run's)"] == "spots_by_resolution"
    )
    own = {
        "program": "mxi_correlation",
        "cluster": 1,
        "data_sets": 12,
        "of_data_sets": 36,
        "cut_from": "scaled.rflx",
    }
    facts = dict(results.facts(own))
    assert (
        facts["Cluster"] == "1" and facts["Data sets"] == "12" and facts["Of"] == "36"
    )
    assert facts["Cut from"] == "scaled.rflx"


def test_correlations_facts_and_plots(tmp_path):
    import json

    from ui_records import correlation

    doc = correlation()
    (tmp_path / "mxi_correlation.json").write_text(
        json.dumps(doc).replace("Infinity", "Infinity")
    )
    path, rec = results.record(str(tmp_path), "correlation")
    assert (
        path.endswith("mxi_correlation.json")
        and rec["program"] == "mxi_correlation_json"
    )
    facts = dict(results.facts(rec))
    assert (
        facts["Data sets"] == "6"
        and facts["Clusters"] == "2"
        and facts["In no cluster"] == "1"
    )
    assert facts["Their data sets"] == "0: 3, 1: 2"
    plots = results.plots(rec)
    assert [p.title for p in plots] == [
        "Correlation, in OPTICS's order",
        "Cosine angles, in OPTICS's order",
        "Reachability, in OPTICS's order",
        "The data sets, principal components",
        "Functional by dimension",
    ]
    matrix = plots[0]
    order = doc["ordering"]
    assert (
        matrix.kind == "matrix" and matrix.matrix[2][3] == doc["r"][order[2]][order[3]]
    )
    reach = plots[2]
    assert reach.bars and [s.name for s in reach.series] == [
        "in no cluster",
        "cluster 0",
        "cluster 1",
    ]
    first = next(s for s in reach.series if s.name == "cluster 0")
    assert (
        first.x[0] == 0 and first.y[0] > 1.5
    ), "the unreachable first drawn, above the rest"
    pcs = plots[3]
    assert (
        pcs.kind == "points" and pcs.xlabel == "PC1 (80%)" and pcs.ylabel == "PC2 (15%)"
    )
    assert results.record(str(tmp_path / "nowhere"), "correlation")[1] is None


def test_integration_a_sweep_at_a_time():
    from ui_records import chain

    rec = chain()["integrate"]
    sweep = {"by_resolution": rec["by_resolution"]}
    assert (
        results.each_sweep_plots(dict(rec, each_sweep=[dict(sweep, sweep=0)])) == []
    ), "one sweep, the whole run's tables: nothing more"
    two = dict(rec, each_sweep=[dict(sweep, sweep=0), dict(sweep, sweep=1)])
    plots = results.each_sweep_plots(two)
    assert [p.title for p in plots] == [
        "<I/σ> by resolution, each sweep",
        "Profile correlation by resolution, each sweep",
    ]
    assert [s.name for s in plots[0].series] == ["sweep 0", "sweep 1"]
    cut = {"program": "mxi_integrate", "each_sweep": [dict(sweep, sweep=0)]}
    assert results.each_sweep_plots(cut), "a cluster's: only the sweeps left"
    assert results.plots(two)[0].title == "<I/σ> by resolution, each sweep"


def test_the_cosine_angles_as_the_report_computes_them():
    import math

    from ui_records import correlation

    doc = correlation()
    n = len(doc["labels"])
    cos = results.cosine_angles(doc["coordinates"], n)
    x = [[doc["coordinates"][k * n + i] for k in range(2)] for i in range(n)]

    def angle(i, j):
        dot = x[i][0] * x[j][0] + x[i][1] * x[j][1]
        return dot / (math.hypot(*x[i]) * math.hypot(*x[j]))

    assert all(math.isclose(cos[i][j], angle(i, j)) for i in range(n) for j in range(n))
    assert cos[0][1] > 0.9 and cos[0][2] < -0.9, "a cluster together, the other opposed"
    assert (
        results.cosine_angles([], 6) == [] and results.cosine_angles([1.0] * 7, 6) == []
    )
    assert results.cosine_angles([0.0] * 12, 6)[0][1] == 0.0, "no length, no angle"
    plot = [p for p in results.correlation_plots(doc) if p.title.startswith("Cosine")][
        0
    ]
    order = doc["ordering"]
    assert (
        plot.matrix[2][3] == cos[order[2]][order[3]] and plot.low < 0
    ), "ordered; its scale from -1"


def test_a_matrix_with_no_value_draws_none_there():
    from ui_records import correlation

    doc = correlation()
    doc["r"][0][1] = None
    plot = results.correlation_plots(doc)[0]
    i, j = doc["ordering"].index(0), doc["ordering"].index(1)
    assert plot.matrix[i][j] != plot.matrix[i][j], "NaN, not drawn"
    assert plot.matrix[i][i] == 1.0, "the diagonal 1"


def test_many_sweeps_grey_under_their_median_and_a_bar_each():
    # 36 crystals' sweeps, a line each with a legend of 36, had buried the
    # plot: each grey and fine, the median bold, and a bar a sweep.
    from ui_records import chain

    rec = chain()["integrate"]
    each = []
    for k in range(36):
        t = dict(rec["by_resolution"])
        t["isigma_prf"] = [v * (1 + k / 100) for v in t["isigma_prf"]]
        overall = {
            "columns": ["isigma_prf", "cc_prf"],
            "isigma_prf": [10 + k],
            "cc_prf": [0.9],
        }
        each.append({"sweep": k, "by_resolution": t, "overall": overall})
    plots = results.each_sweep_plots(dict(rec, each_sweep=each))
    assert [p.title for p in plots] == [
        "<I/σ> by resolution, 36 sweeps",
        "<I/σ>, each sweep",
        "Profile correlation, each sweep",
    ]
    lines = plots[0].series
    assert (
        sum(s.thin for s in lines) == 36
        and lines[-1].name == "median"
        and not lines[-1].thin
    )
    first = sorted(s.y[0] for s in lines[:36])
    assert math.isclose(lines[-1].y[0], (first[17] + first[18]) / 2), "the median of 36"
    bars = plots[1]
    assert (
        bars.bars and bars.series[0].x == list(range(36)) and bars.series[0].y[35] == 45
    )
    few = results.each_sweep_plots(dict(rec, each_sweep=each[:8]))
    assert (
        few[0].title == "<I/σ> by resolution, each sweep" and len(few[0].series) == 8
    ), "eight, a line each, as before"
