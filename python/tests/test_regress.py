"""mxeq regress (docs/regression.md): the figures a run is judged by, read from
a processing summary, and each figure's sense of better and worse. The runner
itself runs real data sets, and is tried by hand."""

from mxeq import regress


def test_the_figures_come_from_each_programs_last_step():
    summary = {
        "steps": [
            {"program": "mxi_index", "indexed": 900, "reflections": 1000},
            {"program": "mxi_index", "indexed": 980, "reflections": 1000},  # again
            {"program": "mxi_refine", "rmsd_px_px_images": [0.3, 0.25, 0.2]},
            {
                "program": "mxi_scale",
                "merging_summary": {
                    "columns": ["shell", "rmerge", "cc_half"],
                    "shell": ["overall", "low", "high"],
                    "rmerge": [0.05, 0.03, 0.4],
                    "cc_half": [0.99, 1.0, 0.5],
                },
                "two_theta_cell": [78.0, 78.0, 78.0, 90.0, 90.0, 90.0],
                "error_model_a": 1.1,
            },
        ]
    }
    f = regress.figures(summary)
    assert f["index.indexed_fraction"] == 0.98  # the later step
    assert f["refine.rmsd_y_px"] == 0.25
    assert f["scale.overall.rmerge"] == 0.05
    assert f["scale.high.cc_half"] == 0.5
    assert f["cell.a"] == 78.0
    assert f["scale.error_model_a"] == 1.1


def test_each_figure_knows_which_way_is_better():
    assert regress.sense_of("scale.overall.rmerge") == -1
    assert regress.sense_of("scale.overall.cc_half") == +1
    assert regress.sense_of("refine.rmsd_x_px") == -1
    assert regress.sense_of("time.integrate") == -1
    assert regress.sense_of("cell.a") == 0  # a change is shown, not judged
    assert regress.tolerance_of("time.find") == (0.10, 0.5)


def test_the_tiers_are_in_order_and_a_named_data_set_runs_whatever_its_tier():
    sets = {
        "a": {"tier": "quick"},
        "b": {},  # full
        "c": {"tier": "slow"},
    }
    assert set(regress.choose(sets, upto="quick")) == {"a"}
    assert set(regress.choose(sets)) == {"a", "b"}  # the default
    assert set(regress.choose(sets, upto="slow")) == {"a", "b", "c"}
    assert set(regress.choose(sets, tier="slow")) == {"c"}
    assert set(regress.choose(sets, names=["c"])) == {"c"}
    import pytest

    with pytest.raises(SystemExit):
        regress.choose({"x": {"tier": "fast"}})  # not a tier
    with pytest.raises(SystemExit):
        regress.choose(sets, names=["nope"])


def test_clean_keeps_the_recent_and_the_tagged_and_every_runs_figures(
    tmp_path, monkeypatch, capsys
):
    import types

    monkeypatch.setenv("MXI_REGRESSION", str(tmp_path))
    db = regress.database()
    for k, tag in enumerate([None, "reference", None, None]):
        run = f"2026101{k}-000000"
        db.execute(
            "INSERT INTO runs VALUES (?, ?, NULL, ?, 0, 'h', 'p', 'b', ?)",
            (run, f"2026-10-1{k}T00:00:00", f"c{k}" * 6, tag),
        )
        db.execute("INSERT INTO metrics VALUES (?, 'x', 'time.total', 1.0)", (run,))
        (tmp_path / "runs" / run / "x").mkdir(parents=True)
        (tmp_path / "runs" / run / "x" / "big").write_bytes(b"0" * 1000)
    db.commit()
    args = types.SimpleNamespace(keep_runs=1, keep_builds=0, yes=False)
    regress.cmd_clean(args)  # says what, removes nothing
    assert len(list((tmp_path / "runs").iterdir())) == 4
    args.yes = True
    regress.cmd_clean(args)
    left = sorted(d.name for d in (tmp_path / "runs").iterdir())
    assert left == ["20261011-000000", "20261013-000000"]  # tagged, newest
    figures = db.execute("SELECT COUNT(*) FROM metrics").fetchone()[0]
    assert figures == 4  # every run's figures stay
