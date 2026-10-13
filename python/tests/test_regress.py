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
