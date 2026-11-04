"""Summary records of a whole chain, made in the shapes mxi's commands write
them (apps/mxi_*.cc's summary.set and summary.table: a table column by
column, its names under "columns") -- invented numbers, insulin-like, for
testing and drawing mxi ui's summaries and plots, never as results."""

import math


def table(**columns):
    names = list(columns)
    out = {"columns": names}
    out.update({k: list(v) for k, v in columns.items()})
    return out


def shells(n=10, d_min=1.6, d_max=40.0):
    """n shells equal in 1/d^2, from d_max to d_min: (d_max, d_min) pairs."""
    lo, hi = 1 / d_max**2, 1 / d_min**2
    edges = [lo + (hi - lo) * i / n for i in range(n + 1)]
    return [(1 / math.sqrt(a), 1 / math.sqrt(b)) for a, b in zip(edges, edges[1:])]


def chain():
    s = shells()
    fall = [math.exp(-2.2 * i / 9) for i in range(10)]
    return {
        "import": {
            "program": "mxi_import",
            "time": "2026-11-04 10:00:00",
            "sweeps": table(sweep=[0], images=[1800], rotation=[180.0]),
        },
        "find": {
            "program": "mxi_find",
            "time": "2026-11-04 10:01:00",
            "spots": 120000,
            "signal_pixels": 900000,
            "first_image": [1],
            "spots_per_image": [[60 + (i % 7) for i in range(1800)]],
            "spots_by_resolution": table(
                d_max=[a for a, _ in s],
                d_min=[b for _, b in s],
                spots=[3000 * f for f in fall],
            ),
        },
        "index": {
            "program": "mxi_index",
            "time": "2026-11-04 10:02:00",
            "cell": [78.1, 78.1, 78.1, 90.0, 90.0, 90.0],
            "indexed_fraction": 0.982,
            "lattices": 1,
            "reflections": 120000,
            "volume": 476379.0,
            "cycles": table(
                cycle=[1, 2, 3],
                refined_on=[50000, 52000, 52100],
                rejected=[400, 120, 90],
                indexed=[117000, 117800, 117900],
                rmsd_x_px=[0.42, 0.31, 0.30],
                rmsd_y_px=[0.40, 0.30, 0.29],
                rmsd_z_images=[0.35, 0.22, 0.21],
            ),
            "a_figure_added_later": 7,
        },
        "refine": {
            "program": "mxi_refine",
            "time": "2026-11-04 10:03:00",
            "reflections": 117900,
            "refined_on": 117000,
            "rejected": 900,
            "rmsd_px_px_images": [0.21, 0.20, 0.15],
            "scan_varying": 1,
            "cycles": table(
                **{
                    "pass": [1, 1, 2, 2],
                    "cycle": [1, 2, 1, 2],
                    "reflections": [117900] * 4,
                    "rejected": [900] * 4,
                    "rmsd_x_px": [0.3, 0.25, 0.22, 0.21],
                    "rmsd_y_px": [0.29, 0.24, 0.21, 0.20],
                    "rmsd_z_images": [0.2, 0.17, 0.16, 0.15],
                }
            ),
            "cell_over_scan": [
                {
                    "experiment": 0,
                    "angle": list(range(0, 181, 10)),
                    **{
                        k: [78.1 + 0.01 * math.sin(i / 3) for i in range(19)]
                        for k in ("a", "b", "c")
                    },
                    **{k: [90.0] * 19 for k in ("alpha", "beta", "gamma")},
                }
            ],
        },
        "integrate": {
            "program": "mxi_integrate",
            "time": "2026-11-04 10:05:00",
            "by_resolution": table(
                d_max=[a for a, _ in s],
                d_min=[b for _, b in s],
                isigma_sum=[40 * f for f in fall],
                isigma_prf=[44 * f for f in fall],
                cc_prf=[0.99 - 0.02 * i for i in range(10)],
                ibg=[5.0 - 0.3 * i for i in range(10)],
            ),
            "overall": table(
                reflections=[310000], isigma_sum=[12.1], isigma_prf=[13.4]
            ),
        },
        "symmetry": {
            "program": "mxi_symmetry",
            "time": "2026-11-04 10:06:00",
            "laue_group": "I m -3",
            "space_group": "I 2 3",
            "reindex": "a,b,c",
            "subgroups": table(
                likelihood=[0.97, 0.02],
                laue_group=["I m -3", "I 4/m m m"],
                z_cc_net=[9.1, 4.2],
            ),
        },
        "scale": {
            "program": "mxi_scale",
            "time": "2026-11-04 10:08:00",
            "observations": 310000,
            "unique": 21000,
            "outliers": 140,
            "resolution_cc_half": 1.62,
            "error_model_a": 1.02,
            "error_model_b": 0.031,
            "merging_by_resolution": table(
                shell=list(range(10)),
                d_max=[a for a, _ in s],
                d_min=[b for _, b in s],
                completeness=[99.8 - 0.3 * i for i in range(10)],
                multiplicity=[15 - 0.4 * i for i in range(10)],
                i_over_sigma=[60 * f for f in fall],
                cc_half=[0.999 - 0.06 * i for i in range(10)],
                rmerge=[0.03 / f for f in fall],
                rmeas=[0.031 / f for f in fall],
                rpim=[0.008 / f for f in fall],
            ),
            "merging_summary": table(
                shell=["overall", "low", "high"],
                d_max=[40.0, 40.0, 1.69],
                d_min=[1.6, 4.3, 1.6],
                completeness=[99.1, 99.9, 97.2],
                multiplicity=[13.2, 14.8, 11.0],
                i_over_sigma=[21.5, 61.0, 3.1],
                rmerge=[0.052, 0.031, 0.61],
                rmeas=[0.054, 0.032, 0.64],
                rpim=[0.015, 0.008, 0.19],
                cc_half=[0.999, 0.999, 0.74],
                observations=[310000, 16000, 25000],
                unique=[21000, 1100, 2300],
            ),
        },
    }


def install(monkeypatch, directory, records=None):
    """Each step's output file in the directory, its summary the steps up to
    it -- read through mxi_tools.read_summary, replaced."""
    # The mxi_tools results.py holds: other tests load mxi_tools afresh (the
    # build's copy), and patching whichever "import mxi_tools" finds then
    # missed the one in use.
    from mxi_ui import results
    from mxi_ui.steps import STEPS

    mxi_tools = results.mxi_tools

    records = records or chain()
    order = [s.key for s in STEPS if s.key in records]
    by_file = {}
    for i, key in enumerate(order):
        step = next(s for s in STEPS if s.key == key)
        (directory / step.output).write_text("")
        by_file[str(directory / step.output)] = [records[k] for k in order[: i + 1]]
    monkeypatch.setattr(
        mxi_tools, "read_summary", lambda path: {"steps": by_file.get(str(path), [])}
    )
    return records
