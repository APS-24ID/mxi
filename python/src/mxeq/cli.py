"""``mxeq`` -- equivalence checks between two runs of an MX pipeline.

    mxeq check strong      a.refl b.refl
    mxeq check indexed     a.refl b.refl -e a.expt
    mxeq check refined     a.expt b.expt
    mxeq check integrated  a.refl b.refl -e a.expt [--operator ...]
    mxeq check scaled      a.refl b.refl -e a.expt [--operator ...]
    mxeq inspect           file.refl

Exit status is zero whenever the check ran.  It is not a verdict: this version
applies no thresholds, and a non-zero status means the comparison could not be
made, not that the two disagreed.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

from . import expt, refl
from .checks import indexed, integrated, refined, scaled, strong

BOUNDARIES = ("strong", "indexed", "refined", "integrated", "scaled")


def _operator(text: str | None) -> np.ndarray | None:
    if not text:
        return None
    parts = [p for p in text.replace(",", " ").split() if p]
    if len(parts) != 9:
        raise SystemExit(f"--operator needs nine integers, got {len(parts)}")
    return np.array([int(p) for p in parts], dtype=np.int64).reshape(3, 3)


def _load_expt(path: str | None) -> expt.ExperimentList | None:
    return None if path is None else expt.load(path)


def _guess(table: refl.ReflectionTable) -> str:
    if "intensity.scale.value" in table or "inverse_scale_factor" in table:
        return "scaled"
    if "intensity.prf.value" in table or "intensity.sum.variance" in table:
        return "integrated"
    if "miller_index" in table:
        return "indexed"
    return "strong"


def _sniff(path):
    """Is this a reflection table or an experiment list? None if neither.

    By the first byte, which is enough and does not require reading either
    format. A msgpack reflection table begins 0x93, the header for a
    three-element array; an experiment list is JSON and begins with whitespace
    or a brace.
    """
    try:
        with open(path, "rb") as handle:
            first = handle.read(1)
    except OSError:
        return None
    if not first:
        return None
    if first[0] == 0x93:
        return "refl"
    if first[0] in b"{ \t\r\n":
        return "expt"
    return None


def build_parser() -> argparse.ArgumentParser:
    """The command line, without running it.

    Separate from main so that what the commands ARE can be asked of the
    program itself -- the documentation test checks every `mxeq` subcommand a
    document names against this -- rather than read out of the source.
    """
    parser = argparse.ArgumentParser(
        prog="mxeq", description="Compare two runs of an MX pipeline, stage by stage."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    c = sub.add_parser("check", help="compare two files at a pipeline boundary")
    c.add_argument("boundary", choices=(*BOUNDARIES, "auto"))
    c.add_argument("a")
    c.add_argument("b")
    c.add_argument("-e", "--expt", help="experiment list, for resolution and symmetry")
    c.add_argument("--expt-b", help="B's experiment list, if it differs from A's")
    c.add_argument(
        "--operator",
        help="nine integers: the reindexing operator from `mxeq check indexed`",
    )
    c.add_argument("--radius", type=float, default=2.0, help="spatial match radius, px")
    c.add_argument("--bins", type=int, default=10, help="number of resolution shells")
    c.add_argument("--json", action="store_true", help="emit JSON instead of text")

    t = sub.add_parser(
        "trend",
        help="where two integrated files disagree, binned against what might explain it",
    )
    t.add_argument("a", help="ours")
    t.add_argument("b", help="theirs, the reference")
    t.add_argument(
        "--value",
        action="append",
        help="a column to compare; repeatable. Defaults to the intensities, "
        "their variances and the background.",
    )
    t.add_argument(
        "--value-b",
        help="the column to compare against, if different. Pass the same file "
        "twice to compare our own summed and fitted intensities.",
    )
    t.add_argument("--bins", type=int, default=10, help="bins per variable")
    t.add_argument(
        "--worst",
        type=int,
        default=0,
        help="also list the N reflections that disagree most, with their "
        "indices and positions",
    )
    t.add_argument(
        "--radius",
        type=float,
        default=5.0,
        help="how many images apart two rows may be and still be the same "
        "observation (5). Generous on purpose: two observations of one "
        "reflection are a whole turn apart, so this guards against nothing and "
        "a tight value loses real pairs",
    )

    hr = sub.add_parser(
        "html",
        help="an HTML report of the same comparison as trend, with graphs",
    )
    hr.add_argument("a", help="ours")
    hr.add_argument("b", help="theirs, the reference")
    hr.add_argument("-o", "--output", default="comparison.html", help="where to write")
    hr.add_argument(
        "--value",
        action="append",
        help="a column to compare; repeatable. Defaults to the intensities, "
        "their variances and the background.",
    )
    hr.add_argument("--bins", type=int, default=12, help="bins per variable")
    hr.add_argument("--title", default="integration comparison")
    hr.add_argument(
        "--radius",
        type=float,
        default=5.0,
        help="how many images apart two rows may be and still be the same "
        "observation (5)",
    )
    hr.add_argument(
        "--plotly",
        default=None,
        help="where the page should load plotly from; by default its CDN, "
        "which needs a network connection once to draw",
    )

    rs = sub.add_parser(
        "residuals",
        help="an HTML report of how well positions were predicted, from "
        "mxi_integrate's xyzres.px columns",
    )
    rs.add_argument("path", help="an integrated.refl written by mxi_integrate")
    rs.add_argument("-o", "--output", default="residuals.html", help="where to write")
    rs.add_argument(
        "--least-signal",
        type=float,
        default=10.0,
        help="only reflections of at least this I/sigma (10). For strong ones "
        "counting noise is a hundredth of a pixel and what is left is the "
        "prediction; lower it to see the weak ones too",
    )
    rs.add_argument("--bins", type=int, default=30, help="bins along each axis")
    rs.add_argument(
        "--cells",
        type=int,
        default=12,
        help="the detector map is this many cells a side",
    )
    rs.add_argument("--title", default="position residuals")

    ex = sub.add_parser(
        "explain",
        help="say why observations of two integrations went unpartnered",
    )
    ex.add_argument("a")
    ex.add_argument("b")
    ex.add_argument("--radius", type=float, default=5.0)
    ex.add_argument("--sample", type=int, default=20000)

    pl = sub.add_parser(
        "profiles", help="draw the reference profiles mxi_integrate learned"
    )
    pl.add_argument("path", help="what --save-profiles wrote")
    pl.add_argument(
        "-o", "--output", default="profiles.png", help="where to write the sections"
    )
    pl.add_argument(
        "--widths", help="also write a plot of profile width by region and axis"
    )
    pl.add_argument(
        "--block",
        type=int,
        help="draw only this block of the scan; with a divided scan there can "
        "be more profiles than fit on a page",
    )

    dg = sub.add_parser(
        "disagree",
        help="write a reflection table of only the reflections two "
        "integrations disagree about, to look at in the image viewer",
    )
    dg.add_argument("ours", help="our integrated.refl")
    dg.add_argument("theirs", help="the reference integrated.refl")
    dg.add_argument("-o", "--output", default="disagree.refl", help="where to write")
    dg.add_argument(
        "--value", default="intensity.sum.value", help="the column to compare"
    )
    dg.add_argument(
        "--value-b",
        help="the column to compare against, if different; pass the same file "
        "twice to compare our own summed and fitted intensities",
    )
    dg.add_argument(
        "--factor",
        type=float,
        help="RELATIVE: a/b outside [1/F, F], either way round. The default "
        "when no criterion is given, at 2",
    )
    dg.add_argument(
        "--difference",
        type=float,
        help="ABSOLUTE: |a - b| in counts. A background biased by a fraction "
        "of a count costs every reflection the same number of counts, and only "
        "this sees that as one thing",
    )
    dg.add_argument(
        "--sigma",
        type=float,
        help="the difference in units of the two variances added: the only "
        "criterion that knows whether a disagreement is larger than the "
        "measurement. Needs the matching variance column",
    )
    dg.add_argument(
        "--floor",
        type=float,
        default=5.0,
        help="with --factor, ignore pairs where both are below this, since a "
        "ratio between two numbers near zero means nothing (5). Does not "
        "apply to --difference or --sigma",
    )
    dg.add_argument(
        "--limit", type=int, default=0, help="keep only the N worst; 0 for all"
    )

    bm = sub.add_parser(
        "background-model",
        help="fit a model of the background -- R(s) . G(phi), smooth, times the "
        "polarisation, solid angle and sensor efficiency -- and write the table with "
        "the reflections whose background does not fit it flagged out of scaling",
    )
    bm.add_argument("expt", help="the integrated experiments")
    bm.add_argument("refl", help="the integrated reflections")
    bm.add_argument(
        "-o",
        "--output",
        default="filtered.refl",
        help="the filtered table (filtered.refl)",
    )
    bm.add_argument("--z-max", type=float, default=5.0, help="flag beyond this |z| (5)")
    bm.add_argument(
        "--knots", type=int, default=30, help="intervals of R's spline in 1/d (30)"
    )
    bm.add_argument(
        "--phi-spacing", type=float, default=10.0, help="degrees between G's knots (10)"
    )
    bm.add_argument(
        "--smoothness",
        type=float,
        default=100.0,
        help="the roughness penalty, in observations (100)",
    )
    bm.add_argument(
        "--reject-shells",
        type=float,
        default=None,
        metavar="SPREAD",
        help="also leave out whole resolution shells whose z spread exceeds this",
    )
    bm.add_argument(
        "--plot", default=None, help="draw the model against the data to this file"
    )
    bm.add_argument(
        "--zoom",
        type=float,
        default=250.0,
        help="pixels about the beam in the map (250)",
    )

    bk = sub.add_parser(
        "background",
        help="each reflection's background against resolution, azimuth, image and "
        "place on the detector, against its running median: for looking at before "
        "modelling it",
    )
    bk.add_argument("expt", help="the integrated experiments, for the beam centre")
    bk.add_argument("refl", help="the integrated reflections, mxi's or DIALS's")
    bk.add_argument(
        "-o", "--output", default="background.png", help="the figure (background.png)"
    )
    bk.add_argument(
        "--shells", type=int, default=20, help="resolution shells in the table (20)"
    )
    bk.add_argument(
        "--range",
        nargs=2,
        type=float,
        default=None,
        metavar=("DMAX", "DMIN"),
        help="limit the table, and the azimuth and image panels, to this resolution range",
    )
    bk.add_argument(
        "--zoom",
        type=float,
        default=250.0,
        help="pixels about the beam in the map (250)",
    )
    bk.add_argument(
        "--no-plot", action="store_true", help="the table only, without matplotlib"
    )

    ob = sub.add_parser(
        "observations",
        help="every observation of given reflections, all their symmetry "
        "equivalents, with what the table records of each",
    )
    ob.add_argument(
        "expt",
        help="the experiments, for the space group: a scaled one for all equivalents",
    )
    ob.add_argument("refl", help="the reflections, DIALS's or mxi's")
    ob.add_argument("hkl", nargs="+", help="reflections, as h,k,l")
    ob.add_argument(
        "--limit",
        type=int,
        default=0,
        help="at most this many observations each (0, all)",
    )

    un = sub.add_parser(
        "unique",
        help="two scaled data sets, DIALS's and mxi's say, compared one unique "
        "reflection at a time: I/sigma split into intensity, multiplicity and "
        "sigma, and the observations' scatter, by resolution shell",
    )
    un.add_argument("expt_a", help="the first data set's scaled experiments")
    un.add_argument("refl_a", help="its scaled reflections")
    un.add_argument("expt_b", help="the second's scaled experiments")
    un.add_argument("refl_b", help="its scaled reflections")
    un.add_argument(
        "--shells", type=int, default=10, help="resolution shells, of equal counts (10)"
    )
    un.add_argument(
        "--csv", default=None, help="also write every matched reflection to this file"
    )
    un.add_argument(
        "--labels", nargs=2, default=None, metavar=("A", "B"), help="names for the two"
    )

    fl = sub.add_parser(
        "failures",
        help="where profile fitting fails: an integrated table's profile.failure, "
        "counted by reason and binned by scan, zeta, box size, detector place",
    )
    fl.add_argument("refl", help="an integrated table from mxi_integrate")
    fl.add_argument(
        "--code", type=int, default=None, help="one reason only, by its code"
    )

    ce = sub.add_parser(
        "compare-expt",
        help="two experiment lists' models side by side: mxi_import's against "
        "dials.import's, say -- each line the same or by how much it differs",
    )
    ce.add_argument("first", help="an experiment list")
    ce.add_argument("second", help="another, to compare it with")

    eq = sub.add_parser(
        "equivalents",
        help="where integration is biased: each observation against its clean "
        "symmetry equivalents, binned by partiality, masked pixels and scan edges",
    )
    eq.add_argument(
        "expt",
        help="the scaled (or integrated) experiments, for the space group and scan",
    )
    eq.add_argument("refl", help="the scaled (or integrated) reflections")
    eq.add_argument(
        "--intensity",
        choices=["prf", "sum", "both"],
        default="both",
        help="profile fitted, summed, or both (both)",
    )
    eq.add_argument(
        "--min-i-sigma",
        type=float,
        default=5.0,
        help="the I/sigma an observation and its reference both need to be counted (5)",
    )
    eq.add_argument(
        "--reference",
        choices=["mean", "weighted"],
        default="mean",
        help="the clean equivalents' unweighted mean (mean), or weighted by "
        "1/variance, which reads unbiased data high where it is weak",
    )
    eq.add_argument(
        "--worst",
        type=int,
        default=0,
        help="also write the N observations furthest from their equivalents to "
        "equivalents_<kind>.refl, for dials.image_viewer",
    )
    eq.add_argument(
        "--all",
        action="store_true",
        help="count every observation in every table, not only those clean in "
        "every other respect",
    )
    eq.add_argument(
        "--anomalous",
        action="store_true",
        help="Friedel mates as different reflections, for data with a strong anomalous signal",
    )

    i = sub.add_parser("inspect", help="describe a file without assuming its layout")
    i.add_argument("path")

    return parser


def main(argv: list[str] | None = None) -> int:
    # Reports are long by design and `mxeq check ... | head` is the normal way
    # to read one. Python's default SIGPIPE handling turns that into a
    # traceback on stderr; restoring the Unix default makes it do nothing,
    # which is what every other command line tool does.
    try:
        import signal

        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except (AttributeError, ValueError):  # not POSIX, or not the main thread
        pass

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "background-model":
        from . import background as bg_tables
        from . import background_model

        result = background_model.run(
            args.expt,
            args.refl,
            args.z_max,
            args.reject_shells,
            knots_s=args.knots,
            phi_spacing=args.phi_spacing,
            smoothness=args.smoothness,
        )
        print(result.report)
        refl.write(args.output, result.table)
        print(f"\nWrote {args.output}")
        if args.plot:
            from .plots import background as background_plot

            centre = bg_tables.load(args.expt, args.refl).centre
            px = np.asarray(result.table.columns["xyzcal.px"]).reshape(-1, 3)
            background_plot.draw_model(
                result.fit,
                np.asarray(result.table.columns["background.mean"], float).ravel(),
                px[:, 0],
                px[:, 1],
                centre,
                args.plot,
                args.z_max,
                args.zoom,
                args.refl,
            )
            print(f"Wrote {args.plot}")
        return 0

    if args.command == "background":
        from . import background

        b = background.load(args.expt, args.refl)
        print(
            background.table(b, args.shells, tuple(args.range) if args.range else None)
        )
        if not args.no_plot:
            from .plots import background as background_plot

            background_plot.draw(
                b,
                args.output,
                tuple(args.range) if args.range else None,
                args.zoom,
                args.refl,
            )
            print(f"\nWrote {args.output}")
        return 0

    if args.command == "observations":
        from . import observations

        print(
            observations.listing(
                args.expt,
                args.refl,
                [observations.parse_index(h) for h in args.hkl],
                args.limit,
            )
        )
        return 0

    if args.command == "unique":
        from . import unique

        names = args.labels or [args.refl_a, args.refl_b]
        a = unique.load(args.expt_a, args.refl_a, names[0])
        b = unique.load(args.expt_b, args.refl_b, names[1])
        print(unique.compare(a, b, args.shells, args.csv))
        return 0

    if args.command == "failures":
        from . import failures

        print(failures.analyse(args.refl, args.code))
        return 0

    if args.command == "compare-expt":
        from . import compare_expt

        result = compare_expt.compare_files(args.first, args.second)
        print("\n".join(result.lines))
        return 0 if result.differences == 0 else 1

    if args.command == "equivalents":
        from . import equivalents

        kinds = ("prf", "sum") if args.intensity == "both" else (args.intensity,)
        results, prepared = equivalents.from_files(
            args.expt,
            args.refl,
            kinds,
            args.min_i_sigma,
            args.anomalous,
            not args.all,
            args.reference,
        )
        print(equivalents.report(results, args.min_i_sigma, not args.all, prepared))
        if args.worst > 0:
            for r in results:
                path = f"equivalents_{r.intensity}.refl"
                table = equivalents.worst(prepared, r, r.intensity, args.worst)
                refl.write(path, table)
                print(
                    f"Wrote the {table.nrows} {r.intensity} observations furthest from "
                    f"their equivalents to {path}, for dials.image_viewer"
                )
        return 0

    if args.command == "residuals":
        from . import residuals

        print(
            residuals.write(
                refl.load(args.path),
                args.output,
                least_signal=args.least_signal,
                n_bins=args.bins,
                cells=args.cells,
                title=args.title,
            )
        )
        return 0

    if args.command == "explain":
        from . import match as matching

        a, b = refl.load(args.a), refl.load(args.b)
        cols = lambda t: (
            t.columns["miller_index"],
            t.columns["entering"].ravel().astype(int),
            t.columns["xyzcal.px"][:, 2],
        )
        ia, ib, unpartnered = matching.match_observations(
            *cols(a), *cols(b), radius=args.radius
        )
        turn = matching.estimate_turn(*cols(a))
        print(
            f"{a.nrows} rows against {b.nrows}, matched {len(ia)}"
            f" ({100.0 * len(ia) / max(min(a.nrows, b.nrows), 1):.1f} per cent "
            "of the smaller)"
        )
        print(
            "one turn is "
            + (f"{turn:.1f} frames" if turn else "not measurable, less than a turn")
            + ", from the repeats in A"
        )
        for label, first, second, matched in (
            (args.a, a, b, ia),
            (args.b, b, a, ib),
        ):
            why = matching.explain_unpartnered(
                *cols(first),
                *cols(second),
                matched,
                args.radius,
                turn,
                sample=args.sample,
            )
            total = sum(why.values())
            print()
            print(f"unpartnered in {label}, a sample of {total}:")
            for reason, n in why.items():
                if n:
                    print(f"  {100.0 * n / max(total, 1):5.1f}%  {reason}")
        return 0

    if args.command == "html":
        from . import htmlreport

        values = args.value or [
            "intensity.sum.value",
            "intensity.prf.value",
            "background.mean",
        ]
        print(
            htmlreport.write(
                refl.load(args.a),
                refl.load(args.b),
                args.output,
                values,
                n_bins=args.bins,
                radius=args.radius,
                title=args.title,
                plotly_src=args.plotly or htmlreport.PLOTLY_CDN,
            )
        )
        return 0

    if args.command == "disagree":
        from . import disagree

        table, report = disagree.select(
            refl.load(args.ours),
            refl.load(args.theirs),
            value=args.value,
            value_b=args.value_b,
            factor=args.factor,
            difference=args.difference,
            sigmas=args.sigma,
            floor=args.floor,
            limit=args.limit,
        )
        print(report)
        if table.nrows == 0:
            print("nothing to write")
            return 0
        refl.write(args.output, table)
        print(f"wrote {args.output}")
        print(f"  dials.image_viewer imported.expt {args.output}")
        return 0

    if args.command == "profiles":
        from .plots import profiles as profile_plots

        reference = profile_plots.read_profiles(args.path)
        print(
            f"{len(reference.profiles)} profiles, {reference.side} a side, "
            f"{reference.divisions}x{reference.divisions} regions per panel"
        )
        thin = [r for r, n in enumerate(reference.spots) if n < 10]
        if thin:
            print(f"  {len(thin)} cells saw fewer than ten spots: {thin[:12]}")
        print(
            f"  {reference.blocks} scan blocks, "
            f"{reference.divisions}x{reference.divisions} cells a panel"
        )
        print(profile_plots.draw_profiles(reference, args.output, args.block))
        if args.widths:
            print(profile_plots.draw_profile_widths(reference, args.widths))
        return 0

    if args.command == "trend":
        from . import trends

        values = args.value or [
            "intensity.sum.value",
            "intensity.prf.value",
            "intensity.sum.variance",
            "intensity.prf.variance",
            "background.mean",
        ]
        print(
            trends.compare(
                refl.load(args.a),
                refl.load(args.b),
                values,
                value_b=args.value_b,
                n_bins=args.bins,
                radius=args.radius,
                worst=args.worst,
            )
        )
        return 0

    if args.command == "inspect":
        with open(args.path, "rb") as f:
            raw = f.read()
        if raw[:1] in (b"{", b"\n", b" "):
            experiments = expt.load(args.path)
            print(f"experiment list, {len(experiments)} experiments")
            for n, e in enumerate(experiments):
                print(f"  [{n}] identifier {e.identifier!r}")
                if e.crystal:
                    cell = " ".join(f"{v:.4f}" for v in e.crystal.cell)
                    print(f"      cell {cell}")
                    print(f"      hall {e.crystal.hall!r}")
                    print(f"      scan varying: {e.crystal.scan_varying}")
                if e.scan:
                    print(f"      images {e.scan.image_range} osc {e.scan.oscillation}")
            return 0
        try:
            table = refl.loads(raw)
        except refl.ReflFormatError as exc:
            print(f"could not read as a reflection table: {exc}\n", file=sys.stderr)
            print("\n".join(refl.structure(raw)))
            return 2
        print(f"reflection table, {table.nrows} rows")
        print(f"identifiers: {table.identifiers}")
        for name in sorted(table.columns):
            values = table.columns[name]
            shape = "" if values.ndim == 1 else f" x{values.shape[1]}"
            print(f"  {name:<34} {table.types[name]}{shape}")
        for name, type_name in sorted(table.opaque.items()):
            print(f"  {name:<34} {type_name}  (not decoded)")
        return 0

    # Which kind of file each boundary wants, checked before anything tries to
    # parse one as the other. `check refined` compares experiment lists and the
    # rest compare reflection tables, and passing the wrong pair used to fail
    # inside a UTF-8 decoder with the offending byte's position -- a message
    # that says nothing about what went wrong or what to do.
    # `check refined` compares models, and models live in .expt. But comparing
    # two refined pipelines usually means asking how their residuals compare,
    # and that lives in the .refl -- so a pair of reflection tables is accepted
    # and routed to the comparison that answers it, rather than refused on a
    # technicality.
    kinds = {_sniff(args.a), _sniff(args.b)}
    if args.boundary == "refined" and kinds == {"refl"}:
        print(
            "mxeq: comparing reflections, since both files are tables; "
            "pass the .expt files to compare the models instead",
            file=sys.stderr,
        )
        args.boundary = "indexed"

    for path in (args.a, args.b):
        kind = _sniff(path)
        wanted = "expt" if args.boundary == "refined" else "refl"
        if kind is not None and kind != wanted:
            print(
                f"mxeq: {path} looks like a{'n' if kind == 'expt' else ''} "
                f"{'experiment list (.expt)' if kind == 'expt' else 'reflection table (.refl)'}, "
                f"and 'check {args.boundary}' compares "
                f"{'experiment lists' if wanted == 'expt' else 'reflection tables'}.",
                file=sys.stderr,
            )
            if wanted == "expt":
                print(
                    "       'check refined' compares the models -- cell, "
                    "orientation, detector -- so it wants the .expt files.",
                    file=sys.stderr,
                )
            else:
                print(
                    f"       'check {args.boundary}' compares reflections, so "
                    "it wants the .refl files; the .expt goes in -e.",
                    file=sys.stderr,
                )
            return 2

    boundary = args.boundary
    experiments_a = _load_expt(args.expt)
    experiments_b = _load_expt(args.expt_b or args.expt)
    operator = _operator(args.operator)

    if boundary == "refined":
        report = refined.check(expt.load(args.a), expt.load(args.b), args.a, args.b)
    else:
        table_a, table_b = refl.load(args.a), refl.load(args.b)
        if boundary == "auto":
            boundary = _guess(table_a)
            print(
                f"mxeq: guessing boundary {boundary!r} from A's columns",
                file=sys.stderr,
            )
        if boundary == "strong":
            report = strong.check(table_a, table_b, args.a, args.b, radius=args.radius)
        elif boundary == "indexed":
            report = indexed.check(
                table_a, table_b, args.a, args.b, experiments_a, radius=args.radius
            )
        elif boundary == "integrated":
            report = integrated.check(
                table_a,
                table_b,
                args.a,
                args.b,
                experiments_a,
                operator=operator,
                n_bins=args.bins,
            )
        else:
            report = scaled.check(
                table_a,
                table_b,
                args.a,
                args.b,
                experiments_a or experiments_b,
                operator=operator,
                n_bins=args.bins,
            )

    print(report.json() if args.json else report.text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
