"""The workflow's steps, in order, and the tools beside them -- the main
sequence as the README gives it, each command's default input and output,
and its main options as its form's fields. Every option here is checked
against the command's own --help (tests/test_ui.py), so a form cannot offer
what mxi does not take."""

from mxi_ui.model import Option, Step

D_MIN = Option("--d-min", "Resolution limit (Å)", "number", "nothing beyond this d")

STEPS = (
    Step(
        "import",
        "Import",
        "the images' experiment list, a sweep for each master file",
        (),
        "imported.rflx",
        (
            Option(
                "--wavelength",
                "Wavelength (Å)",
                "number",
                "for a file without one, or a wrong one",
            ),
            Option(
                "--distance",
                "Distance (mm)",
                "number",
                "the detector distance, if the file's is wrong",
            ),
            Option(
                "--beam-centre",
                "Beam centre (x,y px)",
                "text",
                "where the beam meets the detector",
            ),
            Option(
                "--image-range",
                "Images (first,last)",
                "text",
                "only these, counting from 1",
            ),
            Option(
                "--split",
                "Split each sweep into",
                "number",
                "sweeps of consecutive images",
            ),
        ),
    ),
    Step(
        "find",
        "Find spots",
        "strong spots on every image",
        ("imported.rflx",),
        "strong.rflx",
        (
            Option(
                "--gpu",
                "On the GPU",
                "flag",
                "the threshold on the GPU, where the build has one",
            ),
            Option(
                "--min-spot-size",
                "Smallest spot (pixels)",
                "number",
                "contiguous pixels a spot needs (3)",
            ),
            Option(
                "--max-spot-size",
                "Largest spot (pixels)",
                "number",
                "the most a spot may have (1000)",
            ),
            Option(
                "--max-separation",
                "Peak to centroid (px)",
                "number",
                "0 turns the test off (2)",
            ),
        ),
    ),
    Step(
        "index",
        "Index",
        "the lattice, or lattices, that explain the spots",
        ("strong.rflx",),
        "indexed.rflx",
        (
            D_MIN,
            Option(
                "--max-cell",
                "Longest cell edge (Å)",
                "number",
                "default: from the spots' spacing",
            ),
            Option("--known-cell", "Known cell", "text", "a,b,c,alpha,beta,gamma"),
            Option(
                "--known-space-group",
                "Known space group",
                "text",
                "the known cell's setting: P 21 21 21",
            ),
            Option(
                "--max-lattices",
                "Lattices, at most",
                "number",
                "several crystals in one beam",
            ),
            Option(
                "--no-joint",
                "Each sweep on its own",
                "flag",
                "several sweeps, each a crystal",
            ),
        ),
    ),
    Step(
        "refine",
        "Refine",
        "beam, detector and crystal against the spots",
        ("indexed.rflx",),
        "refined.rflx",
        (
            Option(
                "--static",
                "Static crystal",
                "flag",
                "one orientation for the scan, not scan-varying",
            ),
            Option(
                "--sample-offset",
                "Crystal off the axis",
                "flag",
                "refine its offset from the rotation axis",
            ),
            Option("--no-beam", "Hold the beam", "flag", "its direction fixed"),
            Option("--no-detector", "Hold the detector", "flag", "its position fixed"),
            Option(
                "--outlier-sigma",
                "Outliers beyond (σ)",
                "number",
                "4; 0 keeps every reflection",
            ),
        ),
    ),
    Step(
        "integrate",
        "Integrate",
        "every reflection's intensity, summed and profile-fitted",
        ("refined.rflx",),
        "integrated.rflx",
        (
            D_MIN,
            Option(
                "--gpu",
                "On the GPU",
                "flag",
                "profile fitting on the GPU, where the build has one",
            ),
            Option("--summation-only", "Summation only", "flag", "no profile fitting"),
            Option(
                "--min-zeta",
                "Nearest the axis (zeta)",
                "number",
                "skip reflections nearer than this",
            ),
        ),
    ),
    Step(
        "symmetry",
        "Symmetry",
        "the Laue group and space group, the data reindexed",
        ("integrated.rflx",),
        "symmetrized.rflx",
        (
            Option(
                "--max-delta",
                "Obliquity (degrees)",
                "number",
                "the lattice's symmetry to this (2)",
            ),
            Option(
                "--no-cosym",
                "Several crystals as they are",
                "flag",
                "pooled, not first brought to one setting",
            ),
        ),
    ),
    Step(
        "scale",
        "Scale",
        "scaled intensities, with merging statistics",
        ("symmetrized.rflx",),
        "scaled.rflx",
        (
            Option("--space-group", "Space group", "text", "default: symmetry's"),
            Option("--anomalous", "Anomalous", "flag", "Friedel mates apart"),
            D_MIN,
            Option(
                "--d-min-auto",
                "Limit from CC½",
                "flag",
                "scale everything, then set the limit from CC½",
            ),
            Option(
                "--absorption-level",
                "Absorption",
                "choice",
                "as dials.scale's",
                ("low", "medium", "high"),
            ),
            Option(
                "--exclude", "Leave out", "text", "experiment E's images: E:FIRST:LAST"
            ),
        ),
    ),
    Step(
        "export",
        "Export",
        "an unmerged MTZ, for the programs after",
        ("scaled.rflx",),
        "scaled.mtz",
        (
            D_MIN,
            Option("--crystal-name", "Crystal name", "text", "(XTAL)"),
            Option("--project-name", "Project name", "text", "(mxi)"),
        ),
    ),
)

TOOLS = (
    Step(
        "verify",
        "Verify a lattice",
        "how much of what each lattice predicts was seen",
        ("indexed.rflx",),
        "",
        (
            D_MIN,
            Option("--shells", "Shells", "number", "of equal reciprocal volume (10)"),
            Option(
                "--no-mask",
                "Ignore the mask",
                "flag",
                "count masked pixels as observable",
            ),
        ),
        output_flag="",
        tool=True,
    ),
    Step(
        "correlation",
        "Correlation",
        "crystals clustered by how alike their data are",
        ("symmetrized.rflx",),
        "",
        (
            Option(
                "--write-clusters",
                "Write each cluster",
                "flag",
                "as its own .rflx, to scale apart",
            ),
            Option("--prefix", "Their names", "text", "PREFIX_0.rflx, ... (cluster)"),
            Option("--min-samples", "OPTICS min_samples", "number", "not searched for"),
        ),
        output_flag="",
        tool=True,
    ),
    Step(
        "report",
        "Report",
        "an HTML report of the processing",
        ("scaled.rflx",),
        "mxi_report.html",
        tool=True,
    ),
)


def by_key(key):
    for step in STEPS + TOOLS:
        if step.key == key:
            return step
    raise KeyError(key)
