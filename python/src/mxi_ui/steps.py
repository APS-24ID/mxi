"""The workflow's steps, in order, and the tools beside them -- the main
sequence as the README gives it, each command's default input and output."""

from mxi_ui.model import Step

STEPS = (
    Step(
        "import",
        "Import",
        "the images' experiment list, a sweep for each master file",
        (),
        "imported.rflx",
    ),
    Step(
        "find",
        "Find spots",
        "strong spots on every image",
        ("imported.rflx",),
        "strong.rflx",
    ),
    Step(
        "index",
        "Index",
        "the lattice, or lattices, that explain the spots",
        ("strong.rflx",),
        "indexed.rflx",
    ),
    Step(
        "refine",
        "Refine",
        "beam, detector and crystal against the spots",
        ("indexed.rflx",),
        "refined.rflx",
    ),
    Step(
        "integrate",
        "Integrate",
        "every reflection's intensity, summed and profile-fitted",
        ("refined.rflx",),
        "integrated.rflx",
    ),
    Step(
        "symmetry",
        "Symmetry",
        "the Laue group and space group, the data reindexed",
        ("integrated.rflx",),
        "symmetrized.rflx",
    ),
    Step(
        "scale",
        "Scale",
        "scaled intensities, with merging statistics",
        ("symmetrized.rflx",),
        "scaled.rflx",
    ),
    Step(
        "export",
        "Export",
        "an unmerged MTZ, for the programs after",
        ("scaled.rflx",),
        "scaled.mtz",
    ),
)

TOOLS = (
    Step(
        "verify",
        "Verify a lattice",
        "how much of what each lattice predicts was seen",
        ("indexed.rflx",),
        tool=True,
    ),
    Step(
        "correlation",
        "Correlation",
        "crystals clustered by how alike their data are",
        ("symmetrized.rflx",),
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
