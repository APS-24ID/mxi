"""Each step's notes, for learning what it does: what it does, what to look
at in its Summary and Plots, what goes wrong and which option helps, and the
guides that say more (on GitHub: the release archives carry no documents).
Every option named here is one its command takes -- tests/test_ui.py checks
each against the command's --help, as the forms' are; another command's
written with it, "mxi scale --space-group", and checked against that one's. In mxi's terms, from
its documents (docs/algorithms.md and those linked), not DIALS's.
"""

from dataclasses import dataclass

DOCS = "https://github.com/APS-24ID/mxi/blob/master/"


@dataclass(frozen=True)
class Note:
    does: str
    look: str
    wrong: str
    read: tuple  # (label, path in the repository)


NOTES = {
    "import": Note(
        "Reads each NXmx master file -- the beam, the detector's panels, the "
        "goniometer and the scan -- into an experiment list, one experiment a "
        "sweep, with the detector's own pixel mask. Nothing is measured yet: "
        "everything after depends on this geometry being right.",
        "The Summary's sweeps: the images, the oscillation, the wavelength, the "
        "distance and the beam centre each file gave. View images: the images "
        "themselves, and whether anything shadows the detector -- a backstop, "
        "its arm, a cryo nozzle -- to mark untrusted before finding spots.",
        "A file without a value mxi needs says so: give it with --wavelength, "
        "--distance or --beam-centre. A beam centre far out is the commonest "
        "reason indexing fails later.",
        (("Importing", "docs/import.md"), ("Masking the detector", "docs/mask.md")),
    ),
    "find": Note(
        "Finds the spots on every image by the extended dispersion threshold, "
        "as DIALS's default: a pixel is signal where its neighbourhood's "
        "dispersion and mean say so; signal pixels joined across images into "
        "spots, each with its centroid and box. On the CPU, or with --gpu on a "
        "GPU.",
        "Spots per image: a steady count through the sweep, falling only as the "
        "crystal turns or decays; a sudden drop is the crystal leaving the beam "
        "or a shutter. Spots by resolution: how far the diffraction reaches. "
        "View images: the boxes on the spots, and none on the background.",
        "Ice rings show as spikes in spots by resolution. Many tiny spots on "
        "noise: raise --min-spot-size. Streaked or split spots thrown out: "
        "--max-separation, how far a spot's peak may lie from its centroid "
        "(2 pixels; 0 turns the check off).",
        (("The spot finder", "docs/spots.md"),),
    ),
    "index": Note(
        "Finds the crystal's lattice: the spots carried into reciprocal space, "
        "the 3D FFT for candidate basis vectors, the basis chosen as dials.index "
        "chooses, then assignment and refinement in turn. With several sweeps "
        "of one crystal they are indexed together; of several crystals, each "
        "apart (--no-joint).",
        "The indexed fraction -- most of the spots, for a single good crystal; "
        "the cell, and whether it is the one expected; the RMSDs cycle by cycle "
        "falling to well under a pixel and a fraction of an image. The "
        "candidate cells: whether another came close.",
        "No lattice: check the beam centre first, then give --max-cell, or "
        "--d-min to leave out the weak outer spots. The cell known: give it "
        "with --known-cell, and --known-space-group for its centring. A second "
        "lattice left over: --max-lattices 2.",
        (
            ("Methods, indexing", "docs/algorithms.md"),
            ("Several lattices", "docs/multi-lattice.md"),
        ),
    ),
    "refine": Note(
        "Refines the crystal, the detector and the beam against the spots' "
        "positions, x, y and rotation angle, by weighted least squares -- first "
        "one crystal for the whole scan, then the crystal varying smoothly "
        "through it, control points every 10 degrees.",
        "The RMSDs, falling cycle by cycle; the cell over the scan, smooth, "
        "changing by little -- a jump is a fault in the data, a slow drift the "
        "crystal changing in the beam.",
        "A cell wandering on sparse data: --interval 20 or 30, control points "
        "further apart, or --static, one crystal for the scan. "
        "A crystal off the rotation axis -- spots predicted wrong as it turns: "
        "--sample-offset. Too many outliers rejected: --outlier-sigma.",
        (
            ("Methods, refinement", "docs/algorithms.md"),
            ("Several sweeps", "docs/multi-sweep.md"),
        ),
    ),
    "integrate": Note(
        "Measures every reflection the lattice predicts: its profile modelled "
        "from the strong spots, a box and foreground for each, the background "
        "by a robust Poisson fit, then the intensity summed and fitted against "
        "reference profiles. On the CPU, or with --gpu on a GPU.",
        "<I/σ> by resolution: where it falls to about 1 is near the data's "
        "limit. Profile correlation, high where the profiles fit. With several "
        "sweeps, a line each: a sweep unlike the others is worth a look.",
        "Integrating far past the diffraction costs time and adds noise: "
        "--d-min. Profiles that fit badly -- a split or odd-shaped crystal -- "
        "show as low profile correlation; --summation-only uses the summed "
        "intensities alone.",
        (("Integration", "docs/integration.md"),),
    ),
    "symmetry": Note(
        "Finds the Laue group and space group from the intensities: every "
        "symmetry element of the lattice scored by how well the reflections it "
        "relates agree, the subgroups weighed by likelihood, then the space "
        "group by its absences, and the data reindexed into it. Several crystals "
        "are first brought to one setting by cosym.",
        "The Laue group and space group, and how far ahead the first subgroup's "
        "likelihood is of the next: a close second is worth knowing about. The "
        "reindexing applied, if any.",
        "A lattice nearly of a higher symmetry, taken for it or not: --max-delta, "
        "the obliquity in degrees within which the higher symmetry is allowed "
        "(2). A group known: give it to scaling, mxi scale --space-group.",
        (("Symmetry", "docs/symmetry.md"),),
    ),
    "scale": Note(
        "Puts every observation on one scale -- a smooth scale and decay "
        "through the scan, and absorption -- and merges the symmetry "
        "equivalents, rejecting outliers and fitting an error model, as "
        "dials.scale's physical model does.",
        "The merging statistics: CC½ by resolution, which says how far the data "
        "are good -- with --d-min-auto, the limit from CC½ itself; completeness "
        "and multiplicity; <I/σ>; the R factors, rising outwards. The error "
        "model's a and b.",
        "Anomalous data: --anomalous, so that Friedel mates are kept apart. A "
        "bad stretch of images, or a crystal's whole sweep: --exclude "
        "E:FIRST:LAST, experiment E's images left out. Strong absorption: "
        "--absorption-level high.",
        (("Scaling", "docs/scaling.md"), ("What each step reports", "docs/output.md")),
    ),
    "export": Note(
        "Writes the scaled data as an unmerged MTZ, as dials.export writes it, "
        "for merging and phasing elsewhere: the intensities with their "
        "corrections, a batch an image.",
        "The file, scaled.mtz, beside the others.",
        "--crystal-name and --project-name set the names programs downstream "
        "show; --d-min cuts the resolution.",
        (("Export", "docs/export.md"),),
    ),
    "verify": Note(
        "Checks a lattice against the spots: how much of what the indexed "
        "lattice predicts was actually seen, shell by shell -- a lattice that "
        "explains the spots predicts few that are not there.",
        "The fraction of predictions seen in each shell; low throughout is a "
        "wrong lattice, low only outwards the diffraction's limit.",
        "--d-min limits it to where the diffraction is.",
        (("Methods", "docs/algorithms.md"),),
    ),
    "correlation": Note(
        "Clusters data sets -- crystals, or sweeps of them -- by how alike their "
        "intensities are, as dials.correlation_matrix: the correlation of every "
        "pair, cosym's coordinates, and OPTICS to find the clusters. For many "
        "crystals, to find which belong together before merging them.",
        "The matrix in OPTICS's order: its blocks are the clusters. The "
        "reachability, low within a cluster and high between. The principal "
        "components, each cluster a group apart.",
        "Clusters too fine or too coarse: --min-samples. To scale each cluster "
        "apart: --write-clusters, then mxi process on each cluster's file.",
        (
            ("Correlation", "docs/correlation.md"),
            ("Many crystals: a tutorial", "docs/tutorials/cows-pigs-people.md"),
        ),
    ),
    "report": Note(
        "Writes an HTML report of the whole processing, from the records each "
        "step left in the file -- every step's figures and plots, to read in a "
        "browser or send.",
        "mxi_report.html, opened in a browser.",
        "A step that left no record shows no section: run it with this mxi.",
        (("What each step reports", "docs/output.md"),),
    ),
}


def note_for(key):
    return NOTES.get(key)


def link(path):
    return DOCS + path
