# Several lattices in one data set: a plan

STATUS: agreed, October 2026, with Graeme's answers (below). Built:

* Stage 1, in part: models shared in refinement as the experiment list numbers
  them -- ExperimentList::beam_of and detector_of, read from the file; a
  parameter block a distinct detector and a distinct beam, owned by all that
  share it; written back interned, once. A single sweep refines to the same
  bytes as before. The profile model a lattice is moved to stage 5, where
  integration first uses it.
* Stage 3, before stage 2: mxi_index --max-lattices N (index_lattices in the
  library). On Graeme's ins2_27: lattice 1 indexes 15307 of 21299 spots,
  lattice 2 5906 of the 5992 left, no third in the last 86 -- 99.6 per cent --
  refined together, one detector, to 0.258 and 0.264 px and 0.227 images; one
  detector, one beam, one imageset written, two crystals. mxi_refine takes the
  file as it is, the detector shared as the file numbers it. Unindexed spots
  keep mxi's convention, Miller index 0 0 0 and `id` 0, until stage 2 carries
  the sweep's own experiment through -- done after, on the multi-lattice output
  it must handle. Several sweeps with several lattices: refused for now.
* Stage 0's spots planted in the tests: two lattices 11.6 degrees apart on one
  sweep, found and refined; its images, for integration, at stage 5.

## Decisions (Graeme, October 2026)

1. **Sweeps carried through**, as dials/dials#1029 proposes, in mxi's own files;
   the crystal-less experiments dropped only when exporting for DIALS.
2. **Data**: a real data set, two crystals of I23 symmetry, with DIALS's
   processing of it to judge against. I23 in an I432 lattice has a twofold
   indexing ambiguity, so the lattices' settings need a cosym-type analysis.
3. **First case**: two crystals in the beam.
4. **Models shared where they are the same model**, as the experiment list
   records them: a sweep's lattices share its detector, beam and goniometer;
   separate sweeps share a detector where it stood still between them -- as in
   some of Graeme's data -- and not where it moved, as in the four sweeps,
   whose results so do not change.

dials.index finds several lattices in one data set (`max_lattices`): a split
crystal, a twin that is not merohedral, two crystals in the beam. Each lattice
has its own orientation and cell -- its own crystal model, its own matrix --
and its own spots, with their own mosaicity and so their own profile; the
lattices share everything else, because they share the images: the beam, the
detector, the goniometer and the scan. mxi should do the same, and should do
it without the problem dials/dials#1029 describes.

## What DIALS does, and the problem with it

An experiment in DIALS is a set of images with every model needed to predict
from them, the crystal included. Indexing turns each experiment it is given
into one experiment per lattice found: one lattice, one out; two, two out,
sharing the imageset and the scan; none, none out -- the experiment is deleted,
and with it any way to look at its images and the spots that did not index.
Unindexed spots keep `id` -1, which names no experiment. Graeme's issue,
dials/dials#1029 (2019, open): carry every experiment given through
indexing, and add one for each lattice found, so that the count out is always
the count in plus the lattices, and nothing given is lost. The objection in
the issue is the cost in DIALS, where the assumption that every experiment has
a crystal is everywhere -- `if exp.crystal is None: continue`, sprinkled.

mxi can be built so that it is never sprinkled: the distinction made once, where
experiment lists are read, and nowhere else.

## The model

Two things DIALS calls by one name:

* **A sweep**: images, and the models they need -- beam, detector, goniometer,
  scan, the imageset. Every sweep given is carried through every program,
  whether anything is found in it or not.
* **A lattice**: a crystal model and the profile model of its spots, in one
  sweep. A sweep has none, one or several.

In a file -- an `.expt`, or a `.rflx`'s `/experiments` -- the experiment list is
DIALS's, unchanged in form, so that DIALS and dxtbx read it:

* an experiment for each lattice, referring by index to its sweep's beam,
  detector, goniometer, scan and imageset -- shared, as DIALS shares them, not
  copied -- and to its own crystal and profile;
* and, as #1029 proposes, each sweep's own experiment, crystal-less, carried
  through from import: the experiment its unindexed spots belong to.

In a reflection table, `id` is the experiment -- a lattice's, or for an
unindexed spot its sweep's crystal-less one, never -1 for a spot that has a
sweep -- and `imageset_id` its sweep, always.

In memory, the experiment list says which experiments share which model,
rather than each experiment holding its own copy as now; and holds a profile
model for each lattice, where it now holds one for the whole list. Every
program that predicts takes the lattices; every program that reads images
takes the sweeps -- two views of the one list, made in one place, so that no
program asks of an experiment whether it has a crystal.

For DIALS, mxi_convert writes what DIALS expects: the crystal-less experiments
dropped, the ids renumbered, unindexed spots -1 again -- that conversion, too,
in one place.

## The programs

**Indexing** (mxi_index `--max-lattices N`, default 1). As dials.index:
index; refine; take the spots the lattice indexes out; index what remains;
until N lattices or too few spots are left. A lattice that is the same as one
found -- the same cell and orientation to within the tolerance, or related to
it by the lattice's symmetry -- is not new, and stops the search. A spot that
two lattices both predict within the tolerance goes to the nearer prediction,
and is flagged; a minimum number of spots makes a lattice. Each found lattice
refined with those before it, between cycles, as dials.index does.

**Refinement** (mxi_refine). One parameter block for each distinct model, not
each experiment: the lattices of a sweep share one detector, one beam, one
goniometer and so one set of their parameters; each lattice its own crystal,
static or scan-varying, with the sweep's control points. Now each experiment
gets its own six detector and two beam parameters -- fine for separate sweeps,
wrong for lattices on the same images, which would refine as many detectors as
lattices. Outliers judged within each lattice, as now within each sweep.

**Integration** (mxi_integrate). The lattices of a sweep integrated together,
in one pass over its images -- not each in a run of its own, as several sweeps
are now, which would read every image once a lattice and leave each lattice
blind to the others' spots. Every lattice predicts; where one reflection's
foreground lies over another's -- of the same lattice or another -- those
pixels are no one's background, and both reflections are flagged as
overlapped and kept out of the reference profiles. Each lattice estimates its
own profile model, sigma_b and sigma_m, from its own strong spots, and learns
its own reference profiles: a split crystal's two halves need not be equally
mosaic.

**Symmetry** (mxi_symmetry). Lattices of different crystals -- the first case,
two crystals in the beam -- are separate data sets, each indexed in a basis of
its own choosing among those its lattice's symmetry allows: before they can be
scaled together they must be brought to one setting, and their Patterson group
found from them all. Gildea & Winter (2018), Acta Cryst. D74, 405-410, which
dials.cosym implements, does both at once:

* The lattice group: the highest symmetry the cells allow, its m operations --
  I432's 24, for I23.
* Each data set i under each operation k, its intensities merged in P1 and
  reindexed by it; the correlation r_ik,jl of each such with every other over
  the reflections they share -- a matrix (n m) by (n m), n data sets.
* Each data set under each operation a point x_ik in m dimensions, placed by
  minimising the sum over every pair of (r_ik,jl - x_ik . x_jl)^2, by L-BFGS
  from random starting points in 0 to 1.
* The points' principal components: they fall into as many clusters as the
  lattice group has cosets of the true Patterson group -- two for I23 in I432
  -- and the share of the variance on each component says how many.
* Within one cluster, the operations S_k^-1 S_l between its points are the
  Patterson group; a data set whose point in that cluster carries an operation
  S_k outside the group is reindexed by it -- the twinning operation -- taken
  to the unique representative of its coset of the group (Flack, 1987).

With two data sets the matrix is small and the answer could be had by trying
the ambiguity's operators directly; built as the paper describes, it serves
any number of crystals, and a multi-crystal data set is the next case after
two lattices. The Patterson group so found is then tested, and the space group
chosen by absences, by what mxi_symmetry already does for one data set.

**Scaling** (mxi_scale). A scale model a lattice, as now a sweep; mxi_export a
data set a lattice, with its own batches.

## The plan

Each stage checked before the next, and none changing any result on data with
one lattice: insulin, the four sweeps, the small molecule and Graeme's P 4_1
crystal give the bytes they give now, stage by stage.

0. **Data.** A planted data set: two lattices predicted onto one sweep, of known
   orientations, cells and profiles -- including a pair whose spots overlap.
   And Graeme's two crystals of I23 (ins2_27, October 2026): 600 images of
   0.1 degree on an Eiger 4M at 170 mm, 21299 spots. One lattice indexes 71.9
   per cent of them, insulin's primitive cell, 67.47 A and 109.47 degrees; the
   rest, indexed alone, 98.6 per cent of theirs to a second, 67.51 A, the two
   11.6 degrees apart under cubic symmetry -- two crystals, as said. mxi
   marks a spot unindexed by its Miller index and flags, `id` staying 0: stage
   2's crystal-less experiment is where those rows' `id` will point.
1. **Shared models.** The experiment list's sharing kept, read and written by
   index; a profile model a lattice; the two views, sweeps and lattices.
   Refinement's parameters a block a model, separate sweeps sharing theirs as
   DIALS shares them (decision 4). On every data set of one sweep, the same
   results; on the four sweeps, a change, set against dials.refine's -- which
   shares -- and judged on its own.
2. **Carrying sweeps through** (#1029): the crystal-less experiments, `id` and
   `imageset_id`, every program reading the list through the two views;
   mxi_convert's DIALS form. A sweep that indexes nothing still there, its
   spots on its images.
3. **Indexing several lattices.** On the planted data, both lattices found, to
   their true cells and orientations; one lattice not found twice; on the real
   data, the lattices dials.index finds.
4. **Refining them**, shared detector and beam, separate crystals: on the
   planted data, the true models; on the real, dials.refine's.
5. **Integrating them** together, overlaps flagged, a profile model each: on the
   planted data, each lattice's intensities as they were planted, overlapped
   spots flagged and the rest unharmed; on the real, against dials.integrate.
6. **Symmetry** of several lattices, cosym as above: on planted data -- two data
   sets of an I23 crystal, one reindexed by the I432 operation I23 lacks, and
   then many partial ones, as the paper's examples -- the Patterson group and
   each data set's operation recovered; then on Graeme's two crystals, against
   dials.cosym.
7. **Scaling and export** of several lattices: a scale model and a data set a
   lattice; against dials.scale.
8. The documents.

## Questions settled

Graeme's answers to the plan's questions are the decisions above, and his
two-crystal data is here. Still to come: DIALS's processing of it, to judge
against.
