# Several lattices in one data set: a plan

STATUS: a plan for review, October 2026. Nothing here is built.

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

**Symmetry** (mxi_symmetry). Each lattice's Laue group from its own data, then
one setting for them all: where lattices are of one crystal -- a split crystal,
a twin -- or of one kind, their indexing must agree, and an indexing ambiguity
resolved between them as dials.cosym resolves it between data sets, by
correlating their intensities under each of the ambiguity's operators.

**Scaling** (mxi_scale). A scale model a lattice, as now a sweep; mxi_export a
data set a lattice, with its own batches.

## The plan

Each stage checked before the next, and none changing any result on data with
one lattice: insulin, the four sweeps, the small molecule and Graeme's P 4_1
crystal give the bytes they give now, stage by stage.

0. **Data.** A planted data set: two lattices predicted onto one sweep, of known
   orientations, cells and profiles -- including a pair whose spots overlap.
   And real data with several lattices, processed by DIALS to compare with
   (Graeme: which?).
1. **Shared models.** The experiment list's sharing kept, read and written by
   index; a profile model a lattice; the two views, sweeps and lattices.
   Refinement's parameters a block a model. On every data set now, the same
   results -- the four-sweep data excepted if separate sweeps are then to share
   a detector, as DIALS shares one (question 4).
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
6. **Symmetry and scaling** of several lattices: one setting, ambiguities
   resolved; on the real data, against dials.cosym and dials.scale.
7. The documents.

## Questions for Graeme

1. **#1029 in mxi's own files**: carry the sweeps through as crystal-less
   experiments, as the issue proposes, with mxi_convert dropping them for
   DIALS? Or DIALS's form in mxi's files too?
2. **Data**: a real data set with several lattices, and DIALS's processing of
   it, to judge against.
3. **What first**: two crystals in the beam, a split crystal, or a twin that is
   not merohedral? They differ most in symmetry: lattices of one crystal must
   share a setting, of different crystals need not.
4. **Several sweeps, one detector**: should separate sweeps share the detector
   in refinement, as DIALS shares one model by default, once refinement can share
   at all? It would change the four-sweep results, and is a choice for itself.
5. **Stills**: out of scope, as mxi is rotation data; but the same model -- a
   sweep of one image, lattices in it -- would serve.
