# More than one sweep

STATUS: 6 October 2026; all six steps built, and working well on a genuine
four-sweep data set at different orientations: cubic insulin from Diamond's
i04, four full turns of 3600 images each (https://zenodo.org/records/8376818 --
ins10_1 is the first of them). There the crystal shared through refinement had
left the RMSDs high; with a crystal a sweep, refined apart after indexing with
one matrix, as DIALS does, it works well (Graeme). DIALS processes several sweeps of
one crystal together -- imported as one experiment each, indexed jointly to one
crystal, refined, integrated sweep by sweep, and symmetry and scaling over all
of them -- and mxi should, from `mxi_import` to `mxi_scale`.

## Where each program stands

| program | more than one sweep |
| --- | --- |
| `mxi_import` | built: several masters, one experiment each (step 1) |
| `mxi_find` | built: each experiment from its own images, one table with `id`s (step 2) |
| `mxi_index` | joint indexing with one matrix, then each sweep's crystal refined apart (step 3) |
| `mxi_refine` | several experiments, a crystal each by default |
| `mxi_integrate` | built: each sweep alone, its own profile model, one table with `id`s (step 4) |
| `mxi_symmetry` | built: the sweeps pooled, every experiment reindexed alike (step 5) |
| `mxi_scale` | built: a model a sweep, sharing the merged intensities (step 6) |

## The steps

1. **`mxi_import`, several masters** -- built: one experiment a master, each
   with its own models, as dials.import writes them (`docs/import.md`). With it, an
   offset without `offset_units` read in its transformation's units, and
   `mxeq compare-expt` comparing experiment by experiment.
2. **`mxi_find`** -- built: the spots of every experiment, from its own images,
   in one table with each spot's `id` and each identifier in the map. One sweep
   byte-identical to before; the insulin sweep split at image 150, found as two
   experiments, exactly the two halves found apart, every column.
3. **`mxi_index` and `mxi_refine`** -- built, and a bug found under them. The two
   halves indexed together gave 67.9 per cent where each alone gave 98.1: the
   second half's angles were 30 to 45 degrees where its images are 15 to 30.
   mxi took z, a frame's array index, from image one, where dxtbx takes it from
   the scan's first image -- `osc[0] + (z + 1 - image_range[0]) * width` -- so
   every scan not beginning at image one had its angles wrong by its start, in
   prediction, refinement, integration and the scan-varying models. Unseen while
   every scan began at image one; fixed, the frame offset the scan's first
   image's index, a scan's frames, its fraction and its start and end measured
   from there. Then the halves together index 98.1 per cent, as either alone;
   a scan from image one byte-identical through integration.
   **A crystal a sweep** (Graeme, on the four insulin sweeps: the cell held identical
   across them in refinement, and the RMSDs high). DIALS's protocol: index with
   one matrix, which puts the sweeps in one basis, then in refinement a crystal
   for each scan, refined apart. So now `mxi_index`'s macrocycles refine the one
   matrix together, and with several sweeps the last refines each sweep's
   crystal apart from it, indices then assigned through each reflection's own
   sweep's crystal; `mxi_refine` refines several sweeps a crystal each by
   default (it had `--separate`, never the default). `--shared-crystal` keeps
   one in either. One sweep byte-identical to before through both. Insulin's
   halves: refined apart, RMSD 0.297, 0.228 px, 0.244 images, as each half alone
   (0.288 to 0.302, 0.225 to 0.231, 0.236 to 0.246), where shared gave 0.296,
   0.242, 0.248; the two crystals' orientations within a thousandth of a degree
   -- one basis -- their cells within 0.03 per cent; symmetry and scaling as
   with one crystal. Small on one sweep cut in two, as it must be; four real
   sweeps are where it counts -- and on Graeme's four, it worked well.
4. **`mxi_integrate`** -- built. Each sweep integrated alone, by a whole run of
   the program on a list of that one sweep -- its reflections, by `id`, with
   their shoeboxes for its profile model -- `--postrefine` and all; then the
   tables joined, each row its sweep's `id`, and the lists joined, each sweep
   with its own profile model, and its `imageset_id` its image set's in the
   joined list -- each sweep's run writes 0, and left so every sweep's spots fell
   on the first image set in dials.image_viewer (Graeme, on four sweeps; the
   programs go by `id` and never noticed). A crystal the sweeps came in sharing, and none
   changed, is written once, shared, as DIALS writes it. One sweep never takes
   this path, and is byte-identical to before. `--save-shoeboxes` and
   `--save-profiles` are refused for several, not yet joined. On the insulin
   sweep split at image 150, integrated as two with the profile model given:
   every reflection well inside either half sums to exactly what the sweep
   integrated whole gives it -- summation depends on the box and the background,
   not the reference profiles, so a sweep's angles wrong would show there.
   Indexed, refined and integrated as two from their spots: 21029 reflections
   against the whole sweep's 21031, sigma_m 0.126 and 0.130 against 0.129.
5. **`mxi_symmetry`** -- built. Scaling's data gathering, which merging in P1
   uses, took only id 0, "one sweep, for now": now every sweep's observations,
   each with its sweep and its place in the rotation from its own sweep's scan
   -- and absorption's crystal frame from its own goniometer and beam -- and an
   id naming no experiment refused. Pooled, unscaled; each sweep must have a
   crystal and the cells agree to 2 per cent, the lattice being the first's.
   One sweep's symmetry and scaling byte-identical to before. On insulin as two
   sweeps: the same 20020 observations, I m -3 at NetZcc 8.41 against 8.42, the
   same space group, the crystal still shared. `mxi_scale` still refuses
   several: its model has one set of parameters, and two sweeps each running
   0 to 1 in rotation would share them -- step 6.
6. **`mxi_scale`** -- built. The model holds a block of parameters for each
   sweep -- scale, decay and absorption, from that sweep's own width of rotation
   -- one after another, an observation's sweep choosing its block; the fit sees
   only indices and gradients and needed no change. The one overall scale and
   the one overall B are fixed over every sweep's points together, in the
   normalisation and in the covariance's constraints, not each sweep's apart:
   that would erase the relative scale between sweeps. The restraints sum over
   every sweep, in the order they always were. A sweep with nothing to scale
   is refused. One sweep's scaling and report byte-identical to before.
   The whole sweep's intensities, split in two at scaling alone: I/sigma 19.0
   against 18.9, Rmerge 0.034 and CC1/2 0.987 either way -- the same, as it
   should be -- the two sweeps' scales meeting at the split (1.01, 1.00). The
   halves processed as two from their spots come out a little lower, I/sigma
   18.2 against 18.9: from integration, each 15 degree half's reference
   profiles drawn from fewer spots, as integrating sweeps apart does in DIALS
   too.

## Running it

As for one sweep, each program given the list of several; every reflection
carries its sweep's index as `id` from spot finding on:

```sh
mxi_import    a.nxs b.nxs                          # one experiment a master
mxi_find      imported.expt -o strong.refl         # each sweep from its own images
mxi_index     imported.expt strong.refl            # one matrix, then a crystal a sweep
mxi_refine    indexed.expt indexed.refl
mxi_integrate refined.expt refined.refl            # each sweep alone, then joined
mxi_symmetry  integrated.expt integrated.refl      # pooled
mxi_scale     symmetrized.expt symmetrized.refl    # a model a sweep
```

Integration wants the refined reflections, as for one sweep: its profile model
is estimated from the spots refinement used.

## How each is tested

Here: the insulin sweep, whose images are at hand, split in two -- imported
twice, `--image-range 1,150` and `151,300`, and the two lists joined into two
experiments of one crystal -- which runs every program's path for several sweeps
and has a single-sweep answer to agree with; and the tests built on it. Then the
real test, by Graeme: the four insulin sweeps of https://zenodo.org/records/8376818,
each at its own orientation, the whole chain from their masters.
