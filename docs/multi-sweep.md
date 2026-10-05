# More than one sweep

STATUS: a plan, 5 October 2026; steps 1 to 5 built. DIALS processes several sweeps of
one crystal together -- imported as one experiment each, indexed jointly to one
crystal, refined, integrated sweep by sweep, and symmetry and scaling over all
of them -- and mxi should, from `mxi_import` to `mxi_scale`. The data: Graeme's
threonine, two sweeps from Diamond's I03, with dials.import's experiment list to
reproduce; small-molecule data, whose intermediate files can be shared.

## Where each program stands

| program | more than one sweep |
| --- | --- |
| `mxi_import` | built: several masters, one experiment each (step 1) |
| `mxi_find` | built: each experiment from its own images, one table with `id`s (step 2) |
| `mxi_index` | joint indexing to one crystal: right, once a scan's frames were (step 3) |
| `mxi_refine` | several experiments, the crystal shared, built |
| `mxi_integrate` | built: each sweep alone, its own profile model, one table with `id`s (step 4) |
| `mxi_symmetry` | built: the sweeps pooled, every experiment reindexed alike (step 5) |
| `mxi_scale` | refuses more than one experiment |

## The steps

1. **`mxi_import`, several masters** -- built: one experiment a master, each
   with its own models, as dials.import writes them; on threonine every model
   of both as dials.import's but two fields (`docs/import.md`). With it, an
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
4. **`mxi_integrate`** -- built. Each sweep integrated alone, by a whole run of
   the program on a list of that one sweep -- its reflections, by `id`, with
   their shoeboxes for its profile model -- `--postrefine` and all; then the
   tables joined, each row its sweep's `id`, and the lists joined, each sweep
   with its own profile model. A crystal the sweeps came in sharing, and none
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
6. **`mxi_scale`**: a scale, decay and absorption model for each sweep, sharing
   the merged intensities -- the largest step.

## How each is tested

Here, without the threonine images: the insulin sweep, whose images are at
hand, split in two -- imported twice, `--image-range 1,150` and `151,300`, and the
two lists joined into two experiments of one crystal -- which runs every program's path for several
sweeps and has a single-sweep answer to agree with. Then on threonine, from
Graeme's intermediate files: `strong.refl` onwards, against DIALS's at each
boundary.
