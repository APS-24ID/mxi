# More than one sweep

STATUS: a plan, 5 October 2026; step 1 built. DIALS processes several sweeps of
one crystal together -- imported as one experiment each, indexed jointly to one
crystal, refined, integrated sweep by sweep, and symmetry and scaling over all
of them -- and mxi should, from `mxi_import` to `mxi_scale`. The data: Graeme's
threonine, two sweeps from Diamond's I03, with dials.import's experiment list to
reproduce; small-molecule data, whose intermediate files can be shared.

## Where each program stands

| program | more than one sweep |
| --- | --- |
| `mxi_import` | built: several masters, one experiment each (step 1) |
| `mxi_find` | refuses more than one experiment |
| `mxi_index` | joint indexing to one crystal built, tried once on two sweeps of insulin |
| `mxi_refine` | several experiments, the crystal shared, built |
| `mxi_integrate` | uses only the first experiment -- silently |
| `mxi_symmetry`, `mxi_scale` | refuse more than one experiment |

## The steps

1. **`mxi_import`, several masters** -- built: one experiment a master, each
   with its own models, as dials.import writes them; on threonine every model
   of both as dials.import's but two fields (`docs/import.md`). With it, an
   offset without `offset_units` read in its transformation's units, and
   `mxeq compare-expt` comparing experiment by experiment.
2. **`mxi_find`**: the spots of every experiment, from its own images, in one
   table with each spot's `id` -- the experiment's index, as dials.find_spots
   writes it.
3. **`mxi_index` and `mxi_refine`**: built; to be checked with the `id`s from
   step 2, and on the two sweeps.
4. **`mxi_integrate`**: each sweep from its own images, with its own reference
   profiles and profile model, into one table with `id`s. Silently using only
   the first experiment ends here whatever else does.
5. **`mxi_symmetry`**: the sweeps' reflections pooled for the Laue group, and
   every experiment reindexed alike.
6. **`mxi_scale`**: a scale, decay and absorption model for each sweep, sharing
   the merged intensities -- the largest step.

## How each is tested

Here, without the threonine images: the insulin sweep, whose images are at
hand, split in two -- imported twice, `--image-range 1,150` and `151,300`, and the
two lists joined into two experiments of one crystal -- which runs every program's path for several
sweeps and has a single-sweep answer to agree with. Then on threonine, from
Graeme's intermediate files: `strong.refl` onwards, against DIALS's at each
boundary.
