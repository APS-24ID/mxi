# Processing 🐮🐷🧑 with mxi

DRAFT: the commands and the story are final; the program output marked
*to come* is to be pasted in from a run of
[`cows-pigs-people.sh`](cows-pigs-people.sh) on the tutorial data.

## Introduction

mxi processes data from many crystals with the same steps as one: read the
data, find spots, index, refine, integrate, determine the symmetry, scale.
This tutorial follows Graeme Winter's DIALS tutorial of the same data
([COWS_PIGS_PEOPLE](https://github.com/graeme-winter/dials_tutorials/blob/main/iucr-2026/COWS_PIGS_PEOPLE.md)),
step for step, so the two can be read side by side. It *starts* with data
from many crystals -- first of one kind, then a mixture -- and shows how to
find which data sets belong together when the differences between them are
small: here, the same protein from three species.

## The data

Thirty-six small rotation data sets of insulin, recorded on i24 at Diamond
Light Source: twelve crystals each of bovine (`CIX`), porcine (`PIX`) and
human (`X`) insulin. All are I2₁3 with very similar cells, so they *can* be
merged together -- it will work, and give poor results, a mixture of three
structures.

mxi reads NXmx only, so the tutorial uses the data transcribed from the
original CBF to NXmx: [Tutorial data for MXI](https://zenodo.org/records/23423865),
one `.nxs` file a data set. Download the 36 into `data`:

```sh
mkdir data && cd data
for set in CIX1 CIX2 CIX3 CIX5 CIX6 CIX8 CIX9 CIX10 CIX11 CIX12 CIX14 CIX15 \
           PIX5 PIX6 PIX7 PIX8 PIX9 PIX10 PIX11 PIX12 PIX13 PIX14 PIX15 PIX16 \
           X1 X2 X3 X4 X5 X6 X7 X8 X9 X11 X13 X14; do
  wget https://zenodo.org/records/23423865/files/${set}.nxs
done
cd ..
```

## The workflow

| step | mxi | DIALS |
| --- | --- | --- |
| read the data | `mxi_import` | `dials.import` |
| find spots | `mxi_find` | `dials.find_spots` |
| index | `mxi_index --no-joint` | `dials.index joint=false` |
| refine | `mxi_refine` | `dials.refine` |
| integrate | `mxi_integrate` | `dials.integrate` |
| symmetry, and the crystals into one setting | `mxi_symmetry` | `dials.cosym` |
| scale | `mxi_scale` | `dials.scale` |
| which data sets belong together | `mxi_correlation` | `dials.correlation_matrix` |
| export | `mxi_export` | `dials.export` |

Two things differ from one crystal. The crystals share no orientation, so
each is indexed on its own (`--no-joint`). And the lattice's symmetry (cubic,
I) is higher than the data's (I2₁3, Laue group m-3), so each crystal can be
indexed in either of two settings that the data tell apart: `mxi_symmetry`
brings them to one, with cosym, whenever it is given several crystals --
there is no separate program for it.

Each program reads the last one's `.rflx` and writes its own, and writes what
it prints to `mxi_<program>.log`.

## The cows

First the twelve bovine data sets alone, in a directory of their own:

```sh
mkdir cows && cd cows
mxi_import ../data/CIX*.nxs
```

*to come: mxi_import's summary -- twelve sweeps, 100 images each.*

Check this before anything else: as with DIALS, many problems are solved
here. `mxi_show imported.rflx` shows mxi's understanding of every sweep --
its images, its rotation, the detector and the beam.

### Spot finding

```sh
mxi_find imported.rflx
```

*to come: the spots found on each sweep.*

`mxi_find --gpu` finds them on a GPU, where the build has one: the same
spots.

### Indexing and refinement

```sh
mxi_index strong.rflx --no-joint
mxi_show indexed.rflx | grep experiment
```

*to come: each crystal's cell, all close to 67 A and 109.5 deg -- the
primitive cell of cubic insulin.*

DIALS can hold every crystal to an expected cell (`unit_cell=`) when one is an
outlier; mxi cannot yet (`docs/outstanding.md`), and it is not needed here.
Then refinement:

```sh
mxi_refine indexed.rflx
```

*to come: each sweep's RMSDs, well under a pixel and under a fifth of an
image.*

Each sweep is refined on its own, its beam and detector its own: these are
different crystals, mounted apart.

### Integration

```sh
mxi_integrate refined.rflx
```

*to come: each sweep's profile model and the reflections integrated.*

`mxi_integrate --gpu` fits the profiles on a GPU.

### Symmetry and scaling

```sh
mxi_symmetry integrated.rflx
```

*to come: cosym's one setting for the twelve, the operation each needed, and
the space group, I 2 3 -- indistinguishable by its absences from I 21 3.*

Then scaling, the resolution limit estimated from CC½ and the data scaled
again to it:

```sh
mxi_scale symmetrized.rflx --d-min-auto
```

*to come: the summary table.*

### Do they belong together?

So far we have assumed the twelve are isomorphous -- that they merge well --
without testing it. `mxi_correlation` measures how alike every pair of data
sets is and clusters them:

```sh
mxi_correlation scaled.rflx --html correlation.html
```

*to come: one cluster of twelve, no outliers.*

One cluster: scaling them together was right. Look at `correlation.html` --
the correlation matrix, its dendrogram, and the clusters in the coordinates
cosym found.

## Cows, pigs and people

Now all thirty-six, from three species: isomorphous on a coarse scale,
different at the scale of single residues -- different enough to spoil a
merge, close enough to merge by accident. First, ignorantly, as one:

```sh
cd .. && mkdir all && cd all
mxi_import ../data/CIX*.nxs ../data/PIX*.nxs ../data/X*.nxs
mxi_find imported.rflx
mxi_index strong.rflx --no-joint
mxi_refine indexed.rflx
mxi_integrate refined.rflx
mxi_symmetry integrated.rflx
mxi_scale symmetrized.rflx --d-min-auto
```

*to come: the 36 in one setting; the merging statistics, poorer than the
cows' alone.*

The cells alone do not separate the three -- they are too alike -- but the
intensities do:

```sh
mxi_correlation scaled.rflx --html correlation.html --write-clusters
```

*to come: three clusters of twelve.*

`correlation.html` shows the three blocks plainly in the correlation matrix.
`--write-clusters` writes each cluster as its own file, `cluster_0.rflx`,
`cluster_1.rflx` and `cluster_2.rflx`, ready to scale alone:

```sh
mkdir 0 1 2
cd 0 && mxi_scale ../cluster_0.rflx --d-min-auto && cd ..
cd 1 && mxi_scale ../cluster_1.rflx --d-min-auto && cd ..
cd 2 && mxi_scale ../cluster_2.rflx --d-min-auto && cd ..
```

*to come: each cluster's merging statistics, far better than the three
together.*

## Export

```sh
mxi_export scaled.rflx
```

writes `scaled.mtz`, unmerged, for CCP4 or PHENIX.

## Working alongside DIALS

Any step can be DIALS's. `mxi_convert scaled.rflx` writes the DIALS pair, so
`dials.image_viewer`, `dials.reciprocal_lattice_viewer` and `dials.report`
show mxi's results.

## References

1. Gildea & Winter (2018), cosym. *Acta Cryst.* D74, 405–410.
2. Thompson et al. (2025), intensity-based clustering. *Acta Cryst.* D81,
   278–290.
3. Winter et al. (2018), DIALS. *Acta Cryst.* D74, 85–97.
