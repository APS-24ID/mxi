# Processing 🐮🐷🧑 with mxi

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
| read the data | `mxi import` | `dials.import` |
| find spots | `mxi find` | `dials.find_spots` |
| index | `mxi index --no-joint` | `dials.index joint=false` |
| refine | `mxi refine` | `dials.refine` |
| integrate | `mxi integrate` | `dials.integrate` |
| symmetry, and the crystals into one setting | `mxi symmetry` | `dials.cosym` |
| scale | `mxi scale` | `dials.scale` |
| which data sets belong together | `mxi correlation` | `dials.correlation_matrix` |
| export | `mxi export` | `dials.export` |

Two things differ from one crystal. The crystals share no orientation, so
each is indexed on its own (`--no-joint`). And the lattice's symmetry (cubic,
I) is higher than the data's (I2₁3, Laue group m-3), so each crystal can be
indexed in either of two settings that the data tell apart: `mxi symmetry`
brings them to one, with cosym, whenever it is given several crystals --
there is no separate program for it.

Each program reads the last one's `.rflx` and writes its own, and writes what
it prints to `mxi_<program>.log`.

## The cows

First the twelve bovine data sets alone, in a directory of their own:

```sh
mkdir cows && cd cows
mxi import ../data/CIX*.nxs
```

Each sweep, as mxi reads it:

```
  wavelength 0.999870 A
  panel /entry/instrument/detector/module: 2463 x 2527 pixels of 0.1720 x 0.1720 mm, Si 0.450 mm thick, mu 4.2064/mm
    distance 200.000 mm; the beam meets it at pixel 1292.85, 1311.74; trusted to 91798
  goniometer: omega (scan)
  scan: images 1 to 100, from 0.0000 deg, 0.1000 deg each, 10.0 deg in all
```

Twelve sweeps of 100 images, 10 degrees each: small data sets, as from many
crystals.

Check this before anything else: as with DIALS, many problems are solved
here. `mxi show imported.rflx` shows mxi's understanding of every sweep --
its images, its rotation, the detector and the beam.

### Spot finding

```sh
mxi find imported.rflx
```

```
Experiment 1 of 12: images .../CIX1.nxs
Found 216272 signal pixels on 100 frames
Extracted 7626 spots
Experiment 2 of 12: images .../CIX10.nxs
Found 72711 signal pixels on 100 frames
Extracted 6099 spots
...
Wrote 83544 reflections to strong.rflx (11.6 MB, most of it shoeboxes)
```

The second, CIX10, is the weakest of the twelve -- a third of the first's
signal pixels -- which matters later.

`mxi find --gpu` finds them on a GPU, where the build has one: the same
spots.

### Indexing and refinement

```sh
mxi index strong.rflx --no-joint
mxi show indexed.rflx | grep experiment
```

```
12 sweeps indexed apart: 12 lattices, 0 sweeps not indexed
    experiment 1: 67.406 67.457 67.467 A, 109.465 109.438 109.505 deg; P 1
    experiment 3: 67.317 67.429 67.408 A, 109.560 109.420 109.413 deg; P 1
    experiment 5: 67.425 67.442 67.404 A, 109.465 109.445 109.505 deg; P 1
    ...
```

Every crystal's cell close to 67.4 A and 109.5 deg: the primitive cell of
body-centred cubic insulin, 78 A. (The even-numbered experiments are each
sweep's own, holding its unindexed spots.)

If one crystal's cell were an outlier, the known cell would choose every
crystal's basis, as `unit_cell=` does for dials.index -- given conventional,
with its lattice:

```sh
mxi index strong.rflx --no-joint --known-cell 78,78,78,90,90,90 --known-space-group I213
```

It is not needed here. Then refinement:

```sh
mxi refine indexed.rflx
```

```
Refining 24 experiments sweep by sweep, 12 runs, each sweep's beam and detector its own
RMSD over those 5379: 0.2053 px, 0.2363 px, 0.1599 images
RMSD over those 4188: 0.2078 px, 0.2362 px, 0.1462 images
RMSD over those 5428: 0.2104 px, 0.2252 px, 0.1925 images
...
```

About a fifth of a pixel and of an image for every sweep -- one, 0.4 of an
image in rotation.

Each sweep is refined on its own, its beam and detector its own: these are
different crystals, mounted apart.

### Integration

```sh
mxi integrate refined.rflx
```

```
Profile model: sigma_b 0.0467 deg, sigma_m 0.1262 deg (estimated from the indexed spots); foreground to 3.0 sigma
Integrated 31295 reflections
Profile model: sigma_b 0.0364 deg, sigma_m 0.0624 deg (estimated from the indexed spots); foreground to 3.0 sigma
Integrated 31258 reflections
...
Wrote 375723 reflections of 24 experiments, and the models, to integrated.rflx
```

Each sweep its own profile model: the crystals' reflecting ranges differ.

`mxi integrate --gpu` fits the profiles on a GPU.

### Symmetry and scaling

```sh
mxi symmetry integrated.rflx
```

```
One setting for 12 crystals (cosym): the lattice's 24 rotations; principal components 41.1% 18.1% 1.5% 1.5%, 2 kept
Recommended space group: I 2 3
  indistinguishable by its absences from I 21 3
Reindex operator: b+c,a+c,a+b
Cell 77.962 77.962 77.962 A, 90.000 90.000 90.000 deg
```

The twelve in one setting, and the data in the conventional cubic cell.
I 2 3 and I 2_1 3 cannot be told apart by their absences: insulin is I 2_1 3,
which the structure, not the data, decides.

Then scaling, the resolution limit estimated from CC½ and the data scaled
again to it:

```sh
mxi scale symmetrized.rflx --d-min-auto
```

```
--d-min-auto: CC half falls to 0.30 at 1.28 A with everything scaled; scaling again, as --d-min 1.28 would
                                              Overall     Low    High
High resolution limit                            1.28    3.47    1.28
Completeness                                     99.9   100.0    97.6
Multiplicity                                     13.1    12.8    10.1
I/sigma                                          12.3    32.5     1.4
Rmerge(I)                                       0.091   0.062   2.183
Rpim(I)                                         0.026   0.018   0.696
CC half                                         0.996   0.994   0.244
```

### Do they belong together?

So far we have assumed the twelve are isomorphous -- that they merge well --
without testing it. `mxi correlation` measures how alike every pair of data
sets is and clusters them:

```sh
mxi correlation scaled.rflx --html correlation.html
```

```
Normalised by a Wilson model fitted to each data set (the
isotropic equivalent of its anisotropic B): 17.47 19.81 18.63 19.73 16.25 19.30 19.22 17.63 18.53 18.75 16.86 19.21 A^2
correlations: 0.924 to 0.972, median 0.954
3 dimensions, at the elbow; the functional there 0.000221176
OPTICS (Ankerst et al. 1999), xi 0.05, max distance 0.5: no min_samples gave two clusters; one cluster, as DIALS's single-cluster heuristic finds it, min_samples 5
  cluster 0, 12 data sets: 1 3 5 7 9 11 13 15 17 19 21 23
```

Each data set is first normalised by its own Wilson fall-off: correlated
as scaled, a crystal whose intensities fall off differently -- CIX10 here --
looks less like the rest than it is.

One cluster: scaling them together was right. Look at `correlation.html` --
the correlation matrix, its dendrogram, and the clusters in the coordinates
cosym found.

## Cows, pigs and people

Now all thirty-six, from three species: isomorphous on a coarse scale,
different at the scale of single residues -- different enough to spoil a
merge, close enough to merge by accident. First, ignorantly, as one:

```sh
cd .. && mkdir all && cd all
mxi import ../data/CIX*.nxs ../data/PIX*.nxs ../data/X*.nxs
mxi find imported.rflx
mxi index strong.rflx --no-joint
mxi refine indexed.rflx
mxi integrate refined.rflx
mxi symmetry integrated.rflx
mxi scale symmetrized.rflx --d-min-auto
```

```
One setting for 36 crystals (cosym): the lattice's 24 rotations; principal components 34.8% 14.6% 1.2% 1.2%, 2 kept
Recommended space group: I 2 3
...
High resolution limit                            1.26    3.41    1.26
Multiplicity                                     38.5    36.8    23.8
I/sigma                                          11.7    23.8     2.6
Rmerge(I)                                       0.204   0.164   2.546
CC half                                         0.998   0.991   0.367
```

The merge works, and is poor: the inner shell's R_merge 0.164, against
0.062 for the cows alone, and the error model's b 0.085 where the cows'
was 0.061 -- three structures averaged.

The cells alone do not separate the three -- they are too alike -- but the
intensities do:

```sh
mxi correlation scaled.rflx --html correlation.html --write-clusters
```

```
correlations: 0.571 to 0.972, median 0.807
3 dimensions, at the elbow; the functional there 0.0523707
OPTICS (Ankerst et al. 1999), xi 0.05, max distance 0.5: min_samples 6 (the initial guess 6), silhouette 0.817, score 0.817
  cluster 0, 12 data sets: 1 3 5 7 9 11 13 15 17 19 21 23
  cluster 1, 12 data sets: 25 27 29 31 33 35 37 39 41 43 45 47
  cluster 2, 12 data sets: 49 51 53 55 57 59 61 63 65 67 69 71
```

Three clusters of twelve -- in the order the data sets were imported, the
cows, the pigs and the people -- with no outliers.

`correlation.html` shows the three blocks plainly in the correlation matrix.
`--write-clusters` writes each cluster as its own file, `cluster_0.rflx`,
`cluster_1.rflx` and `cluster_2.rflx`, ready to scale alone:

```sh
mkdir 0 1 2
cd 0 && mxi scale ../cluster_0.rflx --d-min-auto && cd ..
cd 1 && mxi scale ../cluster_1.rflx --d-min-auto && cd ..
cd 2 && mxi scale ../cluster_2.rflx --d-min-auto && cd ..
```

| cluster | d_min (A) | I/sigma, inner | R_merge, inner | CC half |
| --- | --- | --- | --- | --- |
| 0, cows | 1.28 | 32.5 | 0.062 | 0.996 |
| 1, pigs | 1.32 | 22.0 | 0.095 | 0.994 |
| 2, people | 1.17 | 36.2 | 0.051 | 0.998 |
| all 36 together | 1.26 | 23.8 | 0.164 | 0.998 |

Each species alone merges far better than the three together; the cows'
cluster is exactly the cows processed alone.

## Export

```sh
mxi export scaled.rflx
```

writes `scaled.mtz`, unmerged, for CCP4 or PHENIX.

## Working alongside DIALS

Any step can be DIALS's. `mxi convert scaled.rflx` writes the DIALS pair, so
`dials.image_viewer`, `dials.reciprocal_lattice_viewer` and `dials.report`
show mxi's results.

## References

1. Gildea & Winter (2018), cosym. *Acta Cryst.* D74, 405–410.
2. Thompson et al. (2025), intensity-based clustering. *Acta Cryst.* D81,
   278–290.
3. Winter et al. (2018), DIALS. *Acta Cryst.* D74, 85–97.
