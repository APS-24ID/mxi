# mxi triage: what the data say of twinning

`mxi triage` is mxi's phenix.xtriage (Zwart, Grosse-Kunstleve & Adams; CCP4
newsletters 42 and 43, 2005), begun: what the scaled intensities say of
twinning, run after scaling.

```sh
mxi triage scaled.rflx
mxi triage scaled.rflx --d-min 3      # the twin tests to 3 A
```

Every crystal's scaled observations are merged, as `mxi correlation` merges
them -- scaled, not outliers, partiality at least 0.4, cctbx's merging of
equivalents -- in the file's space group, Friedel mates together, the cell
regularised to the group. Its output is also in `mxi_triage.log`.

## The data used

As xtriage's, from its source (cctbx's mmtbx/scaling: data_statistics.py,
twin_analyses.py, absolute_scaling.py, twinning.h):

- **The range**: 10 A to a limit chosen from the data's strength -- the data's
  whole range in 15 shells of equal volume, the last left out; of those finer
  than 3.5 A, the finest in which more than 85 per cent of the reflections
  possible have I > 3 sigma; 3.5 A where none is (`--d-max`, `--d-min`).
- **Normalisation**: <I/eps> by a Gaussian kernel -- as wide as the d*^2 of
  the 50th reflection from low resolution, less the least -- at 23 Chebyshev
  nodes, a 13-term Chebyshev series fitted to its log; E^2 = I over that.
  As xtriage's, eps is not divided out of E^2: the few acentric reflections on
  a rotation axis come out eps times too strong, which raises the moments a
  little (perfectly twinned P 4 data: <I^2>/<I>^2 1.58, not 1.5; with eps
  divided, 1.5). Kept, to compare with xtriage.
- **The positive E^2 only**, for the moments, N(z) and the L test.

## Statistics that do not depend on twin laws

| | acentric untwinned | perfect twin | centric untwinned | perfect twin |
| --- | --- | --- | --- | --- |
| <I^2>/<I>^2 | 2.000 | 1.500 | 3.000 | 2.000 |
| <F>^2/<F^2> | 0.785 | 0.885 | 0.637 | 0.785 |
| <\|E^2-1\|> | 0.736 | 0.541 | 0.968 | 0.736 |

**N(z)**, the cumulative distribution of E^2 against Wilson's -- 1 - exp(-z)
acentric, erf(sqrt(z/2)) centric: twinning makes the acentric curve
sigmoidal.

**The L test** (Padilla & Yeates 2003), as xtriage's: L = (E1^2 - E2^2) /
(E1^2 + E2^2) for each positive acentric reflection and one neighbour drawn
for it, h + (dh, dk, dl), each even, |dh|+|dk|+|dl| 2, 4 or 6 -- none where
that neighbour was not measured. Untwinned, <|L|> 1/2 and <L^2> 1/3;
perfectly twinned 3/8 and 1/5. Its **multivariate Z score** is xtriage's:
the distance of (<|L|>, <L^2>) from the mean of some 1200 PDB data sets,
(0.4851, 0.3164), by their covariance -- how far from typical good data, not
a test against the ideal values, so not growing with the number of pairs.

## Twin laws

The rotations the lattice's metric allows within 3 degrees' obliquity (Le
Page; `--max-delta`) beyond the space group's point group, one for each
coset, a two-fold where there is one -- merohedral where it keeps the cell's
metric exactly. For each, on the intensities as merged: **R obs**, sum |I1 -
I2| / sum |I1 + I2| over every reflection with a mate, centric too, each from
both sides (about 0.5 for unrelated intensities, small for equal ones); and
**the H test's twin fraction** (Yeates 1997), from the strongest half of the
acentric pairs by I1 + I2, <|H|> = 1/2 - alpha.

## The summary

As xtriage's (twin_results_interpretation): the intensity statistics are
abnormal where Z is 3.5 or more **and** <|L|> is below 1/2, as twinning moves
it; the data are called twinned only where they are and a twin law is
possible.

- Abnormal, a twin law possible: twinning could be why; the largest H-test
  fraction is given.
- A twin law's pairs nearly equal (H-test fraction over 0.4) while the
  statistics are not abnormal: a symmetry element missed, more likely than a
  twin.
- Abnormal with no twin law possible: something else.
- Z 3.5 or more with <|L|> 1/2 or more: a spread wider than good data's, which
  twinning never makes -- translational NCS can.
- Otherwise: no twinning is suspected.

## Where it still differs from xtriage

- **The basic filters before the twin analysis**: xtriage first corrects
  anisotropy (a maximum-likelihood anisotropic Wilson fit) and removes Wilson
  outliers; not yet here.
- **The H test's fraction**: from <|H|>; xtriage reports the fraction from a
  fit to the cumulative distribution of H, the two agreeing to 0.001 on the
  insulin data compared.
- **The L test's draw of neighbours**: a different random generator, so the
  pairs differ; the statistics agree to the third figure.
- **A twin law's name**: the same law, perhaps another member of its coset --
  -h,-l,-k here where xtriage writes -l,-k,-h for I 2 3, the two related by
  the three-fold.

## Tested, and to come

Tested on intensities drawn as Wilson's statistics have them, where every
figure has a known value: untwinned, twinned at 0.3 and perfectly by a
merohedral law, and P 4 2 2 data taken as P 4 (`tests/test_triage.cc`); and
end to end, on scaled files of such intensities, merged as the command merges
them (`python/tests/test_triage.py`). Compared with xtriage on insulin, I 2 3,
DIALS's processing against mxi's: the moments, N(z) and the L test within
0.03, the same twin law, both untwinned. A comparison on the same data -- mxi's
scaled.mtz given to xtriage -- is still to make.

To come: the Britton plot and xtriage's maximum-likelihood twin fraction;
translational NCS from the Patterson; anisotropy corrected first; Wilson
statistics and outliers; the data's health -- low-resolution completeness,
ice rings, anomalous measurability.
