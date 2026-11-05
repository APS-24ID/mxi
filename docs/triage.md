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

The twin tests use 10 A to a high-resolution limit chosen as xtriage's: the
finest d, shell by shell from 10 A, to which at least 85 per cent of the
reflections possible have I/sigma of 3 or more -- never coarser than 3.5 A,
never finer than the data (`--d-max`, `--d-min` to choose). The intensities
are normalised by resolution, E^2 = I / (epsilon <I/epsilon>) in shells of
about 400 reflections.

## Statistics that do not depend on twin laws

| | acentric untwinned | perfect twin | centric untwinned | perfect twin |
| --- | --- | --- | --- | --- |
| <I^2>/<I>^2 | 2.000 | 1.500 | 3.000 | 2.000 |
| <F>^2/<F^2> | 0.785 | 0.885 | 0.637 | 0.785 |
| <\|E^2-1\|> | 0.736 | 0.541 | 0.968 | 0.736 |

**N(z)**, the cumulative distribution of E^2 against Wilson's -- 1 - exp(-z)
acentric, erf(sqrt(z/2)) centric: twinning makes the acentric curve
sigmoidal.

**The L test** (Padilla & Yeates 2003): L = (E1^2 - E2^2) / (E1^2 + E2^2)
for each acentric reflection and a neighbour h + (2h', 2k', 2l'), 2 <=
|dh|+|dk|+|dl| <= 8. Untwinned, <|L|> 1/2 and <L^2> 1/3; perfectly twinned
3/8 and 1/5. Its multivariate Z score is the distance of the two from the
untwinned values, as xtriage's: over 3.5 is abnormal.

## Twin laws

The rotations the lattice's metric allows within 3 degrees' obliquity (Le
Page; `--max-delta`) beyond the space group's point group, one for each
coset, a two-fold where there is one -- merohedral where it keeps the cell's
metric exactly. For each: **R obs**, sum |I1 - I2| / sum (I1 + I2) over its
pairs of acentric reflections (about 0.5 for unrelated intensities, small for
equal ones); and **the H test's twin fraction** (Yeates 1997), from the
strongest half of the pairs, <|H|> = 1/2 - alpha.

## The summary

As xtriage's: the data are called twinned only where both kinds of
statistic say so.

- The L test abnormal as twinning makes it, and a twin law possible:
  twinning could be why; the largest H-test fraction is given.
- A twin law's pairs nearly equal (H-test fraction over 0.4) while the L test
  shows no twinning: a symmetry element missed, more likely than a twin.
- The L test abnormal with no twin law possible: something else.
- Neither: no twinning is suspected.

## Where it differs from xtriage

- **The L test's direction.** xtriage's Z score is two-sided. Here it is
  reported as xtriage reports it, but taken for twinning only where <|L|> is
  below 1/2, as twinning moves it. Above, the spread is wider than untwinned
  data's, which twinning never makes: the zones a missed symmetry makes
  centric, taken as acentric, do -- xtriage's own example of a symmetry
  missed has <|L|> 0.512 -- and translational NCS can.
- **The L test's pairs.** A neighbour the lattice's own symmetry relates to
  the reflection is not used: where a symmetry was missed, such pairs are
  equal, L is 0, and they had made the test look like twinning.
- **L from E^2, not I,** so that a neighbour's being a little further out
  does not count against it.
- **Normalisation in resolution shells only.** xtriage first corrects
  anisotropy, by a maximum-likelihood Wilson fit; that is to come.

## Tested, and to come

Tested on intensities drawn as Wilson's statistics have them, where every
figure has a known value: untwinned, twinned at 0.3 and perfectly by a
merohedral law, and P 4 2 2 data taken as P 4 (`tests/test_triage.cc`); and
end to end, on scaled files of such intensities, merged as the command merges
them (`python/tests/test_triage.py`). Not yet against xtriage itself on the
same real data.

To come: the Britton plot and xtriage's maximum-likelihood twin fraction;
translational NCS from the Patterson; anisotropy corrected first; Wilson
statistics and outliers; the data's health -- low-resolution completeness,
ice rings, anomalous measurability.
