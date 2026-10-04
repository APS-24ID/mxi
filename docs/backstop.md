# The backstop: its shadow, its flare, and the outliers they make

STATUS: a plan, 4 October 2026, not yet built. Items 52 and 53 in
`docs/outstanding.md` are the findings; this is what is to be done about them,
in what order, and how it will be judged.

## Terms

* **The backstop shadow**: the region the backstop hides, beside the beam
  centre. Its pixels are valid -- not module gaps, not marked bad -- and see
  almost no X-rays: low counts.
* **The flare**: the bright region around the backstop's edge, high background
  with a steep gradient: X-rays scattered there, real, but not from the crystal.
* A **module gap** is masked already (its pixels carry the tile-join marker) and
  is not part of this.

## What was found

On ferritin (zenodo 8376818), DIALS's inner shell had the better CC1/2 (0.9998
against mxi's 0.9992), and the two data sets' merged intensities correlated
there at only 0.971 against 0.99 and above everywhere else (`mxeq unique`). The
reflections responsible are the lowest-resolution ones, beside the beam
(`mxeq observations`). Of (3,1,1)'s 48 observations, 21 fall in the horizontal
module gap, read 0 with no valid pixels, and are left out, rightly. The other
27, by scaled intensity:

| scaled intensity | observations | backgrounds | in scaling |
| --- | --- | --- | --- |
| 934 to 1078 | 13 | 22 to 27 | every one an outlier |
| 146 to 588 | 7 | 3.8 to 15 | every one an outlier |
| 23 to 99 | 5 | 0.5 to 2.8 | every one an outlier |
| 10.5 and 10.7 | 2 | 0.50 and 0.53 | **kept**: their mean, 10.6, merged |

**The intensity follows the background.** Scaled intensity over background is
some 40 for every one of them, the 13 unharmed and the rest alike -- 1000/25,
588/14.9, 274/7.6, 146/3.8. The backstop shadow's edge is not sharp: it
attenuates, the reflection's signal and the background under it by the same
factor. So an observation's background, against the background expected at its
distance from the beam, is its transmission -- measured, not guessed.

Two faults, then:

1. **Integration measures the shadow as if it were sample**, wholly or in part.
   A pixel in the shadow is valid, so a reflection there is integrated low, or
   as nothing, with a sigma from its own few counts.
2. **Scaling's outlier rejection can be outvoted.** Each observation is tested
   against the mean of the others in its group, weighted by 1/sigma^2, sigma its
   own; the worst beyond 6 is removed and the rest retested. An attenuated
   observation's sigma is small, the deepest smallest, so the weighted mean
   sinks towards them, the unharmed look deviant and go one at a time, and what
   is left is the deepest shadow.

No consensus can rescue (3,1,1): 14 of its 27 are attenuated, and their median,
588, is among them. The first fault has to be found where it happens, from the
background; robust rejection is the guard for the rarer case, a few bad among
many good.

Insulin, here, shows the first: its (1,1,0) is measured once, in the shadow, and
scaled at -0.04. `--d-max` in both programs leaves such reflections out by hand;
on ferritin it is being tested now, and Graeme will judge by refining a
structure against each data set.

## 1. The backstop shadow, found from the background

**A map of the expected background.** Before the one pass, read a sample of
frames spread over the scan -- every k-th, at least 50; on ferritin 1 in 20 of
3600, some 5 per cent more reading -- and sum each pixel's counts over them,
markers excluded. Group pixels in rings about the beam centre, a pixel or two
wide, the centre from the experiment's geometry; each ring's median sum, over its
valid pixels, is the background expected at that distance, unshadowed. The
ring is mostly not shadowed, so its median is not.

**Each reflection's transmission.** Integration measures each reflection's
background. Its transmission is that background over the map's expectation for
its pixels, on the same per-frame scale. A reflection whose transmission is
below t -- 0.9 to start, judged on the data -- and below by more than the
counting noise of the expectation, is attenuated: flagged, not profile fitted,
not scaled. This is the test for an abnormal background: it catches the deep
shadow and the partly attenuated edge alike, and on (3,1,1) it would leave the
13 and drop the 14.

**Pixels in deep shadow** -- a sum below a fifth of its ring's, by more than five
Poisson standard deviations of the ring's -- are also masked as a module gap is,
so that a reflection partly in deep shadow is fitted from its illuminated
pixels. A dead pixel the detector did not mark is caught the same way, rightly.

**The flare** raises a ring's pixels above its median, not below: it does not
look like attenuation, and is left to part 3.

**Assumed, and to be checked.** The backstop is fixed to the detector, so the
shadow does not move through the scan: one map from the sum. A shadow that
moves -- the goniometer's, the cryostream's -- would need a map per scan block;
not in the first version, which reports how much it flagged and where, so a
moving shadow shows as too little.

**Interface.** `mxi_integrate --backstop`, off at first; on by default once it
has been judged. It reports the reflections flagged and the pixels masked, and
can write the transmission map (`--write-transmission`) to view beside a frame.
The flag goes into the table, so `mxeq` can count it.

## 2. Outlier rejection that the few cannot win

The same test -- an observation against the others in its group, the worst
beyond zmax removed and the rest retested -- with two changes, so that a few
observations with small, wrong sigmas cannot outvote many consistent ones:

* **The reference is robust**: the median of the others' scaled intensities,
  not their 1/sigma^2-weighted mean, so that two at 10 cannot drag twenty at
  1000.
* **The sigma is the one the observation would have if it were right**: its own
  variance plus its counting part for the difference between the reference and
  what it read, in its own counts-to-intensity units, through the error model
  as now. A deeply shadowed observation of a strong reflection is then many
  sigma low.

The weighted mean stays what scaling and merging use; only the rejection's
reference and denominator change. Weak reflections, where the median and the
weighted mean agree within their errors, should lose about what they lose now:
the number rejected, overall and by shell, will be reported against the present
test's.

## 3. The flare: measured before anything is done

Reflections in the flare have high backgrounds with steep gradients, and the
background model -- a constant over each box -- is the wrong shape there. But
(3,1,1)'s observations in it, backgrounds 22 to 27, agree among themselves. So,
before excluding them or fitting a sloped background: measure each reflection's
background gradient across its box, binned in `mxeq equivalents` against its
equivalents. Only if those in the flare are biased does something change -- a
plane for the background there, or a flag and leaving them out.

## Order, and how each is judged

1. **The backstop shadow** first, the fault where it happens. Judged by
   (3,1,1)'s 14 attenuated observations flagged and its 13 kept, insulin's
   (1,1,0) leaving the scaled data, the transmission map drawn beside a frame
   matching the shadow and its edge, nothing flagged away from the backstop,
   `--d-max`'s results matched without it, `mxeq unique`'s inner shell against
   DIALS, and Graeme's refinement.
2. **Outlier rejection** second, the guard: judged by a planted test (a group of
   20 consistent and 3 near zero with small sigmas: the 3 go, the 20 stay), by
   the rejected counts against the present test's, and by refinement.
3. **The flare** third, measured, and changed only if the measurement says so.

Each default stays as it is until the change has been judged, and each is
byte-identical when off.
