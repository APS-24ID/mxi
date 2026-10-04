# The backstop: its shadow, its flare, and the outliers they make

STATUS: a plan, 4 October 2026; revised the same day to Graeme's metric, the
background's own dispersion, in place of a model of the images. Step 1, the
measurement, is built, and on ferritin it says the dispersion cannot find the
attenuation: see "Measured so far" and "What can, then". Items 52 and 53 in
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

## 1. The shadow's edge, from each reflection's own background

**Not a model of the images.** An earlier draft built a map of the background
expected at each distance from the beam, from a sample of frames, and compared
each reflection with it. Graeme's objection stands: anything that must look at
images to model the background is a heuristic, and will fail on the data it
did not anticipate. Instead, a test on what each reflection's shoebox already
holds.

**The dispersion of its background.** On a flat background a photon-counting
detector's pixels are Poisson: their variance over their mean, the dispersion,
is about 1. A reflection on the shadow's edge sits on a ramp from nearly nothing
to the full background, and its dispersion is far above 1: across a box whose
background ramps by D with mean m it is about 1 + D^2 / 12m, so on (3,1,1)'s
edge -- a ramp of some 25 under a mean of some 5 -- about 11, where a flat
background gives 1 at any level.

**A test, not a threshold.** For n Poisson pixels with a flat mean, (n - 1)
times the dispersion is chi-squared on n - 1 degrees of freedom. A reflection is
flagged when its background is incompatible with that at a stated false-alarm
probability -- 1e-6, a few false flags among millions -- and is then not profile
fitted, and not scaled. A detector that is not photon counting has its gain,
`--gain`, divided out.

**What it does not catch.** A reflection wholly in the shadow sits on a flat,
nearly empty background: dispersion about 1, intensity about 0. That one is
left to scaling (part 2), which must reject it -- once the edge is flagged,
(3,1,1) has 13 unharmed observations against 2 wholly shadowed, and a test the
many cannot lose to the few rejects the 2.

**What may raise it besides.** The flare has a gradient too, though under a
mean of 25 a ramp must be large to double the dispersion; and on dense data a
neighbour's spot in a box's background region adds variance. Which background
pixels the dispersion is taken over -- all the valid ones, or those the robust
background fit kept -- decides both, and is to be measured, not chosen in
advance.

**First, measurement only.** Each reflection's background dispersion written to
`integrated.refl` as `background.dispersion`, with how many pixels it was taken
over, so that on ferritin `mxeq observations` and `mxeq equivalents` can show
whether the edge stands out, where the flare falls, and whether dense regions
flag falsely. Then the test, `mxi_integrate --flag-background`, off until it has
been judged; the flag goes into the table, so `mxeq` can count it.

## 2. Outlier rejection that the few cannot win

Essential, not only a guard: a reflection wholly in the backstop shadow passes
part 1's test, and must be rejected here.

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

## Measured so far

`mxi_integrate --save-background-parameters` writes `background.dispersion`
-- only when asked, so that by default the table is what it was -- over every
valid background pixel, and `background.dispersion_trimmed`, leaving out those beyond five
Poisson standard deviations of the robust mean, with
`num_pixels.background_trimmed`. `mxeq observations` lists both for each
observation, and `mxeq equivalents` bins the bias by the trimmed one.

On the 300 image insulin sweep, which has no flare:

| | median | 1st to 99th percentile | largest |
| --- | --- | --- | --- |
| over every pixel | 1.001 | 0.90 to 1.20 | 184 |
| trimmed | 0.996 | 0.90 to 1.07 | 1.3 |

Poisson, as a photon counter should be. The largest untrimmed values, 5 to 21,
are strong reflections beside the beam whose boxes hold a few bright pixels --
a neighbour's spot, their own tails -- and trimmed they are about 1: trimming
does what it is for. And insulin's shadowed (1,1,0), background 0.43, has a
dispersion of 1.09: wholly in the shadow, on a flat empty background, as part
1 says it would be, and left to part 2. Ferritin's flare and shadow edge are
the measurement still to make.

### On ferritin: the dispersion does not find the attenuation

(3,1,1)'s observations, trimmed dispersion against scaled intensity:

| observations | scaled intensity | background | trimmed dispersion |
| --- | --- | --- | --- |
| 13 unharmed | 900 to 1032 | 22 to 27 | 3.5 to 5.4 |
| 7 partly attenuated | 142 to 572 | 3.8 to 15 | 1.9 to 2.5 |
| 5 deeply attenuated | 22 to 96 | 0.5 to 2.8 | 1.3 to 2.1 |
| 2 wholly shadowed, kept | 10 | 0.5 | 1.35 and 1.43 |

The opposite of the hypothesis. The unharmed lie in the flare, whose gradient
under a background of 25 gives them the highest dispersion; and the attenuation
is uniform across each box -- the shadow's edge is wider than a box -- so an
attenuated box is a dimmed copy of an unharmed one, signal, background and the
flare's gradient alike reduced by its transmission T. A gradient adds to the
dispersion its square over the mean, so the excess falls as T: lowest where the
attenuation is deepest. With the intensity-to-background ratio constant, some
40, it is one picture: **nothing measured inside one box can see its
attenuation.**

Over the whole data set (`mxeq equivalents`, profile fitted, observations clean
otherwise): 597776 between 0.9 and 1.1, 14366 between 1.1 and 1.3, 786 between
1.3 and 2, 2 above; the higher bins 1 to 2 per cent below the rest. So the
dispersion flags nothing widespread -- and finds nothing at the backstop.

And in passing, ferritin's partials do not repeat insulin's: integrated with
sigma_m 0.0920 degrees, the partials suggest 0.0931, and are a few per cent low
rather than high -- insulin's sigma_m (items 1 and 3) does not generalise.

### What can, then

A box cannot see its own attenuation, so the answer must come from outside it:

1. **A mask of the backstop, given by the user**: a circle or polygon over its
   shadow and attenuated edge, as `dials.generate_mask` gives untrusted regions
   -- explicit, no heuristic, and the practice already. Its pixels treated as a
   module gap's are.
2. **`--d-max`**, which exists: cruder, a whole ring, and on ferritin it helped.
3. **In scaling, the physics: attenuation only lowers an intensity.** The true
   value is then the largest set of observations consistent with each other --
   (3,1,1)'s 13 near 950 -- where the attenuated spread from 10 to 590 and no
   value among them gathers as many. It would catch shadows nobody masked, but
   it assumes a low outlier is physical, to be tested on other data first.

Part 2, the outlier rejection, still stands: a wholly shadowed observation,
masked or not, must not outvote unharmed ones.

## Order, and how each is judged

1. **The background dispersion, measured**: the column, and on ferritin where
   (3,1,1)'s 27 observations, the flare's, and those on dense regions fall.
2. **Its test** if the measurement bears it out: judged by whether (3,1,1)'s
   12 attenuated observations are flagged -- those on the ramp should be; any
   that are not lie in shadow flat enough to be left to part 2 -- and its 13
   kept, by nothing flagged on flat backgrounds beyond the false-alarm rate, by
   `mxeq unique`'s inner shell against DIALS, and by Graeme's refinement.
3. **Outlier rejection**, with it: judged by a planted test (a group of 20
   consistent and 3 near zero with small sigmas: the 3 go, the 20 stay), by the
   wholly shadowed observations of (3,1,1) and insulin's (1,1,0) rejected, by
   the rejected counts against the present test's, and by refinement.
4. **The flare**, measured -- the dispersion will show where it is -- and
   changed only if the measurement says so.

Each default stays as it is until the change has been judged, and each is
byte-identical when off.
