# Scaling

STATUS: `mxi_scale` scales one sweep, and several of one crystal together
(5 October 2026). Not yet: free-set validation. `docs/outstanding.md` has the
whole list.

**Several sweeps** -- or several crystals, each an experiment of its own
(`docs/multi-lattice.md`) -- each have a model of their own -- scale, decay and
absorption from each sweep's own width of rotation, as dials.scale's physical
model has one a dataset -- sharing the merged intensities; an observation's id
chooses its sweep's parameters. The one overall scale and the one overall B the
merged intensities can absorb are fixed over every sweep's points together, not
each sweep's apart, so the sweeps' relative scale is kept: it is what scaling
them together is for. The report gives each sweep's scale and relative B. A
sweep with no observations is refused (`docs/multi-sweep.md`).

`mxi_scale` puts every observation of a sweep on one scale, refines an error
model, and reports merging statistics -- the physical model of Beilsten-Edmands
et al. (2020), Acta Cryst. D76, 385-399, with cubic B-splines for the scale and
decay and the paper's spherical harmonics for the absorption surface.

## Running it

    mxi_scale symmetrized.rflx      # scaled.rflx

or from and to the DIALS pair, as `docs/rflx.md` describes:

    mxi_scale symmetrized.expt symmetrized.refl -o scaled.refl --output-expt scaled.expt

The space group is the crystal's, as `mxi_symmetry` wrote it (`docs/symmetry.md`).
`--space-group` and `--change-of-basis` override it for data that have not been
through `mxi_symmetry`; `--d-min`, `--d-max` (a low resolution limit: beside a
backstop, see `docs/backstop.md`), `--no-absorption`, `--profile-only` and
`--shells N` (20) do what they say. `mxi_scale --help` is the authority.

## The model

Inverse scale g = C(r) . T(t, d) . S(s0, s1), as the paper's physical model.

* **Scale C(phi).** A cubic B-spline in rotation angle, where DIALS uses a
  Gaussian-weighted average of the nearest three parameters. Control points 10
  to 20 degrees apart, as the paper finds, and fewer on a narrow sweep: DIALS
  used six on 30 degrees, and most were poorly determined. It is the crystal model's
  clamped uniform B-spline, `spline_weights` in `src/derivatives.cc`.
* **Decay T = exp(B(t) / 2d^2).** B(t) a cubic B-spline in t, image number for
  now as a stand-in for dose, with DIALS' weak restraint of sum B_i^2 toward
  zero. The zero-dose work would belong here: t as dose, and the extrapolation a
  use of the fitted B(t).
* **Absorption S = 1 + sum P_lm [Y_lm(s1) + Y_lm(s0)] / 2** in the crystal frame,
  as the paper's equation 7: l_max = 4 by default, 24 parameters, odd terms
  included since they absorb miscentring, with a restraint of w sum P_lm^2 toward
  zero. `--absorption-level` chooses the degree and w as dials.scale's does: low,
  the default, for about 1 per cent relative absorption -- degree 4, w 5e5;
  medium, about 5 per cent -- 6, 5e4; high, over 25 per cent, long wavelengths or
  heavy absorbers -- 6, 5e3. (w was 1000 until 7 October 2026, 500 times weaker
  than dials.scale's default, for no recorded reason. Forced on the 30 degree
  sweep with `--l-max 4`, the weaker the restraint the better the internal
  statistics -- I/sigma 20.2 at 1000, 19.8 at 5e3, 19.1 at 5e5 -- as more freedom
  always gives; which is right is a question for a wide sweep and an absolute
  measure.) Off for narrow sweeps -- below 60 degrees -- as DIALS turns it off,
  since there is not the angular coverage to determine it. `--l-max L` sets the
  degree, dials.scale's `lmax`: L (L + 2) terms, 48 at 6; 0 for none; asked
  for, it is used whatever the sweep's width, and it contradicts
  `--no-absorption`. More terms fit better and can fit noise: on the 30 degree
  sweep degree 4 and 6 both take Rmeas from 0.041 to 0.039, and degree 6 the
  error model's a to 0.90 -- the signature of a model absorbing noise, which a
  narrow sweep cannot rule out. Free-set validation, not yet built (item 22), is
  what would say which degree a sweep supports.

One overall scale is degenerate with the merged intensities and must be fixed:
the mean of C, or one control point.

### Fitting

The target is the paper's equation 2, sum w (I - g <I>)^2 plus restraints,
with <I_h> in closed form (equation 3) for the current g. At about 70
parameters the normal matrix is tiny, so full-matrix Levenberg-Marquardt with
analytic derivatives runs throughout, where DIALS uses L-BFGS until the last
cycle. The Jacobian is the variable-projection one, in Kaufman's form, which
lets <I_h> follow the parameters; holding <I_h> fixed gives the right gradient
but too much curvature, and fits of fifty steps ended still creeping. The
cost is building J^T J, linear in observations; its inversion gives the
parameters' uncertainties, propagated to the scale factors' variances (below).

Around the fit, as the paper's figure 2:

1. **Reflection selection.** A random subset of symmetry-equivalent groups, at
   least 2000 groups and 50000 reflections (section 4.4), for optimisation; the
   model is applied to everything.
2. **Outlier rejection.** A reflection's normalised deviation from the weighted
   mean of its group, excluding itself, above 6 sigma (Evans 2006). Repeated as
   the model improves, retesting everything.
3. **Intensity combination.** I = w I_prf + (1 - w) I_sum with w = 1 / (1 +
   (I_sum / I_mid)^3), I_mid chosen from logarithmically spaced values by Rmeas.
4. **Error model.** sigma'^2 = a^2 [sigma^2 + (b I)^2], fitted by regression
   across bins of merged intensity, as dials.scale's default, on groups with
   <I> above 25 and <I / sigma^2> above 0.85 ("The error model's fit, by
   regression", below); then, where integration wrote each variance with no
   signal, split, a_bg for that part and a for the rest ("The error model,
   split"). It had been a from the slope of the central normal probability
   plot and b by the paper's equation 17, alternated -- until 6S0R. The deviations take the EXACT
   variance of a difference from a mean the observation is part of, not the
   paper's eqn 12 prefactor, and sigma^2 already carries the scale's uncertainty
   (both below).

Not yet built: **free-set validation**, Rmeas on a tenth of the groups held out,
the check for overfitting. dials.scale warns, on the 300 image sweep, that 64
per cent of its eleven parameters have sigma over half their value.

## What it writes

A scaled, unmerged reflection table with dials.scale's columns --
`inverse_scale_factor` and its variance, `intensity.scale.value` and variance --
and its flags (scaled, bit 26; an outlier in scaling, bit 23; excluded, bit 24)
ADDED to those integration set, and dials.merge takes it -- as the DIALS pair,
from `scaled.rflx` by `mxi_convert scaled.rflx --drop-sweeps`. Without
`--drop-sweeps`, and in the pair `-o scaled.refl --output-expt scaled.expt`
writes, each sweep's crystal-less experiment is carried
(`docs/multi-lattice.md`, stage 2), which DIALS's programs past indexing may
not expect: not yet tried.
The intensity's variance already carries the scale's uncertainty, as dials.scale's
does, so nothing downstream should add `inverse_scale_factor_variance` to it
again. The models, reindexed, go with it, into `scaled.rflx` (`scaled.expt` in the
pair). It prints the model, the error
model, merging statistics in twenty shells and dials.scale's summary. HTML
reports are for a later mxi_report, MTZ and mmCIF for a later mxi_export.

## Against dials.scale, on the 300 image sweep

Integrated here and scaled here, against dials.scale on DIALS' own integration,
overall:

    ours   20007 obs  7183 uniq  mult 2.79  71.4 %  <I> 131.7  <I/sI> 18.9
           Rmerge 0.034  Rmeas 0.041  Rpim 0.022  CC half 0.987
    DIALS  19737 obs  7268 uniq  mult 2.72  71.0 %  <I> 130.0  <I/sI> 31.4
           Rmerge 0.038  Rmeas 0.045  Rpim 0.024  CC half 0.999

* **The error model.** a = 1.017, b = 0.024 here; dials.scale reported 0.531 and
  0.0246. dials.scale normalises its deviations with eqn 12's sqrt((n-1)/n)
  (`calc_deltahl`), which leaves their spread at (n-1)/n and a low by that
  factor -- a planted 1.3 comes back 1.158 in groups of eight and 0.662 in
  pairs -- and this sweep's multiplicity is 2.7. The deviations here take the
  exact variance of a difference from a mean the observation is part of.
  That one factor accounts for most of the difference in <I/sI>: 18.9 scaled by
  1.017 / 0.531 is about 36.
* **The fit.** Levenberg-Marquardt with the variable-projection Jacobian in
  Kaufman's form converges in three or four steps. Holding <I> fixed instead --
  right gradient, curvature overstated -- left three fits of fifty steps each
  still creeping. With the curvature right, a constant offset in B, degenerate
  with the merged intensities but for small differences of d within a group,
  wandered to -2 A^2; B is now centred on its mean, as the scale is on its.
* **CC half in the lowest shell,** 0.975 against 0.999: a few groups, and
  without the three widest 0.9990. What they share is two things scaling has
  found in the integration, both judged by symmetry equivalents: observations
  with partiality below 0.9 are 10.5 per cent high, as are, by 4 to 5 per cent,
  those in the first two and last five images -- partials over-corrected -- and
  those with less than 95 per cent of the profile measured, near a module gap,
  3.4 per cent low. These are integration's to fix.

## The error model against dials.scale's, judged on the same observations

Two crystals of insulin (`docs/multi-lattice.md`), mxi's integration scaled
by each program: dials.scale's error model a = 0.812, b = 0.0301 on some 4670
reflections; mxi's a = 1.100, b = 0.0223 on 4791. Merged mean I/sigma 12.5
and 9.4. The parameters are not comparable across programs -- each applies
its own to what its own combination and scaling hand it -- so the test is of
each program's final sigma on the same observations: the variance of each
observation's EXACT deviation from its group's mean, which is one for sigma
that are right, binned by mean intensity -- not by I/sigma, which depends on
sigma and biases the bins:

    mean intensity sextile, weak to strong    1     2     3     4     5     6
    mxi                                      0.26  0.58  0.84  0.95  0.95  0.98
    dials.scale                              0.51  1.11  1.55  1.71  1.69  1.70

For the stronger two-thirds mxi's sigma are right and dials.scale's some 30
per cent small: its lower a, which is what its higher mean I/sigma comes from,
does not calibrate them better -- by its own eqn 12 deviations still 1.43 on
the strongest. For the weakest third both are too large, mxi's by up to a
factor of two in variance, and no a and b can mend both ends at once. That is
integration's: matched reflection by reflection, mxi's integration sigma for
I/sigma under one are 0.68 to 0.83 of DIALS's own integration's, and the two
agree above ten -- so DIALS's own pipeline starts from larger ones still. The
weak reflections' variance, its background term above all, is the thing to
look at; mxi's error model is kept as it is.

It was the fitted variance (`docs/integration.md`, "The fitted variance"): the
background's term taken as a sum's, and a pixel's least variance above these
crystals' whole background. Both mended, the same two crystals scale with a =
1.132, b = 0.0189, and by mean intensity sextile 0.42 0.80 1.06 1.05 0.93
0.98; insulin's sweep 0.95 1.06 1.09 1.13 1.23 1.20, where it was 0.66 0.92
1.04 1.13 1.27 1.23. CC half as it was.

## Partial reflections (15 October)

Two corrections, both as dials.scale, found on Graeme's 36 insulin crystals
(dials.symmetry's output, scaled by dials.scale and by mxi_scale):

* **Partiality divides the summed intensity only.** A profile-fitted
  intensity estimates the whole reflection whatever part of it was recorded;
  mxi divided it by the partiality too, making every partial one too large
  by 1 / p. dials.scale (combine_intensities): Ipr lp / qe, Isum lp / (qe p).
* **A reflection split by integration is one observation.** DIALS's
  integration cuts a reflection crossing its blocks of images into pieces
  sharing a partial_id; dials.scale combines them first
  (sum_partial_reflections): pieces of partiality under 0.99 with positive
  variances, within an experiment; the profile-fitted value a mean weighted by
  I^2 / sigma^2, the summed values and variances added, the partialities
  added; the partiality cutoff applied after. mxi_scale scaled each piece
  apart -- 23170 of them there, 8 per cent of the observations of raw I/sigma
  over 100. It now combines them and says how many: 23170, exactly the rows
  dials.scale's output has fewer.

Effect there: the error model's b 0.067 to 0.056 (dials.scale's 0.026), R
merge 0.237 to 0.225, at low resolution 0.190 to 0.176 (0.192, 0.127), CC half
0.992 to 0.997 (0.999); data set 12's relative B -10.1 to -5.7 A^2 became -5.8
to -4.0. On insulin as mxi integrates it, with no pieces: R merge 0.035 to
0.034, <I/sigma> 17.5 to 17.7. The correlations between those crystals,
sigma-weighted, rose only a little (within a species 0.965 to 0.969, between
0.863 to 0.870, dials.scale's 0.982, 0.916); unweighted they are near
dials.scale's (0.991, 0.956 against 0.994, 0.967). So what remains is the
sigmas, by the error model's larger b, which caps I / sigma at 1 / (a b):
outstanding item 69.

## The error model's fit, by regression

sigma'^2 = a^2 (sigma^2 + b^2 I^2) is fitted as dials.scale's default does it:
in each of twenty intensity bins -- by the reflection's merged intensity, not
the observation's own -- the deviations' mean square should be a^2 (P + b^2 Q),
P and Q the variance as given and the merged intensity squared, each with the
exact variance of a deviation from a mean it is part of; linear in a^2 and
a^2 b^2, so fitted across the bins by weighted least squares, each bin's weight
its count over its mean square squared; b at zero if the fit wants it below.
a had been the slope of every observation's normal probability plot, b then
fitted about it: on Graeme's 6S0R set, its strong reflections out by a tenth of
their intensity (raw calibration 86 for the strongest sixth, 1.9 for the
weakest), that gave a = 4.81 and b = 0.023 -- the weak reflections' sigmas
three and a half times too large (final calibration 0.08; the anomalous slope
0.33 and its control 0.34). Binned by each observation's own intensity, a
planted a = 1.0, b = 0.15 came back a = 1.18; by the merged intensity, 1.0 and
0.15 (tests/test_scale.cc). Insulin's sweep: a 1.13 to 1.23, b 0.023 to 0.018,
the calibration by resolution 0.95 to 1.07 for the most part, the strongest
sixth by intensity 1.10 to 1.21; the two crystals as they were.

## The error model, split

Where integration wrote each variance with no signal, the error model is
split: a_bg^2 V_bg + a^2 (V_signal + the scale's term) + a^2 b^2 I^2 -- one a
for a background's variance, another for the signal's -- fitted by the same
regression, three columns, over every group of two or more in twenty bins of
merged intensity, equal counts, the weak groups (which the single model's fit
leaves out below <I> 25, as dials.scale) deciding a_bg. Kept if all three come
out positive; a_bg = a is the single model. One a for both had carried 6S0R's
strong reflections' excess, a tenth of their intensity, into its weak ones'
sigmas: final calibration 0.05 to 0.08 (outstanding 62, 64). Planted, a_bg 1
and a 3 come back as such, and one a of 1.3 as one a (tests/test_scale.cc).
Insulin's sweep: a_bg 1.08, a 1.21; calibration by intensity 0.85 to 1.21
before, 0.96 to 1.08; by resolution 0.84 to 1.25, 0.91 to 1.12. The two
crystals: a_bg 0.91, a 1.11; the weakest sixth 0.84 to 1.01, the anomalous slope
and its control 0.97 and 0.96 to 1.00 and 0.99. Under --verbose the fit is
printed a bin at a time.

## The scale model, written

`mxi_scale` writes its model into the scaled experiments, as dials.scale writes
its `scaling_model` (outstanding 26): an entry a lattice, the experiment's
`scaling_model` its index. mxi's own type, `mxi_physical`: the same three
corrections as DIALS's physical model, but cubic B-splines where DIALS smooths
with Gaussians, so its numbers would be read wrongly as DIALS's. Each entry:
`scale`, `decay` and `absorption` (with `lmax`), each `n_parameters`,
`parameters` and `est_standard_devs` (from the covariance); the run's
`error_model` (`a`, `a_background`, `b`, `split`); and `configuration`, which
says in words how to evaluate g: t the observation's place in its scan from
`xyzcal.px`, the spline's segments and weights, the harmonics' order and frame.
Evaluated from the entry alone, g matches every scaled reflection's
`inverse_scale_factor` to 1e-9 (tests/test_scaling_model.py). `mxi_convert`
leaves it out of DIALS files, which DIALS could not read, and a program told
to write the DIALS pair writes none.

## The calibration

After the merging statistics, how well the sigmas describe the scatter: each
observation's deviation from the mean of its equivalents without it, over the
sigma of the difference -- 1 where the sigmas are right -- its variance in ten
shells of resolution and six of mean intensity, equal counts a bin; with the
variances before the error model ("raw") and after ("final"); and each as a
mean square and robustly, (1.4826 MAD)^2, so a few wild observations show as
the two apart. Printed above the summary of merging statistics; in the
.rflx's summary (`calibration_by_resolution`, `calibration_by_intensity`);
drawn by mxi_report, plain and robust, raw and final, on a log scale. The two insulin crystals: final 0.95 to 1.03
nearly everywhere, plain and robust alike. Insulin's sweep before the error
model: the strongest sixth 7.8 plain, 2.3 robust -- a tail of a few -- 1.10 and
0.76 after. Written for Graeme's 6S0R set, whose error model came out a = 4.8
while its anomalous slope's control said the sigmas were three times too large.

## The anomalous slope and its control

The anomalous slope -- the normal probability plot's slope of dI/sigma(dI) over
the Friedel pairs, |x| < 0.9, as dials.scale -- is above 1 when there is
anomalous signal or when the sigmas are too small, and alone cannot say which.
Graeme's ten-turn insulin went from 1.04 to 1.52 when scaling's variances came
from the model. So beside it, its control: the same slope over random halves
of each pair's observations, I+ and I- pooled and dealt out again in their
numbers -- the same counts and weights, no signal. Near 1 when the sigmas are
right; the slope above it is the signal; both above 1, the sigmas too small.
Printed under the slope ("random halves"), in the summary
(`anom_slope_random`), the report and the regression figures. Insulin's sweep:
1.003 and 1.006; the two crystals: 0.969 and 0.961 -- no signal in either, and
the sigmas right.

## The normal equations, a group at a time

An observation's row in the fit is d_i - g_i G, d_i its own few parameters'
derivatives (scale, decay, absorption: some thirty) and G those of its group's
merged intensity, which touches every parameter any of the group's observations
does. Summed over a group, sum w m^2 (d_i - g_i G)(d_i - g_i G)^T = m^2 (sum w
d_i d_i^T - den G G^T), sum w g_i d_i being den G: each observation's own
parameters, then G G^T once a group. Formed an observation at a time, every row
was dense at high multiplicity -- Graeme's 36000 images, multiplicity 374 over
ten turns, 448 parameters: three fits of 22 to 30 s each, twice over with
--d-min-auto, of mxi_scale's 188. The same sum rearranged: the inverse scales
agree with the old to 2e-15, insulin replicated twentyfold and the two crystals.
Graeme's 36000 images then: each fit 0.3 to 0.4 s where it had been 22 to 30,
mxi_scale 188 s to 22.5. The scale's uncertainty, 3 s of what remained, inverted
its constrained normal matrix a column at a time with a factorisation each, m^4
/ 3 at 446 parameters; one Cholesky factorisation now, the columns solved from
it, the same numbers to the bit; and each observation's variance from the
covariance with one gradient buffer a thread, not a vector an observation.

## The cell from two-theta

The cell's e.s.d.s are scaling's, as dials.two_theta_refine gives them: the
cell alone, refined against the scattering angle of every scaled reflection's
observed centroid (I/sigma 3 or more, not an outlier) -- an angle that depends
on the cell, the wavelength and where the spot fell, not on the orientation --
one cell over every crystal, constrained to the space group's crystal system,
its e.s.d.s from (J^T W J)^-1 chi^2 / (n - p) (`src/two_theta.cc`). Each
centroid's variance at least a pixel's own, size^2 / 12: counting statistics
alone put strong spots' at a tenth of a micron, and those, mostly at low
resolution where an angle says least about the cell, took 43 per cent of the
weight -- insulin's a 77.931 A, chi^2 89 times its degrees of freedom; floored,
77.9109(3), the median of each reflection's own (77.9106) and chi^2 2.6. The
e.s.d.s are statistical: an error in the distance or the wavelength is not in
them. Reported, printed and in the summary and report; the crystals' own cells
are not changed by it. Refinement's covariance e.s.d.s, which fitted the
orientation, the detector and the cell together on the strong spots, are gone
(Graeme). Planted cells, cubic and monoclinic, are recovered, and the cubic
e.s.d. is its scatter over 200 trials to 20 per cent (`tests/test_two_theta.cc`).

## Variances from the model

An observation's variance from integration moves with its own noise -- a low
one gets a small variance and a large weight -- so an inverse-variance
weighted mean leans low, and at a high multiplicity of weak data the merged
I/sigma goes negative: Graeme's ten-turn insulin, -5.7 in the outer shell at a
multiplicity of 210 (`docs/integration.md`, "The fitted variance"). Once a
model exists, each observation's variance is built on it instead: its
variance with no signal, which integration writes, plus the signal's Poisson
variance at the model's intensity g I_h, corrected as the observation is --
the same model for every observation of a reflection. The profile and summed
variances with no signal are mixed as the variances are; the error model's b
term and its intensity bins take the model's intensity too. Done after the
first fit and again after each, the last before the final error model. A
table integrated before, without the columns, scales as it did.

Each reflection's weighted mean against its plain mean, in its merged sigma,
by resolution sextile: the two crystals -0.89 -1.73 -2.37 -2.23 -1.66 -0.96
before, -0.10 -0.08 -0.06 -0.02 -0.01 +0.03 after, merged I/sigma in the outer
two 0.4 and 1.4 to 1.0 and 2.4; insulin -0.05 to -0.17 to within 0.06. The
calibration by mean intensity sextile, the two crystals 0.42 0.80 1.06 1.05
0.93 0.98 to 0.79 1.01 1.02 0.98 0.98 0.99; insulin 0.99 to 1.21 as it was.

`mxi_scale --no-model-variances` keeps each variance as integrated, for
comparison: on insulin the two agree, R merge 0.035 and CC half 0.987 either
way.

## The scale's uncertainty, before the error model

The covariance of the parameters is the inverse of the normal matrix of the
variable-projection Jacobian and the restraints, under the two constraints the
model's normalisation fixes -- the scale's mean at one, the relative B's at
zero -- and scaled by the goodness of fit, whose degrees of freedom count every
merged intensity as a parameter. Tested against the thing it claims: the
geometry fixed and the noise drawn forty times, the spread of the fitted g at
five places across the scan over the variance predicted is 1.023, each between
0.89 and 1.19.

Each observation's var(g) goes into its variance as I^2 var(g) / g^2 BEFORE the
error model is refined and applied, so that a and b correct what remains after
it. dials.scale does the reverse, and applies it as a factor (1 + sigma_g / g)
on the variance after the error model: linear in the fractional error, and the
same at every intensity. On the 300 image sweep the scale is well determined --
sigma(g)/g 0.22 per cent on average, 0.68 at most, three times as large at the
ends of the scan -- and b moves only from 0.0237 to 0.0236: the term goes as
I^2, as (b I)^2 does, and a b fitted without it spreads the uncertainty of the
worst-determined parts of a scan over every observation.

## The summary, and what the error model does to it

mxi_scale ends with dials.scale's summary -- overall, lowest shell, highest --
with iotbx.merging_statistics' definitions, which dials.scale reports: R factors
with Friedel mates merged and apart, anomalous completeness against the complete
acentric set, anomalous multiplicity over the groups with mates apart, the
correlation of I+ - I- between random halves, the slope of the normal
probability plot of dI / sigma(dI) over |x| < 0.9, dF/F as
sqrt(2 <(F+ - F-)^2> / <F+^2 + F-^2>), and dI/s(dI) as mean |dI| over mean
sigma(dI). On the 300 image sweep, overall:

                              here     dials.scale
    Completeness              71.4     71.0
    Rmerge(I), (I+/-)         0.034    0.038, 0.032
                              0.031
    Rpim(I), (I+/-)           0.022    0.024, 0.029
                              0.027
    Anomalous completeness    51.2     50.2
    Anomalous multiplicity    1.7      1.7
    Anomalous slope           0.955    1.890
    dI/s(dI)                  0.787    1.476

What depends on the intensities alone agrees. The two that differ are ratios to
sigma, and dials.scale's sigmas are shrunk by its a of 0.531: ours scaled by
1.017 / 0.531 are 1.83 and 1.51. An anomalous slope of 1.89 reads as anomalous
signal; this sweep's is 0.955 -- none measurable, as insulin's sulphur at
0.95 A would lead one to expect -- and the 1.89 is the error model's bias.

## The resolution limit

As `dials.estimate_resolution` finds it, read from its source
(`dials/util/resolution_analysis.py`): CC half in 50 bins of equal count of
unique reflections, Wilson outliers (E^2 of 16 or more) left out; a tanh in
d*^2, CC half = (1 - tanh((d*^2 - s0) / r)) / 2, fitted by weighted least
squares with sigma 1/sqrt(n - 3) from r = 0.2 and s0 = 0.4; and the limit where
the fitted values cross 0.3. The second limit is where CC half stops being
significant: a bin is significant where CC half exceeds t / sqrt(n - 2 + t^2),
t Student's at the upper 0.1 on n - 2 degrees of freedom, and a logistic
1 - expit(r (d*^2 - res)) through which bins are gives 1/sqrt(res). The summary
gains dials.scale's "Suggested" column: the same statistics cut at the CC half
limit -- unless `--d-min` was given, when the cut is already chosen and the
column is left out; the estimate itself is still printed.

`--d-min-auto` makes the cut itself: everything is scaled once, the limit is
estimated, and the sweep is scaled again to it -- again, not filtered, since the
outliers, the error model and the fit all depend on what is in. The limit is
rounded to a hundredth of an angstrom and reported as the `--d-min` it is
equivalent to, and the table is that run's, byte for byte
(`python/tests/test_d_min_auto.py`). It cuts only if the fitted CC half comes down
to the limit within the data: where it never does, as on the 300 image sweep at
0.3, the estimate is only the last bin with pairs, a cut there would drop good
data, and nothing is applied. `--cc-half-limit C` sets the CC half both the
estimate and the cut use (0.3); `--d-min` and `--d-min-auto` together are
refused.

One difference, deliberately: only bins with pairs enough to fit (n > 3) decide
whether CC half is above the limit everywhere. A detector's corners hold
reflections seen once, with no halves and a CC half of zero, and under
dials.estimate_resolution's rule a sweep whose CC half stays near one to its
last bin was given no limit at all. On the 300 image sweep that is the case: the
limit is the last bin with pairs, 1.73 A, where the data are 83.5 per cent
complete with CC half 0.987, and dials.scale suggested 1.63, the edge. The
significance limit is 1.68.

On a sweep of 1082994 reflections, dials.estimate_resolution run on this
program's own scaled.refl and this program agree exactly: 1.25 A by CC half and
1.17 by significance.

## Threads

`--threads N` (the machine's by default). The fits' target and normal
equations, outlier rejection, the error model's deviations, the scales'
uncertainties and the absorption harmonics run in parallel; the merging
statistics, whose random half-splits depend on order, do not. The work is cut
into 64 blocks whatever the count and combined in block order, so the answer is
the same on any number of threads: scaled.refl is byte-identical on one and four,
and a test holds a fit, outlier rejection, the error model and the covariance to
the bit across counts, failing if the blocks follow the threads.

Only the normal equations' order of summation differs from the serial version's,
by rounding: on the 300 image sweep the scale factors move by 3.4e-11 relative
and their variances by 3.7e-9, the flags and the error model not at all. The
error model's deviations are computed a group at a time, each group's spread in
its own order, and are bit for bit what one thread gives; the normal quantiles
of its probability plot are computed once rather than every round.

## The error model's search, made cheaper

The error model was most of scaling even with threads: 1.35 to 1.49 s of each
pass on the MacBook with 16, against 0.19 to 0.22 for the three fits. Three
changes, the answer moving only in its last digits. The thread pool keeps its
workers rather than making them every call (`src/parallel.cc`). A deviation's
variance is a^2 (P + b^2 Q) with P and Q summed once a fit -- the weights and the
merged means do not change within it -- so an evaluation is a square root an
observation, where every group's spread was summed again. And b is found by
Brent's method to 1e-8, where a golden section of sixty steps went to 1e-13 and
the rounds stop at a change of 1e-6. On one thread on the 300 image sweep the
three error models went from 0.61 s to 0.14; the error model (a 1.0092, b
0.0244), Rmerge, Rmeas, Rpim, CC half, I/sigma, completeness and the anomalous
correlation are the same to every digit printed.

On the MacBook, the 16M sweep with --d-min-auto: 5.3 s to 2.9, and the chain
from images to scaled data 48 s to 43. Gathering the observations was then the
largest phase, 0.40 to 0.45 s a pass, serial: eight lookups of a column by name
a reflection, the group's symmetry operators three times a row, and an ordered
map. The columns are looked up once, a row's symmetry is worked out in parallel,
and the groups numbered in row order through a hash -- scaled.refl
byte-identical, and a third of the time on one thread.

## Friedel mates apart, --anomalous

By default an acentric reflection's I(+) and I(-) are one group, as they are
in dials.scale. With a strong anomalous signal they genuinely differ, and the
difference is then scaling's to explain -- which it does in the worst way.
Planted, with noise at a = 1.3, b = 0.03, four observations of each mate and
3000 reflections:

| Friedel difference | merged: a | merged: outliers | apart: a | apart: outliers |
| --- | --- | --- | --- | --- |
| none | 1.35 | 397 | 1.36 | 338 |
| +-5 % | 1.63 | 1688 | 1.37 | 334 |
| +-10 % | 3.24 | 4509 | 1.36 | 323 |
| +-20 % | 9.80 | 8102 | 1.36 | 330 |

of 24000 observations. Merged, the error model takes the difference as noise --
through a, not b, every sigma inflated -- and outlier rejection throws away the
observations with the largest anomalous differences: a fifth of them at +-10 per
cent. Kept apart, the planted model comes back whatever the signal.

`--anomalous` makes an acentric reflection's I(-) a group of its own, so
scaling, the error model and outlier rejection all see the mates separately,
as dials.scale's `anomalous=True`. Centric reflections stay one group:
`friedel_plus` calls every centric reflection I(+), their mates equal by
symmetry. The merging statistics and the resolution estimate keep the Friedel
pairs -- `ScaleData::pair`, numbered as the groups are without the flag -- since
they split each pair into I(+) and I(-) themselves: the merging table reports the
same reflections either way, its anomalous columns what such data wants read.
The report says how many groups were scaled from how many reflections. Without
the flag, every table is byte-identical to before it.

## Decisions

* **gemmi** for space groups, pinned at a release as a submodule, behind
  `src/symmetry.hh`; only its symmetry source is compiled.
* **Reports** of merging statistics go in a later mxi_report; **MTZ and
  mmCIF** in a later mxi_export. `mxi_scale` prints its tables, as the other
  programs do.
