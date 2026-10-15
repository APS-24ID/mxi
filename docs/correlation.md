# Clustering data sets by their intensities: mxi_correlation

STATUS: stage 1 built, 15 October 2026 -- the correlation matrix. Stages 2
to 4 to follow.

`mxi_correlation` is to do what dials.correlation_matrix does (Thompson et al.
2025, Acta Cryst. D81, 278-290): from a scaled file whose lattices are each a
data set -- many crystals, as in Graeme's COWS_PIGS_PEOPLE tutorial -- find
the groups of data sets that belong together, by their intensities, where
unit cells cannot tell them apart.

## Stage 1: the correlation matrix (built)

Each data set is a lattice experiment. Its scaled reflections -- scaled, not
outliers, partiality at least 0.4 (`--min-partiality`), intensity.scale.value
and its variance -- are merged by symmetry-unique index, Friedel mates
together, as cctbx's merge_equivalents by default: one observation as it is;
several, the inverse-variance weighted mean, its sigma sqrt(max(V_w / n,
1 / sum w)), V_w GSL's weighted variance (use_internal_variance, cctbx's
default). Every pair's correlation over the reflections both have is DIALS's
weighted_cchalf: each reflection weighted 1 / (s_i^2 + s_j^2), normalised; r
the weighted covariance over the weighted standard deviations; Kish's
effective sample size n_eff = 1 / sum w^2. The cosym target's weights are
DIALS's default `count` weights: n_eff, rescaled to sum to the number of pairs
that have one; none on the diagonal or for a pair with fewer than three
reflections in common (`--min-pairs`). Pairs a row a thread.

Checked: an independent Python implementation -- gemmi's asymmetric unit,
cctbx's merging formula and DIALS's correlation written afresh -- gives the
same merged counts and every r to 3e-15, on insulin as four sweeps. Insulin's
four sweeps of one crystal: r 0.993 to 0.994; the two insulin crystals of
ins2_27: 0.980. Unit tests: the pair correlation against the formula written
out, the weights' rescaling, two planted groups told apart
(tests/test_correlation.cc).

## Stage 2: the coordinates (next)

Minimise sum w_ij (r_ij - x_i . x_j)^2 for coordinates x_i in d dimensions.
The dimension as DIALS has found it since September 2026 (#3248): one
minimisation in the largest dimension tested, min(n, 50), the coordinates
rotated to their principal axes, each smaller dimension's functional read off
the leading axes, with the outlier rejection DIALS uses for this; the elbow of
that curve the dimension. The coordinates then rotated to principal axes, the
variance each explains reported.

## Stage 3: OPTICS (to follow)

OPTICS on the coordinates, as scikit-learn's with xi = 0.05 and min_samples
s_min = max(5, b n / d), b = 0.5 (the paper's eqn 9). DIALS now also chooses
min_samples by optimising a silhouette score; to follow if the simple rule
proves not enough. Written in C++, checked against scikit-learn in the Python
tests.

## Stage 4: output

The clusters and outliers printed and kept in the summary; each cluster
optionally written as its own .rflx for mxi_scale; mxi_report draws the
correlation heat map ordered by cluster, the reachability plot, the principal
components and the dimension functional, and Ward dendrograms of r and of the
cosine angles, uncut.
