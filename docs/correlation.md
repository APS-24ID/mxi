# Clustering data sets by their intensities: mxi_correlation

STATUS: stages 1 and 2 built, 15 October 2026 -- the correlation matrix and the
cosym coordinates. Stages 3 and 4 to follow.

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

## Stage 2: the coordinates (built)

Minimise f = 1/2 sum_ij w_ij (r_ij - x_i . x_j)^2 over the whole matrix, its
gradient -2 X (W.R - W.(X^T X)), for coordinates x_i in d dimensions, by
L-BFGS: ten pairs remembered, a backtracking line search, and as the starting
inverse Hessian at every iteration 1 / the curvature, c_ki = 2 sum_j x_kj^2
w_ji, as DIALS's minimiser (scitbx's lbfgs, diag_mode "always"); stopped at 100
iterations or |g| under 1e-5 of max(1, |x|). The start is DIALS's exactly:
numpy's np.random.seed(230); np.random.rand(...), drawn the same way from the
Mersenne twister.

The dimension as dials.correlation_matrix chooses it -- not by DIALS #3248's
projection, which only dials.cosym uses (`--by-projection` has it): a
minimisation at each dimension from 1 to min(n, 50), each scored with outlier
rejection (the elements w (r - x.x)^2 below q2 + (q3 - q1), halved), and the
elbow of that curve, as distl's resolution method 1, at least 2. Then the
minimisation at that dimension, and the coordinates' principal components,
centred as sklearn's PCA, with the variance each explains.

Two mistakes on the way, both in reading DIALS: its scan's
`min(20, max_calls) if max_calls else max_calls` limits a dimension's
minimisation to 20 calls only where max_calls is set, which by default it is
not -- read as 20, the curve was rough and the elbow 14 or 16 of 18 for three
planted groups; and without the curvatures a minimisation got almost nowhere
in those calls. Now three planted groups of six: the curve 0.0072, 0.0044,
2.7e-7 and flat, the elbow 3, the variance 0.51, 0.49, 0.00, each group its own
place. Checked: the start numpy's to the bit; the gradient finite
differences'; the minima scipy's L-BFGS-B's from the same start, 0.00973797
against mxi's 0.00973796 in 2 dimensions, 6.78e-7 against 6.72e-7 in 3.
Insulin's four sweeps of one crystal: 2 dimensions, the first component
everything -- one group.

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
