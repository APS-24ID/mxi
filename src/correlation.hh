// Intensity-based clustering of isomorphous data sets, as
// dials.correlation_matrix (Thompson et al. 2025, Acta Cryst. D81, 278-290):
// the correlation of every pair of data sets, each weighted by its intensities'
// uncertainties, and the effective number of reflections behind it. Stage 1 of
// mxi_correlation; the cosym coordinates, their dimension and OPTICS follow
// (docs/correlation.md).

#ifndef MXI_CORRELATION_HH
#define MXI_CORRELATION_HH

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include "expt.hh"
#include "refl.hh"
#include "symmetry.hh"

namespace mxi {

//: One data set's merged intensities: symmetry-unique indices, Friedel mates
//: merged, ascending; the inverse-variance weighted mean of its observations
//: and that mean's sigma, as cctbx's merge_equivalents.
struct MergedSet {
  std::vector<Miller> hkl;
  std::vector<double> intensity, sigma;
};

//: Each lattice experiment's scaled reflections merged: scaled, not outliers,
//: partiality at least `min_partiality` (dials.correlation_matrix's 0.4),
//: intensity.scale.value and its variance, in `group`. One set an experiment
//: with a crystal, in experiment order; `which` gets their indices.
std::vector<MergedSet> merged_by_experiment(const ExperimentList &experiments,
                                            const Table &reflections,
                                            const SpaceGroup &group,
                                            double min_partiality,
                                            std::vector<std::size_t> *which);

//: Two data sets' correlation over the reflections they share, as DIALS's
//: weighted_cchalf with no variance floor: each pair weighted 1 / (s_i^2 +
//: s_j^2), the weights summing to one; r the weighted covariance over the
//: weighted standard deviations; n_eff = 1 / sum w^2, Kish's effective sample
//: size. Not valid with fewer than `min_pairs` in common, or no spread.
struct PairCorrelation {
  bool valid = false;
  double r = 0.0;
  double n_eff = 0.0;
  std::size_t common = 0;
};
PairCorrelation weighted_correlation(const MergedSet &a, const MergedSet &b,
                                     std::size_t min_pairs = 3);

//: Every pair's correlation, r_ij, and its weight in the cosym target, w_ij:
//: DIALS's default `count` weights -- n_eff, rescaled so that the weights sum
//: to the number of pairs that have one -- zero on the diagonal and where a
//: pair has no correlation. Row major, n x n, symmetric. Pairs a thread at a
//: time.
struct CorrelationMatrix {
  std::size_t n = 0;
  std::vector<double> r, w;
};
CorrelationMatrix correlation_matrix(const std::vector<MergedSet> &sets,
                                     std::size_t min_pairs = 3);

//: Stage 2: coordinates for the data sets, x_i in d dimensions, minimising the
//: cosym target f = 1/2 sum_ij w_ij (r_ij - x_i . x_j)^2 over the whole
//: matrix (Gildea & Winter 2018), its gradient -2 X (W.R - W.(X^T X)), by
//: L-BFGS from coordinates uniform on [0, 1) drawn as numpy's
//: np.random.seed(seed); np.random.rand(...) draws them -- DIALS's start.
//: Coordinates dimension-major, as DIALS's: x[k * n + i] is data set i's
//: k-th coordinate. `max_calls` 0 for no limit.
double cosym_functional(const CorrelationMatrix &m,
                        const std::vector<double> &x, std::size_t dim);
std::vector<double> minimise_cosym(const CorrelationMatrix &m, std::size_t dim,
                                   int max_iterations, int max_calls,
                                   unsigned seed, double *functional = nullptr);

//: The functional at x for judging a dimension, as DIALS's: with outlier
//: rejection, only the elements w (r - x.x)^2 below q2 + (q3 - q1), quantiles
//: over the whole matrix, halved.
double dimension_score(const CorrelationMatrix &m, const std::vector<double> &x,
                       std::size_t dim, bool outlier_rejection);

//: The elbow of a decreasing curve, as distl's resolution method 1 (Zhang et
//: al. 2006) and DIALS's dimension choice: from the point with the least slope
//: to the last, the point furthest from the line joining them.
std::size_t elbow(const std::vector<double> &x, const std::vector<double> &y);

struct CosymOptions {
  std::size_t maximum_dimensions = 50; //: dimensions tested: min(n, this)
  std::size_t dimensions = 0;          //: 0: chosen by the elbow
  bool outlier_rejection = true;       //: in judging a dimension
  bool by_projection = false; //: one minimisation projected (DIALS #3248)
  int max_iterations = 100;
  //: A dimension's minimisation in the scan: no limit, as DIALS's by default --
  //: its min(20, max_calls) applies only where max_calls is set, which by
  //: default it is not (read once as 20 a dimension: the curve rough, the
  //: elbow wrong).
  int scan_max_calls = 0;
  unsigned seed = 230;
};

//: Stage 2 whole, as dials.correlation_matrix: the functional for each
//: dimension tested, the elbow (at least 2), the minimisation there, and the
//: principal components of the coordinates -- centred, as sklearn's PCA -- with
//: the variance each explains.
struct CosymCoordinates {
  std::size_t n = 0, dim = 0;
  std::vector<double> dimensions, functional; //: the curve the elbow is of
  std::vector<double> x;                      //: dimension-major, as minimised
  std::vector<double> components; //: n x dim, row major: each data set's
  std::vector<double> explained_variance, explained_variance_ratio;
  double final_functional = 0.0;
};
CosymCoordinates cosym_coordinates(const CorrelationMatrix &m,
                                   const CosymOptions &options);

//: Stage 3: OPTICS on the coordinates, as scikit-learn's (Ankerst et al. 1999),
//: Euclidean: core distances to the min_samples-th nearest point (itself the
//: first), infinite beyond max_eps; the ordering from the lowest index, each
//: next the unprocessed point of least reachability; reachability
//: max(distance, core distance) to points within max_eps; both rounded as
//: np.around to 15 decimals. Clusters by the xi method, with predecessor
//: correction, the minimum cluster size min_samples, labels -1 for noise.
//: Points row major, n x d.
struct OpticsResult {
  std::vector<std::size_t> ordering;
  std::vector<double> reachability, core_distance; //: by point
  std::vector<long> predecessor;                   //: by point, -1 for none
  std::vector<int> labels;                         //: by point, -1 noise
  std::vector<std::pair<std::size_t, std::size_t>> clusters; //: in ordering
};
OpticsResult optics(const std::vector<double> &points, std::size_t n,
                    std::size_t d, std::size_t min_samples, double max_eps,
                    double xi);

//: The mean silhouette of the points labelled (not -1), as scikit-learn's
//: silhouette_score: (b - a) / max(a, b), a point alone in its cluster 0.
double silhouette(const std::vector<double> &points, std::size_t d,
                  const std::vector<int> &labels);

struct ClusterOptions {
  double xi = 0.05;
  double max_eps = 0.5;            //: DIALS's max_distance
  double min_points_buffer = 0.5;  //: the paper's b
  bool optimise = true;            //: DIALS's optimise_input
  double alpha = 1.0, gamma = 1.5; //: the noise penalty
};

//: Stage 3 whole, as dials.correlation_matrix's cluster_cosine_coords: with
//: under 5 data sets OPTICS with min_samples n; otherwise the paper's
//: min_samples max(5, int(n / d b)), and, optimising (DIALS's default), every
//: min_samples from 5 to min(n, 50) scored silhouette x (1 - alpha f^gamma), f
//: the fraction labelled noise, -1 for one cluster or none -- the best kept
//: (DIALS samples the same range by optuna, 20 random and 30 TPE trials) --
//: and, the best below zero, DIALS's single-cluster heuristic.
struct SignificantClusters {
  std::vector<int> labels; //: by data set, -1 an outlier
  std::size_t min_samples = 0, initial_min_samples = 0;
  double silhouette = 0.0, score = -1.0;
  bool single_cluster = false;
  std::vector<std::pair<std::size_t, double>> scores; //: min_samples, score
  OpticsResult model;
};
SignificantClusters significant_clusters(const CosymCoordinates &c,
                                         const ClusterOptions &options);

//: Stage 4, the report: hierarchical clustering by Ward's method, as
//: scipy.cluster.hierarchy.linkage(method="ward") on a distance matrix
//: (n x n, row major): the Lance-Williams update on the distances, the two
//: nearest clusters merged each step, ties to the lowest pair. Each row as
//: scipy's: the two clusters merged (a leaf below n, a merge n + its row), the
//: distance, the size -- the smaller id first.
std::vector<std::array<double, 4>> ward_linkage(const std::vector<double> &d,
                                                std::size_t n);

//: A dendrogram of that linkage, as scipy's dendrogram draws it: the leaves
//: left to right, and for each merge the four points of its bracket, x at the
//: leaves' places 5, 15, 25 ..., y the heights.
struct Dendrogram {
  std::vector<std::size_t> leaves;
  std::vector<std::array<double, 4>> x, y;
};
Dendrogram dendrogram(const std::vector<std::array<double, 4>> &linkage,
                      std::size_t n);

} // namespace mxi

#endif // MXI_CORRELATION_HH
