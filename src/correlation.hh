// Intensity-based clustering of isomorphous data sets, as
// dials.correlation_matrix (Thompson et al. 2025, Acta Cryst. D81, 278-290):
// the correlation of every pair of data sets, each weighted by its intensities'
// uncertainties, and the effective number of reflections behind it. Stage 1 of
// mxi_correlation; the cosym coordinates, their dimension and OPTICS follow
// (docs/correlation.md).

#ifndef MXI_CORRELATION_HH
#define MXI_CORRELATION_HH

#include <cstddef>
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

} // namespace mxi

#endif // MXI_CORRELATION_HH
