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

} // namespace mxi

#endif // MXI_CORRELATION_HH
