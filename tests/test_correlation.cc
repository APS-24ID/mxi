// The correlation matrix of mxi_correlation (src/correlation.hh): DIALS's
// weighted_cchalf, its effective sample size, the count weights, and two
// planted groups of data sets told apart.

#include <cmath>
#include <random>
#include <vector>

#include "../src/correlation.hh"
#include "check.hh"

namespace mxi {

namespace {

// n reflections, (h, 0, 0) for h from 1, the given intensities and sigmas.
MergedSet set_of(const std::vector<double> &intensity,
                 const std::vector<double> &sigma) {
  MergedSet s;
  for (std::size_t k = 0; k < intensity.size(); ++k)
    s.hkl.push_back(Miller{static_cast<int>(k) + 1, 0, 0});
  s.intensity = intensity;
  s.sigma = sigma;
  return s;
}

} // namespace

TEST(the_pair_correlation_is_dials_weighted_cchalf) {
  const std::vector<double> x = {10, 20, 5, 40, 33, 12, 7},
                            sx = {1, 2, 1, 4, 3, 1, 1},
                            y = {11, 18, 6, 42, 30, 15, 6},
                            sy = {2, 1, 1, 3, 3, 2, 1};
  const PairCorrelation c = weighted_correlation(set_of(x, sx), set_of(y, sy));
  // The same, written out: w = 1 / (sx^2 + sy^2), normalised.
  double sw = 0.0;
  std::vector<double> w;
  for (std::size_t k = 0; k < x.size(); ++k) {
    w.push_back(1.0 / (sx[k] * sx[k] + sy[k] * sy[k]));
    sw += w.back();
  }
  double xb = 0, yb = 0, s2 = 0;
  for (std::size_t k = 0; k < x.size(); ++k) {
    xb += x[k] * w[k] / sw;
    yb += y[k] * w[k] / sw;
    s2 += (w[k] / sw) * (w[k] / sw);
  }
  double cxy = 0, cxx = 0, cyy = 0;
  for (std::size_t k = 0; k < x.size(); ++k) {
    cxy += (x[k] - xb) * (y[k] - yb) * w[k] / sw;
    cxx += (x[k] - xb) * (x[k] - xb) * w[k] / sw;
    cyy += (y[k] - yb) * (y[k] - yb) * w[k] / sw;
  }
  check::is_true(c.valid, "valid");
  check::close(c.r, cxy / std::sqrt(cxx * cyy), 1e-12, "r");
  check::close(c.n_eff, 1.0 / s2, 1e-12, "Kish's effective sample size");
  check::equal(static_cast<long long>(c.common), 7LL, "all seven in common");
}

TEST(the_matrix_weights_sum_to_the_pairs_that_have_one) {
  std::mt19937 rng(3);
  std::normal_distribution<double> noise(0.0, 1.0);
  std::vector<MergedSet> sets;
  for (int d = 0; d < 5; ++d) {
    std::vector<double> I, s;
    for (int k = 0; k < 200; ++k) {
      I.push_back(100.0 + 50.0 * std::sin(k) + noise(rng) * 5.0);
      s.push_back(5.0);
    }
    sets.push_back(set_of(I, s));
  }
  const CorrelationMatrix m = correlation_matrix(sets);
  double total = 0.0;
  std::size_t nonzero = 0;
  for (std::size_t i = 0; i < m.n; ++i) {
    check::is_true(m.w[i * m.n + i] == 0.0, "nothing on the diagonal");
    for (std::size_t j = 0; j < m.n; ++j) {
      check::is_true(m.r[i * m.n + j] == m.r[j * m.n + i], "symmetric");
      if (m.w[i * m.n + j] > 0.0) {
        total += m.w[i * m.n + j];
        ++nonzero;
      }
    }
  }
  check::equal(static_cast<long long>(nonzero), 20LL, "every pair has one");
  check::close(total, static_cast<double>(nonzero), 1e-9, "rescaled");
}

TEST(two_planted_groups_correlate_within_more_than_between) {
  // Common intensities, each group perturbed its own way, each set noisy.
  std::mt19937 rng(5);
  std::normal_distribution<double> noise(0.0, 1.0);
  std::vector<double> truth, shift_a, shift_b;
  for (int k = 0; k < 400; ++k) {
    truth.push_back(200.0 * std::exp(noise(rng) * 0.8));
    shift_a.push_back(noise(rng) * 40.0);
    shift_b.push_back(noise(rng) * 40.0);
  }
  std::vector<MergedSet> sets;
  for (int d = 0; d < 8; ++d) {
    const std::vector<double> &shift = d < 4 ? shift_a : shift_b;
    std::vector<double> I, s;
    for (int k = 0; k < 400; ++k) {
      I.push_back(truth[k] + shift[k] + noise(rng) * 10.0);
      s.push_back(10.0);
    }
    sets.push_back(set_of(I, s));
  }
  const CorrelationMatrix m = correlation_matrix(sets);
  double within = 0.0, between = 0.0;
  int nw = 0, nb = 0;
  for (std::size_t i = 0; i < 8; ++i)
    for (std::size_t j = i + 1; j < 8; ++j) {
      if ((i < 4) == (j < 4)) {
        within += m.r[i * 8 + j];
        ++nw;
      } else {
        between += m.r[i * 8 + j];
        ++nb;
      }
    }
  check::is_true(within / nw > between / nb + 0.02,
                 "within " + std::to_string(within / nw) + ", between " +
                     std::to_string(between / nb));
}

namespace {

// Three groups of data sets: common intensities, each group perturbed its own
// way, every set noisy -- as bovine, porcine and human insulin.
CorrelationMatrix three_groups(std::size_t each, unsigned seed) {
  std::mt19937 rng(seed);
  std::normal_distribution<double> noise(0.0, 1.0);
  std::vector<double> truth;
  std::vector<std::vector<double>> shifts(3);
  for (int k = 0; k < 400; ++k) {
    truth.push_back(200.0 * std::exp(noise(rng) * 0.8));
    for (auto &s : shifts)
      s.push_back(noise(rng) * 40.0);
  }
  std::vector<MergedSet> sets;
  for (std::size_t d = 0; d < 3 * each; ++d) {
    std::vector<double> I, s;
    for (int k = 0; k < 400; ++k) {
      I.push_back(truth[k] + shifts[d / each][k] + noise(rng) * 10.0);
      s.push_back(10.0);
    }
    sets.push_back(set_of(I, s));
  }
  return correlation_matrix(sets);
}

} // namespace

TEST(the_cosym_start_is_numpys_at_seed_230) {
  // np.random.seed(230); np.random.rand(5)
  const double numpy[5] = {0.6309999826057662, 0.8823043197860562,
                           0.258551660640161, 0.3141223801202194,
                           0.9358842864512585};
  const CorrelationMatrix m = three_groups(2, 1);
  const std::vector<double> x = minimise_cosym(m, 1, 0, 0, 230);
  for (int k = 0; k < 5; ++k)
    check::is_true(x[k] == numpy[k], "draw " + std::to_string(k));
}

TEST(the_cosym_gradient_is_the_functionals) {
  const CorrelationMatrix m = three_groups(3, 2);
  const std::size_t dim = 3;
  std::vector<double> x = minimise_cosym(m, dim, 0, 0, 7), grad;
  // The gradient through a step of the minimiser's own: compare with finite
  // differences of cosym_functional.
  const double f0 = cosym_functional(m, x, dim);
  (void)f0;
  for (std::size_t k : {0u, 5u, 11u, 20u}) {
    std::vector<double> up = x, down = x;
    up[k] += 1e-6;
    down[k] -= 1e-6;
    const double numeric =
        (cosym_functional(m, up, dim) - cosym_functional(m, down, dim)) / 2e-6;
    // the analytic gradient, -2 X (W.R - W.G), written out for coordinate k
    const std::size_t n = m.n, c = k / n, i = k % n;
    double analytic = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      double g = 0.0;
      for (std::size_t a = 0; a < dim; ++a)
        g += x[a * n + i] * x[a * n + j];
      analytic += x[c * n + j] * m.w[j * n + i] * (m.r[j * n + i] - g);
    }
    analytic *= -2.0;
    check::close(analytic, numeric, 1e-5 * std::max(1.0, std::fabs(numeric)),
                 "coordinate " + std::to_string(k));
  }
}

TEST(the_elbow_is_where_the_curve_turns) {
  const std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8};
  const std::vector<double> y = {1.0,   0.5,  0.05,  0.04,
                                 0.035, 0.03, 0.028, 0.025};
  check::equal(static_cast<long long>(elbow(x, y)), 3LL, "the elbow at 3");
}

TEST(three_planted_groups_point_three_ways) {
  const CorrelationMatrix m = three_groups(6, 9);
  CosymOptions options;
  const CosymCoordinates c = cosym_coordinates(m, options);
  check::is_true(c.dim >= 2, "dimensions " + std::to_string(c.dim));
  // Cosines between data sets' coordinates: within a group near 1.
  const std::size_t n = c.n;
  const auto cosine = [&](std::size_t i, std::size_t j) {
    double ab = 0, aa = 0, bb = 0;
    for (std::size_t k = 0; k < c.dim; ++k) {
      ab += c.x[k * n + i] * c.x[k * n + j];
      aa += c.x[k * n + i] * c.x[k * n + i];
      bb += c.x[k * n + j] * c.x[k * n + j];
    }
    return ab / std::sqrt(aa * bb);
  };
  double within = 1.0, between = -1.0;
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j) {
      if (i / 6 == j / 6)
        within = std::min(within, cosine(i, j));
      else
        between = std::max(between, cosine(i, j));
    }
  check::is_true(within > between, "least within " + std::to_string(within) +
                                       ", most between " +
                                       std::to_string(between));
  double ratio = 0.0;
  for (double r : c.explained_variance_ratio)
    ratio += r;
  check::close(ratio, 1.0, 1e-9, "the explained variance ratios sum to one");
}

} // namespace mxi
