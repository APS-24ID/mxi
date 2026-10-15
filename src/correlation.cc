// Intensity-based clustering of isomorphous data sets: see correlation.hh.

#include "correlation.hh"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <tuple>

#include "parallel.hh"
#include "scale.hh"

namespace mxi {

std::vector<MergedSet> merged_by_experiment(const ExperimentList &experiments,
                                            const Table &reflections,
                                            const SpaceGroup &group,
                                            double min_partiality,
                                            std::vector<std::size_t> *which) {
  std::vector<std::size_t> lattices;
  std::vector<long> slot(experiments.size(), -1);
  for (std::size_t e = 0; e < experiments.size(); ++e)
    if (experiments[e].crystal) {
      slot[e] = static_cast<long>(lattices.size());
      lattices.push_back(e);
    }
  if (which)
    *which = lattices;
  std::vector<MergedSet> out(lattices.size());
  for (const char *c : {"id", "miller_index", "flags", "intensity.scale.value",
                        "intensity.scale.variance"})
    if (!reflections.has(c))
      return out;
  const Column &id = reflections.at("id"),
               &hkl = reflections.at("miller_index"),
               &flags = reflections.at("flags"),
               &I = reflections.at("intensity.scale.value"),
               &V = reflections.at("intensity.scale.variance");
  const Column *partiality =
      reflections.has("partiality") ? &reflections.at("partiality") : nullptr;
  // Each set's observations by symmetry-unique index, in the table's order.
  struct Obs {
    double value, sigma;
  };
  std::vector<std::map<Miller, std::vector<Obs>>> groups(lattices.size());
  for (std::size_t r = 0; r < reflections.nrows; ++r) {
    const std::int64_t e = id.integer(r);
    if (e < 0 || static_cast<std::size_t>(e) >= slot.size() || slot[e] < 0)
      continue;
    const std::int64_t f = flags.integer(r);
    if ((f & flag::kScaled) == 0 || (f & flag::kOutlierInScaling) != 0)
      continue;
    if (partiality && !(partiality->reals[r] >= min_partiality))
      continue;
    const double v = V.reals[r];
    if (!(v > 0.0) || !std::isfinite(I.reals[r]))
      continue;
    const Miller m{static_cast<int>(hkl.ints[r * 3]),
                   static_cast<int>(hkl.ints[r * 3 + 1]),
                   static_cast<int>(hkl.ints[r * 3 + 2])};
    groups[static_cast<std::size_t>(slot[e])][group.unique(m)].push_back(
        {I.reals[r], std::sqrt(v)});
  }
  // Merged as cctbx's merge_equivalents does by default (gaussian, with the
  // internal variance): one observation as it is; several, the inverse-
  // variance weighted mean, its sigma sqrt(max(V_w / n, 1 / sum w)), V_w the
  // weighted variance as GSL's: sum w / ((sum w)^2 - sum w^2) sum w (x -
  // mean)^2.
  for (std::size_t k = 0; k < lattices.size(); ++k) {
    MergedSet &set = out[k];
    for (const auto &[m, obs] : groups[k]) {
      double mean = 0.0, sigma = 0.0;
      if (obs.size() == 1) {
        mean = obs[0].value;
        sigma = obs[0].sigma;
      } else {
        double sw = 0.0, swx = 0.0, sw2 = 0.0;
        for (const Obs &o : obs) {
          const double w = 1.0 / (o.sigma * o.sigma);
          sw += w;
          swx += w * o.value;
          sw2 += w * w;
        }
        mean = swx / sw;
        double swd2 = 0.0;
        for (const Obs &o : obs) {
          const double w = 1.0 / (o.sigma * o.sigma);
          swd2 += w * (o.value - mean) * (o.value - mean);
        }
        const double wvar = sw * sw > sw2 ? sw / (sw * sw - sw2) * swd2 : 0.0;
        sigma = std::sqrt(
            std::max(wvar / static_cast<double>(obs.size()), 1.0 / sw));
      }
      set.hkl.push_back(m);
      set.intensity.push_back(mean);
      set.sigma.push_back(sigma);
    }
  }
  return out;
}

PairCorrelation weighted_correlation(const MergedSet &a, const MergedSet &b,
                                     std::size_t min_pairs) {
  PairCorrelation out;
  // The reflections both have: the indices are ascending in each.
  std::vector<double> x, y, w;
  std::size_t i = 0, j = 0;
  while (i < a.hkl.size() && j < b.hkl.size()) {
    if (a.hkl[i] < b.hkl[j]) {
      ++i;
    } else if (b.hkl[j] < a.hkl[i]) {
      ++j;
    } else {
      x.push_back(a.intensity[i]);
      y.push_back(b.intensity[j]);
      w.push_back(1.0 / (a.sigma[i] * a.sigma[i] + b.sigma[j] * b.sigma[j]));
      ++i;
      ++j;
    }
  }
  out.common = x.size();
  if (x.size() < std::max<std::size_t>(min_pairs, 2))
    return out;
  double sw = 0.0;
  for (double v : w)
    sw += v;
  double xbar = 0.0, ybar = 0.0, sum_w2 = 0.0;
  for (std::size_t k = 0; k < x.size(); ++k) {
    const double nw = w[k] / sw;
    xbar += x[k] * nw;
    ybar += y[k] * nw;
    sum_w2 += nw * nw;
  }
  double sxy = 0.0, sx = 0.0, sy = 0.0;
  for (std::size_t k = 0; k < x.size(); ++k) {
    const double nw = w[k] / sw, dx = x[k] - xbar, dy = y[k] - ybar;
    sxy += dx * dy * nw;
    sx += dx * dx * nw;
    sy += dy * dy * nw;
  }
  if (!(sx > 0.0) || !(sy > 0.0))
    return out;
  out.valid = true;
  out.r = sxy / std::sqrt(sx * sy);
  out.n_eff = 1.0 / sum_w2;
  return out;
}

CorrelationMatrix correlation_matrix(const std::vector<MergedSet> &sets,
                                     std::size_t min_pairs) {
  CorrelationMatrix out;
  const std::size_t n = sets.size();
  out.n = n;
  out.r.assign(n * n, 0.0);
  out.w.assign(n * n, 0.0);
  // A row a task: each pair once, i < j, both halves written.
  for_each_index(n, [&](std::size_t i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      const PairCorrelation c =
          weighted_correlation(sets[i], sets[j], min_pairs);
      if (!c.valid)
        continue;
      out.r[i * n + j] = out.r[j * n + i] = c.r;
      out.w[i * n + j] = out.w[j * n + i] = c.n_eff;
    }
  });
  // DIALS's count weights: n_eff, scaled so that they sum to the number of
  // entries that have one.
  double total = 0.0;
  std::size_t nonzero = 0;
  for (double v : out.w)
    if (v > 0.0) {
      total += v;
      ++nonzero;
    }
  if (total > 0.0) {
    const double scale = static_cast<double>(nonzero) / total;
    for (double &v : out.w)
      v *= scale;
  }
  return out;
}

namespace {

// numpy's np.random.seed(seed); np.random.rand(count): its legacy seeding is
// the Mersenne twister's init_genrand, and each double is two draws, 53 bits:
// (a >> 5, b >> 6) -> (a 2^26 + b) / 2^53.
std::vector<double> numpy_rand(unsigned seed, std::size_t count) {
  std::mt19937 twister(seed);
  std::vector<double> out(count);
  for (double &v : out) {
    const double a = static_cast<double>(twister() >> 5);
    const double b = static_cast<double>(twister() >> 6);
    v = (a * 67108864.0 + b) / 9007199254740992.0;
  }
  return out;
}

// The Gram matrix x^T x, n x n, of coordinates dimension-major.
std::vector<double> gram(const std::vector<double> &x, std::size_t n,
                         std::size_t dim) {
  std::vector<double> g(n * n, 0.0);
  for_each_index(n, [&](std::size_t i) {
    for (std::size_t j = 0; j < n; ++j) {
      double v = 0.0;
      for (std::size_t k = 0; k < dim; ++k)
        v += x[k * n + i] * x[k * n + j];
      g[i * n + j] = v;
    }
  });
  return g;
}

// The functional and its gradient.
double value_and_gradient(const CorrelationMatrix &m,
                          const std::vector<double> &x, std::size_t dim,
                          std::vector<double> *grad) {
  const std::size_t n = m.n;
  const std::vector<double> g = gram(x, n, dim);
  double f = 0.0;
  std::vector<double> d(n * n);
  for (std::size_t k = 0; k < n * n; ++k) {
    const double e = m.r[k] - g[k];
    f += m.w[k] * e * e;
    d[k] = m.w[k] * e;
  }
  if (grad) {
    grad->assign(dim * n, 0.0);
    for_each_index(dim, [&](std::size_t k) {
      for (std::size_t i = 0; i < n; ++i) {
        double v = 0.0;
        for (std::size_t j = 0; j < n; ++j)
          v += x[k * n + j] * d[j * n + i];
        (*grad)[k * n + i] = -2.0 * v;
      }
    });
  }
  return 0.5 * f;
}

// A symmetric matrix's eigenvalues and eigenvectors (columns of v), by Jacobi
// rotations: for the small matrices here, d x d, d at most fifty.
void symmetric_eigen(std::vector<double> a, std::size_t d,
                     std::vector<double> *values, std::vector<double> *v) {
  v->assign(d * d, 0.0);
  for (std::size_t i = 0; i < d; ++i)
    (*v)[i * d + i] = 1.0;
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0;
    for (std::size_t p = 0; p < d; ++p)
      for (std::size_t q = p + 1; q < d; ++q)
        off += a[p * d + q] * a[p * d + q];
    if (off < 1e-30)
      break;
    for (std::size_t p = 0; p < d; ++p)
      for (std::size_t q = p + 1; q < d; ++q) {
        const double apq = a[p * d + q];
        if (std::fabs(apq) < 1e-300)
          continue;
        const double theta = (a[q * d + q] - a[p * d + p]) / (2.0 * apq);
        const double t = (theta >= 0 ? 1.0 : -1.0) /
                         (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0), sn = t * c;
        for (std::size_t k = 0; k < d; ++k) {
          const double akp = a[k * d + p], akq = a[k * d + q];
          a[k * d + p] = c * akp - sn * akq;
          a[k * d + q] = sn * akp + c * akq;
        }
        for (std::size_t k = 0; k < d; ++k) {
          const double apk = a[p * d + k], aqk = a[q * d + k];
          a[p * d + k] = c * apk - sn * aqk;
          a[q * d + k] = sn * apk + c * aqk;
        }
        for (std::size_t k = 0; k < d; ++k) {
          const double vkp = (*v)[k * d + p], vkq = (*v)[k * d + q];
          (*v)[k * d + p] = c * vkp - sn * vkq;
          (*v)[k * d + q] = sn * vkp + c * vkq;
        }
      }
  }
  values->resize(d);
  for (std::size_t i = 0; i < d; ++i)
    (*values)[i] = a[i * d + i];
}

// Eigen-pairs sorted by eigenvalue, largest first.
std::vector<std::size_t> descending(const std::vector<double> &values) {
  std::vector<std::size_t> order(values.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(
      order.begin(), order.end(),
      [&](std::size_t a, std::size_t b) { return values[a] > values[b]; });
  return order;
}

} // namespace

double cosym_functional(const CorrelationMatrix &m,
                        const std::vector<double> &x, std::size_t dim) {
  return value_and_gradient(m, x, dim, nullptr);
}

std::vector<double> minimise_cosym(const CorrelationMatrix &m, std::size_t dim,
                                   int max_iterations, int max_calls,
                                   unsigned seed, double *functional) {
  const std::size_t n = m.n, size = dim * n;
  std::vector<double> x = numpy_rand(seed, size), g, x_new, g_new;
  double f = value_and_gradient(m, x, dim, &g);
  int calls = 1;
  // L-BFGS, ten pairs remembered, a backtracking line search to Armijo's
  // condition; stopped at scitbx's test, |g| below 1e-5 of max(1, |x|), or at
  // the iterations or calls allowed.
  const std::size_t memory = 10;
  std::vector<std::vector<double>> s_hist, y_hist;
  std::vector<double> rho;
  const auto dot = [](const std::vector<double> &a,
                      const std::vector<double> &b) {
    double v = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k)
      v += a[k] * b[k];
    return v;
  };
  for (int it = 0; it < max_iterations; ++it) {
    const double gnorm = std::sqrt(dot(g, g));
    if (gnorm <= 1e-5 * std::max(1.0, std::sqrt(dot(x, x))))
      break;
    // The two-loop recursion: the direction -H g.
    std::vector<double> q = g;
    std::vector<double> alpha(s_hist.size());
    for (std::size_t k = s_hist.size(); k-- > 0;) {
      alpha[k] = rho[k] * dot(s_hist[k], q);
      for (std::size_t j = 0; j < size; ++j)
        q[j] -= alpha[k] * y_hist[k][j];
    }
    // The starting inverse Hessian: 1 / the curvature, coordinate by
    // coordinate, every iteration, as DIALS's minimiser (scitbx's lbfgs,
    // diag_mode "always"): c_ki = 2 sum_j x_kj^2 w_ji (target.curvatures). A
    // single scale instead got almost nowhere in the scan's 20 calls a
    // dimension -- its curve not even falling, the elbow at 16 of 18 where 3
    // fitted to 7e-7.
    for (std::size_t k = 0; k < dim; ++k)
      for (std::size_t i = 0; i < n; ++i) {
        double c = 0.0;
        for (std::size_t j = 0; j < n; ++j)
          c += x[k * n + j] * x[k * n + j] * m.w[j * n + i];
        c *= 2.0;
        q[k * n + i] *= c > 1e-300 ? 1.0 / c : 1.0 / std::max(gnorm, 1e-300);
      }
    for (std::size_t k = 0; k < s_hist.size(); ++k) {
      const double beta = rho[k] * dot(y_hist[k], q);
      for (std::size_t j = 0; j < size; ++j)
        q[j] += s_hist[k][j] * (alpha[k] - beta);
    }
    double slope = -dot(g, q);
    if (!(slope < 0.0)) { // not a descent direction: steepest descent
      for (std::size_t j = 0; j < size; ++j)
        q[j] = g[j] / std::max(gnorm, 1e-300);
      slope = -dot(g, q);
      s_hist.clear();
      y_hist.clear();
      rho.clear();
    }
    double step = 1.0, f_new = f;
    bool moved = false;
    for (int tries = 0; tries < 30; ++tries) {
      if (max_calls > 0 && calls >= max_calls)
        break;
      x_new = x;
      for (std::size_t j = 0; j < size; ++j)
        x_new[j] -= step * q[j];
      f_new = value_and_gradient(m, x_new, dim, &g_new);
      ++calls;
      if (f_new <= f + 1e-4 * step * slope) {
        moved = true;
        break;
      }
      step *= 0.5;
    }
    if (!moved)
      break;
    std::vector<double> sk(size), yk(size);
    for (std::size_t j = 0; j < size; ++j) {
      sk[j] = x_new[j] - x[j];
      yk[j] = g_new[j] - g[j];
    }
    const double sy = dot(sk, yk);
    if (sy > 1e-300) {
      s_hist.push_back(std::move(sk));
      y_hist.push_back(std::move(yk));
      rho.push_back(1.0 / sy);
      if (s_hist.size() > memory) {
        s_hist.erase(s_hist.begin());
        y_hist.erase(y_hist.begin());
        rho.erase(rho.begin());
      }
    }
    x.swap(x_new);
    g.swap(g_new);
    f = f_new;
    if (max_calls > 0 && calls >= max_calls)
      break;
  }
  if (functional)
    *functional = f;
  return x;
}

double dimension_score(const CorrelationMatrix &m, const std::vector<double> &x,
                       std::size_t dim, bool outlier_rejection) {
  if (!outlier_rejection)
    return cosym_functional(m, x, dim);
  const std::size_t n = m.n;
  const std::vector<double> g = gram(x, n, dim);
  std::vector<double> e(n * n);
  for (std::size_t k = 0; k < n * n; ++k) {
    const double d = m.r[k] - g[k];
    e[k] = m.w[k] * d * d;
  }
  // numpy's quantile, linear between order statistics.
  std::vector<double> sorted = e;
  std::sort(sorted.begin(), sorted.end());
  const auto quantile = [&](double q) {
    const double at = q * static_cast<double>(sorted.size() - 1);
    const auto lo = static_cast<std::size_t>(std::floor(at));
    const std::size_t hi = std::min(lo + 1, sorted.size() - 1);
    return sorted[lo] +
           (at - static_cast<double>(lo)) * (sorted[hi] - sorted[lo]);
  };
  const double q1 = quantile(0.25), q2 = quantile(0.5), q3 = quantile(0.75);
  double sum = 0.0;
  for (double v : e)
    if (v < q2 + (q3 - q1))
      sum += v;
  return 0.5 * sum;
}

std::size_t elbow(const std::vector<double> &x, const std::vector<double> &y) {
  const std::size_t n = x.size();
  if (n < 2)
    return n == 1 ? static_cast<std::size_t>(x[0]) : 0;
  std::size_t p_m = 0;
  double least = HUGE_VAL;
  for (std::size_t i = 0; i + 1 < n; ++i) {
    const double slope = (y[n - 1] - y[i]) / (x[n - 1] - x[i]);
    if (slope < least) {
      least = slope;
      p_m = i;
    }
  }
  const double x1 = x[p_m], y1 = y[p_m], x2 = x[n - 1], y2 = y[n - 1];
  double vx = y2 - y1, vy = -(x2 - x1);
  const double norm = std::sqrt(vx * vx + vy * vy);
  vx /= norm;
  vy /= norm;
  std::size_t best = p_m;
  double gap = -1.0;
  for (std::size_t i = p_m; i < n; ++i) {
    const double g = std::fabs(vx * (x1 - x[i]) + vy * (y1 - y[i]));
    if (g > gap) {
      gap = g;
      best = i;
    }
  }
  return static_cast<std::size_t>(x[best]);
}

CosymCoordinates cosym_coordinates(const CorrelationMatrix &m,
                                   const CosymOptions &options) {
  CosymCoordinates out;
  const std::size_t n = m.n;
  out.n = n;
  if (n < 2)
    return out;
  const std::size_t tested = std::min(n, options.maximum_dimensions);
  std::size_t dim = options.dimensions;
  if (dim == 0) {
    if (options.by_projection) {
      // One minimisation at the most dimensions, rotated to its principal
      // axes (from the small Gram matrix x x^T), each smaller dimension's
      // functional from the leading axes (DIALS #3248).
      const std::vector<double> x =
          minimise_cosym(m, tested, options.max_iterations,
                         options.scan_max_calls, options.seed);
      std::vector<double> c(tested * tested, 0.0);
      for (std::size_t a = 0; a < tested; ++a)
        for (std::size_t b = 0; b < tested; ++b)
          for (std::size_t i = 0; i < n; ++i)
            c[a * tested + b] += x[a * n + i] * x[b * n + i];
      std::vector<double> values, vectors;
      symmetric_eigen(c, tested, &values, &vectors);
      const std::vector<std::size_t> order = descending(values);
      std::vector<double> rotated(tested * n, 0.0);
      for (std::size_t a = 0; a < tested; ++a)
        for (std::size_t i = 0; i < n; ++i) {
          double v = 0.0;
          for (std::size_t b = 0; b < tested; ++b)
            v += vectors[b * tested + order[a]] * x[b * n + i];
          rotated[a * n + i] = v;
        }
      for (std::size_t d = 1; d <= tested; ++d) {
        const std::vector<double> lead(
            rotated.begin(), rotated.begin() + static_cast<long>(d * n));
        out.dimensions.push_back(static_cast<double>(d));
        out.functional.push_back(
            dimension_score(m, lead, d, options.outlier_rejection));
      }
    } else {
      // A minimisation at each dimension, from the same start, as
      // dials.correlation_matrix.
      for (std::size_t d = 1; d <= tested; ++d) {
        const std::vector<double> x = minimise_cosym(
            m, d, options.max_iterations, options.scan_max_calls, options.seed);
        out.dimensions.push_back(static_cast<double>(d));
        out.functional.push_back(
            dimension_score(m, x, d, options.outlier_rejection));
      }
    }
    dim = std::max<std::size_t>(elbow(out.dimensions, out.functional), 2);
  }
  dim = std::min(dim, n);
  out.dim = dim;
  out.x = minimise_cosym(m, dim, options.max_iterations, 0, options.seed,
                         &out.final_functional);
  // Principal components, as sklearn's PCA: centred, the covariance over
  // n - 1, its eigenvectors largest first; each data set's coordinates on them.
  std::vector<double> mean(dim, 0.0);
  for (std::size_t k = 0; k < dim; ++k) {
    for (std::size_t i = 0; i < n; ++i)
      mean[k] += out.x[k * n + i];
    mean[k] /= static_cast<double>(n);
  }
  std::vector<double> cov(dim * dim, 0.0);
  for (std::size_t a = 0; a < dim; ++a)
    for (std::size_t b = 0; b < dim; ++b) {
      double v = 0.0;
      for (std::size_t i = 0; i < n; ++i)
        v += (out.x[a * n + i] - mean[a]) * (out.x[b * n + i] - mean[b]);
      cov[a * dim + b] = v / static_cast<double>(n - 1);
    }
  std::vector<double> values, vectors;
  symmetric_eigen(cov, dim, &values, &vectors);
  const std::vector<std::size_t> order = descending(values);
  double total = 0.0;
  for (double v : values)
    total += std::max(v, 0.0);
  out.components.assign(n * dim, 0.0);
  for (std::size_t c = 0; c < dim; ++c) {
    const double ev = std::max(values[order[c]], 0.0);
    out.explained_variance.push_back(ev);
    out.explained_variance_ratio.push_back(total > 0.0 ? ev / total : 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      double v = 0.0;
      for (std::size_t k = 0; k < dim; ++k)
        v += (out.x[k * n + i] - mean[k]) * vectors[k * dim + order[c]];
      out.components[i * dim + c] = v;
    }
  }
  return out;
}

namespace {

// np.around(x, 15): x times 10^15, rounded half to even, divided back.
double around15(double x) {
  if (!std::isfinite(x))
    return x;
  return std::nearbyint(x * 1e15) / 1e15;
}

double distance(const std::vector<double> &p, std::size_t d, std::size_t i,
                std::size_t j) {
  double s = 0.0;
  for (std::size_t k = 0; k < d; ++k) {
    const double e = p[i * d + k] - p[j * d + k];
    s += e * e;
  }
  return std::sqrt(s);
}

// scikit-learn's _extend_region.
std::size_t extend_region(const std::vector<bool> &steep,
                          const std::vector<bool> &xward, std::size_t start,
                          std::size_t min_samples) {
  const std::size_t n = steep.size();
  std::size_t non_xward = 0, index = start, end = start;
  while (index < n) {
    if (steep[index]) {
      non_xward = 0;
      end = index;
    } else if (!xward[index]) {
      ++non_xward;
      if (non_xward > min_samples)
        break;
    } else {
      return end;
    }
    ++index;
  }
  return end;
}

struct SteepDown {
  std::size_t start, end;
  double mib;
};

// scikit-learn's _update_filter_sdas.
std::vector<SteepDown> update_filter_sdas(const std::vector<SteepDown> &sdas,
                                          double mib, double xi_complement,
                                          const std::vector<double> &plot) {
  if (std::isinf(mib))
    return {};
  std::vector<SteepDown> out;
  for (SteepDown sda : sdas)
    if (mib <= plot[sda.start] * xi_complement) {
      sda.mib = std::max(sda.mib, mib);
      out.push_back(sda);
    }
  return out;
}

// scikit-learn's _xi_cluster, with predecessor correction.
std::vector<std::pair<std::size_t, std::size_t>>
xi_cluster(std::vector<double> plot, const std::vector<long> &predecessor_plot,
           const std::vector<std::size_t> &ordering, double xi,
           std::size_t min_samples, std::size_t min_cluster_size) {
  plot.push_back(HUGE_VAL);
  const double xi_complement = 1.0 - xi;
  const std::size_t n = plot.size() - 1;
  std::vector<bool> steep_up(n), steep_down(n), down(n), up(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double ratio = plot[i] / plot[i + 1]; // inf/inf NaN: all false
    steep_up[i] = ratio <= xi_complement;
    steep_down[i] = ratio >= 1.0 / xi_complement;
    down[i] = ratio > 1.0;
    up[i] = ratio < 1.0;
  }
  std::vector<SteepDown> sdas;
  std::vector<std::pair<std::size_t, std::size_t>> clusters;
  std::size_t index = 0;
  double mib = 0.0;
  const auto correct_predecessor = [&](std::size_t s, std::size_t e, bool *ok) {
    while (s < e) {
      if (plot[s] > plot[e]) {
        *ok = true;
        return std::make_pair(s, e);
      }
      const long p_e = predecessor_plot[e];
      for (std::size_t i = s; i < e; ++i)
        if (p_e == static_cast<long>(ordering[i])) {
          *ok = true;
          return std::make_pair(s, e);
        }
      --e;
    }
    *ok = false;
    return std::make_pair(s, e);
  };
  for (std::size_t steep_index = 0; steep_index < n; ++steep_index) {
    if (!(steep_up[steep_index] || steep_down[steep_index]))
      continue;
    if (steep_index < index)
      continue;
    for (std::size_t k = index; k <= steep_index; ++k)
      mib = std::max(mib, plot[k]);
    if (steep_down[steep_index]) {
      sdas = update_filter_sdas(sdas, mib, xi_complement, plot);
      const std::size_t d_start = steep_index;
      const std::size_t d_end =
          extend_region(steep_down, up, d_start, min_samples);
      sdas.push_back({d_start, d_end, 0.0});
      index = d_end + 1;
      mib = plot[index];
    } else {
      sdas = update_filter_sdas(sdas, mib, xi_complement, plot);
      const std::size_t u_start = steep_index;
      const std::size_t u_end =
          extend_region(steep_up, down, u_start, min_samples);
      index = u_end + 1;
      mib = plot[index];
      std::vector<std::pair<std::size_t, std::size_t>> u_clusters;
      for (const SteepDown &sda : sdas) {
        std::size_t c_start = sda.start, c_end = u_end;
        if (plot[c_end + 1] * xi_complement < sda.mib)
          continue;
        const double d_max = plot[sda.start];
        if (d_max * xi_complement >= plot[c_end + 1]) {
          while (plot[c_start + 1] > plot[c_end + 1] && c_start < sda.end)
            ++c_start;
        } else if (plot[c_end + 1] * xi_complement >= d_max) {
          while (plot[c_end - 1] > d_max && c_end > u_start)
            --c_end;
        }
        bool ok = true;
        std::tie(c_start, c_end) = correct_predecessor(c_start, c_end, &ok);
        if (!ok)
          continue;
        if (c_end - c_start + 1 < min_cluster_size)
          continue;
        if (c_start > sda.end)
          continue;
        if (c_end < u_start)
          continue;
        u_clusters.emplace_back(c_start, c_end);
      }
      std::reverse(u_clusters.begin(), u_clusters.end());
      clusters.insert(clusters.end(), u_clusters.begin(), u_clusters.end());
    }
  }
  return clusters;
}

} // namespace

OpticsResult optics(const std::vector<double> &points, std::size_t n,
                    std::size_t d, std::size_t min_samples, double max_eps,
                    double xi) {
  OpticsResult out;
  out.reachability.assign(n, HUGE_VAL);
  out.predecessor.assign(n, -1);
  out.core_distance.assign(n, HUGE_VAL);
  out.labels.assign(n, -1);
  if (n == 0)
    return out;
  min_samples = std::max<std::size_t>(std::min(min_samples, n), 1);
  for_each_index(n, [&](std::size_t i) {
    std::vector<double> dist(n);
    for (std::size_t j = 0; j < n; ++j)
      dist[j] = distance(points, d, i, j);
    std::nth_element(dist.begin(),
                     dist.begin() + static_cast<long>(min_samples - 1),
                     dist.end());
    const double core = dist[min_samples - 1];
    out.core_distance[i] = core > max_eps ? HUGE_VAL : around15(core);
  });
  std::vector<bool> processed(n, false);
  for (std::size_t o = 0; o < n; ++o) {
    std::size_t point = n;
    for (std::size_t i = 0; i < n; ++i)
      if (!processed[i] &&
          (point == n || out.reachability[i] < out.reachability[point]))
        point = i;
    processed[point] = true;
    out.ordering.push_back(point);
    if (std::isinf(out.core_distance[point]))
      continue;
    for (std::size_t j = 0; j < n; ++j) {
      if (processed[j])
        continue;
      const double dist = distance(points, d, point, j);
      if (dist > max_eps)
        continue;
      const double r = around15(std::max(dist, out.core_distance[point]));
      if (r < out.reachability[j]) {
        out.reachability[j] = r;
        out.predecessor[j] = static_cast<long>(point);
      }
    }
  }
  std::vector<double> plot(n);
  std::vector<long> predecessor_plot(n);
  for (std::size_t k = 0; k < n; ++k) {
    plot[k] = out.reachability[out.ordering[k]];
    predecessor_plot[k] = out.predecessor[out.ordering[k]];
  }
  out.clusters = xi_cluster(plot, predecessor_plot, out.ordering, xi,
                            min_samples, min_samples);
  // scikit-learn's _extract_xi_labels: the clusters in order, each labelled
  // where none of it is yet.
  std::vector<int> in_order(n, -1);
  int label = 0;
  for (const auto &[a, b] : out.clusters) {
    bool free = true;
    for (std::size_t k = a; k <= b; ++k)
      free = free && in_order[k] == -1;
    if (free) {
      for (std::size_t k = a; k <= b; ++k)
        in_order[k] = label;
      ++label;
    }
  }
  for (std::size_t k = 0; k < n; ++k)
    out.labels[out.ordering[k]] = in_order[k];
  return out;
}

double silhouette(const std::vector<double> &points, std::size_t d,
                  const std::vector<int> &labels) {
  const std::size_t n = labels.size();
  std::map<int, std::size_t> sizes;
  for (int l : labels)
    if (l != -1)
      ++sizes[l];
  if (sizes.size() < 2)
    return 0.0;
  double total = 0.0;
  std::size_t counted = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (labels[i] == -1)
      continue;
    std::map<int, double> sums;
    for (std::size_t j = 0; j < n; ++j)
      if (labels[j] != -1 && j != i)
        sums[labels[j]] += distance(points, d, i, j);
    ++counted;
    const std::size_t own = sizes[labels[i]];
    if (own <= 1)
      continue; // a point alone in its cluster: 0
    const double a = sums[labels[i]] / static_cast<double>(own - 1);
    double b = HUGE_VAL;
    for (const auto &[l, size] : sizes)
      if (l != labels[i])
        b = std::min(b, sums[l] / static_cast<double>(size));
    const double denominator = std::max(a, b);
    total += denominator > 0.0 ? (b - a) / denominator : 0.0;
  }
  return counted > 0 ? total / static_cast<double>(counted) : 0.0;
}

SignificantClusters significant_clusters(const CosymCoordinates &c,
                                         const ClusterOptions &options) {
  SignificantClusters out;
  const std::size_t n = c.n, d = c.dim;
  if (n == 0 || d == 0)
    return out;
  // The coordinates as minimised, row major: what DIALS clusters.
  std::vector<double> points(n * d);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t k = 0; k < d; ++k)
      points[i * d + k] = c.x[k * n + i];
  const auto clusters_in = [](const std::vector<int> &labels) {
    std::set<int> s;
    for (int l : labels)
      if (l != -1)
        s.insert(l);
    return s.size();
  };
  if (n < 5) {
    out.min_samples = out.initial_min_samples = n;
    out.model = optics(points, n, d, n, options.max_eps, options.xi);
    out.labels = out.model.labels;
    if (clusters_in(out.labels) > 1)
      out.silhouette = silhouette(points, d, out.labels);
    return out;
  }
  const std::size_t initial = std::max<std::size_t>(
      5,
      static_cast<std::size_t>(static_cast<double>(n) / static_cast<double>(d) *
                               options.min_points_buffer));
  out.initial_min_samples = initial;
  if (!options.optimise) {
    out.min_samples = initial;
    out.model = optics(points, n, d, initial, options.max_eps, options.xi);
    out.labels = out.model.labels;
    if (clusters_in(out.labels) > 1)
      out.silhouette = silhouette(points, d, out.labels);
    return out;
  }
  // Every min_samples DIALS's search ranges over, scored as its objective;
  // the best kept, a tie to the initial guess (DIALS's first trial), then the
  // smallest.
  bool any = false;
  double best = -HUGE_VAL;
  for (std::size_t ms = 5; ms <= std::min<std::size_t>(n, 50); ++ms) {
    OpticsResult model = optics(points, n, d, ms, options.max_eps, options.xi);
    bool all_inf = true;
    for (double r : model.reachability)
      all_inf = all_inf && std::isinf(r);
    if (all_inf)
      continue; // pruned
    any = true;
    std::size_t noise = 0;
    for (int l : model.labels)
      noise += l == -1;
    const double fraction = static_cast<double>(noise) / static_cast<double>(n);
    const double penalty =
        1.0 - options.alpha * std::pow(fraction, options.gamma);
    double score = -1.0, sil = 0.0;
    if (clusters_in(model.labels) > 1) {
      std::vector<double> kept;
      std::vector<int> kept_labels;
      for (std::size_t i = 0; i < n; ++i)
        if (model.labels[i] != -1) {
          kept_labels.push_back(model.labels[i]);
          for (std::size_t k = 0; k < d; ++k)
            kept.push_back(points[i * d + k]);
        }
      sil = silhouette(kept, d, kept_labels);
      score = sil * penalty;
    }
    out.scores.emplace_back(ms, score);
    const bool better = score > best || (score == best && ms == initial &&
                                         out.min_samples != initial);
    if (better) {
      best = score;
      out.min_samples = ms;
      out.score = score;
      out.silhouette = sil;
      out.labels = model.labels;
      out.model = std::move(model);
    }
  }
  if (!any) {
    // DIALS's default_optics: no max_eps.
    out.min_samples = initial;
    out.model = optics(points, n, d, initial, HUGE_VAL, options.xi);
    out.labels = out.model.labels;
    return out;
  }
  if (out.score >= 0.0)
    return out;
  // DIALS's single_cluster_optics, line for line.
  out.single_cluster = true;
  out.min_samples = initial;
  const OpticsResult model =
      optics(points, n, d, initial, options.max_eps, options.xi);
  out.model = model;
  std::vector<double> reach(n), core(n);
  for (std::size_t k = 0; k < n; ++k) {
    reach[k] = model.reachability[model.ordering[k]];
    core[k] = model.core_distance[model.ordering[k]];
  }
  std::vector<std::size_t> finite;
  for (std::size_t k = 0; k < n; ++k)
    if (std::isfinite(reach[k]))
      finite.push_back(k);
  std::vector<double> gradients(n, 0.0);
  for (std::size_t q = 1; q < finite.size(); ++q)
    gradients[finite[q]] = reach[finite[q]] - reach[finite[q - 1]];
  long start = 0;
  long end = -1; // None
  {
    long last_negative = -1, first_large = -1;
    for (std::size_t q = 0; q < finite.size(); ++q) {
      const double g = gradients[finite[q]];
      if (g < -options.xi)
        last_negative = static_cast<long>(q);
      if (g > options.xi && first_large < 0)
        first_large = static_cast<long>(q);
    }
    if (last_negative >= 0)
      start =
          static_cast<long>(finite[static_cast<std::size_t>(last_negative)]);
    if (first_large >= 0)
      end = static_cast<long>(finite[static_cast<std::size_t>(first_large)]);
  }
  const bool end_truthy = end > 0; // Python's `if end`: None and 0 false
  if (end_truthy && start > end)
    start = 0;
  start = start > 2 ? start - 2 : 0;
  std::vector<int> labels(n);
  for (std::size_t k = 0; k < n; ++k)
    labels[k] = model.labels[model.ordering[k]];
  for (long k = 0; k < start; ++k)
    labels[static_cast<std::size_t>(k)] = -1;
  if (end_truthy) {
    for (long k = start; k < end; ++k)
      labels[static_cast<std::size_t>(k)] = 0;
    for (long k = end; k < static_cast<long>(n); ++k)
      labels[static_cast<std::size_t>(k)] = -1;
  } else {
    for (long k = start; k < static_cast<long>(n); ++k)
      labels[static_cast<std::size_t>(k)] = 0;
  }
  for (std::size_t k = 0; k < n; ++k)
    if (std::isinf(core[k]))
      labels[k] = -1;
  long mystery = -1;
  for (std::size_t k = 0; k < n; ++k)
    if (std::isinf(reach[k]) && !std::isinf(core[k])) {
      mystery = static_cast<long>(k);
      break;
    }
  if (mystery >= 0) {
    std::vector<double> verified;
    double mystery_density = 0.0;
    for (std::size_t k = 0; k < n; ++k)
      if (labels[k] == 0) {
        if (static_cast<long>(k) != mystery)
          verified.push_back(core[k]);
        else
          mystery_density = core[k];
      }
    if (mystery_density != 0.0 && !verified.empty()) { // Python's truthiness
      double avg = 0.0;
      for (double v : verified)
        avg += v;
      avg /= static_cast<double>(verified.size());
      double var = 0.0;
      for (double v : verified)
        var += (v - avg) * (v - avg);
      const double sd = std::sqrt(var / static_cast<double>(verified.size()));
      labels[static_cast<std::size_t>(mystery)] =
          mystery_density <= avg + 3.0 * sd && mystery_density >= avg - 3.0 * sd
              ? 0
              : -1;
    }
  }
  out.labels.assign(n, -1);
  for (std::size_t k = 0; k < n; ++k)
    out.labels[model.ordering[k]] = labels[k];
  out.score = -1.0;
  out.silhouette = 0.0;
  return out;
}

std::vector<std::array<double, 4>> ward_linkage(const std::vector<double> &d,
                                                std::size_t n) {
  std::vector<std::array<double, 4>> rows;
  if (n < 2)
    return rows;
  // The working distances between active clusters; each cluster's id (a leaf
  // its index, a merge n + its row) and size.
  std::vector<double> dist = d;
  std::vector<std::size_t> id(n), size(n, 1);
  std::vector<bool> active(n, true);
  for (std::size_t i = 0; i < n; ++i)
    id[i] = i;
  for (std::size_t step = 0; step + 1 < n; ++step) {
    std::size_t bi = 0, bj = 0;
    double best = HUGE_VAL;
    for (std::size_t i = 0; i < n; ++i) {
      if (!active[i])
        continue;
      for (std::size_t j = i + 1; j < n; ++j)
        if (active[j] && dist[i * n + j] < best) {
          best = dist[i * n + j];
          bi = i;
          bj = j;
        }
    }
    const std::size_t ni = size[bi], nj = size[bj];
    rows.push_back({static_cast<double>(std::min(id[bi], id[bj])),
                    static_cast<double>(std::max(id[bi], id[bj])), best,
                    static_cast<double>(ni + nj)});
    // Lance-Williams for Ward, on distances: the merged cluster in bi's slot.
    for (std::size_t k = 0; k < n; ++k) {
      if (!active[k] || k == bi || k == bj)
        continue;
      const double nk = static_cast<double>(size[k]);
      const double dki = dist[k * n + bi], dkj = dist[k * n + bj];
      const double v = std::sqrt(((static_cast<double>(ni) + nk) * dki * dki +
                                  (static_cast<double>(nj) + nk) * dkj * dkj -
                                  nk * best * best) /
                                 (static_cast<double>(ni + nj) + nk));
      dist[k * n + bi] = dist[bi * n + k] = v;
    }
    size[bi] = ni + nj;
    id[bi] = n + step;
    active[bj] = false;
  }
  return rows;
}

Dendrogram dendrogram(const std::vector<std::array<double, 4>> &linkage,
                      std::size_t n) {
  Dendrogram out;
  if (n == 0)
    return out;
  if (linkage.empty()) {
    out.leaves.push_back(0);
    return out;
  }
  // Each node's x and height, by a walk from the root, left (the first of the
  // row) before right, as scipy's.
  const std::size_t root = n + linkage.size() - 1;
  std::vector<double> x(n + linkage.size(), 0.0), h(n + linkage.size(), 0.0);
  std::function<void(std::size_t)> walk = [&](std::size_t node) {
    if (node < n) {
      x[node] = 10.0 * static_cast<double>(out.leaves.size()) + 5.0;
      out.leaves.push_back(node);
      return;
    }
    const auto &row = linkage[node - n];
    const auto a = static_cast<std::size_t>(row[0]),
               b = static_cast<std::size_t>(row[1]);
    walk(a);
    walk(b);
    x[node] = 0.5 * (x[a] + x[b]);
    h[node] = row[2];
    out.x.push_back({x[a], x[a], x[b], x[b]});
    out.y.push_back({h[a], row[2], row[2], h[b]});
  };
  walk(root);
  return out;
}

} // namespace mxi
