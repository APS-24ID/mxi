// Intensity-based clustering of isomorphous data sets: see correlation.hh.

#include "correlation.hh"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <random>

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

} // namespace mxi
