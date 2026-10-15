// cosym: see cosym.hh.

#include "cosym.hh"
#include "parallel.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace mxi {

namespace {

using Key = std::array<int, 3>;

// A reflection and its Friedel mate one: the one whose first non-zero index
// is positive.
Key friedel(const Miller &h) {
  Key k{h[0], h[1], h[2]};
  for (int x : k) {
    if (x > 0)
      break;
    if (x < 0) {
      k = {-k[0], -k[1], -k[2]};
      break;
    }
  }
  return k;
}

std::map<Key, double> by_reflection(const P1Intensities &d) {
  std::map<Key, std::pair<double, int>> sum;
  for (std::size_t i = 0; i < d.size(); ++i) {
    auto &s = sum[friedel(d.hkl[i])];
    s.first += d.i[i];
    s.second += 1;
  }
  std::map<Key, double> out;
  for (const auto &[k, s] : sum)
    out[k] = s.first / s.second;
  return out;
}

Rotation multiply(const Rotation &a, const Rotation &b) {
  Rotation c{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      int v = 0;
      for (int k = 0; k < 3; ++k)
        v += a[3 * i + k] * b[3 * k + j];
      c[3 * i + j] = v;
    }
  return c;
}

// Pearson's correlation of a's values with b's at t applied: over the
// reflections of a whose image under t b has.
std::pair<double, std::size_t> correlation(const std::map<Key, double> &a,
                                           const std::map<Key, double> &b,
                                           const Rotation &t) {
  double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  std::size_t n = 0;
  for (const auto &[k, x] : a) {
    const Miller h{k[0], k[1], k[2]};
    const auto it = b.find(friedel(apply(t, h)));
    if (it == b.end())
      continue;
    const double y = it->second;
    sx += x;
    sy += y;
    sxx += x * x;
    syy += y * y;
    sxy += x * y;
    ++n;
  }
  if (n < 3)
    return {0.0, n};
  const double cov = sxy - sx * sy / static_cast<double>(n);
  const double vx = sxx - sx * sx / static_cast<double>(n);
  const double vy = syy - sy * sy / static_cast<double>(n);
  if (!(vx > 0.0) || !(vy > 0.0))
    return {0.0, n};
  return {cov / std::sqrt(vx * vy), n};
}

// Each data set's mean intensities by a packed index, for looking up in the
// correlations of every pair under every rotation: a std::map's lookup, 31104
// correlations of 20000 reflections for 36 crystals, took two minutes of
// mxi_symmetry's time on Graeme's cows, pigs and people.
std::uint64_t pack(const Key &k) {
  const auto field = [](int v) {
    return static_cast<std::uint64_t>(static_cast<std::int64_t>(v) +
                                      (1 << 20)) &
           ((std::uint64_t{1} << 21) - 1);
  };
  return (field(k[0]) << 42) | (field(k[1]) << 21) | field(k[2]);
}

using Lookup = std::unordered_map<std::uint64_t, double>;

Lookup lookup_of(const std::map<Key, double> &m) {
  Lookup out;
  out.reserve(2 * m.size());
  for (const auto &[k, v] : m)
    out.emplace(pack(k), v);
  return out;
}

// The same as correlation, b looked up by its packed index: the same pairs, in
// the same order, so the same sums.
std::pair<double, std::size_t> correlation(const std::map<Key, double> &a,
                                           const Lookup &b, const Rotation &t) {
  double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  std::size_t n = 0;
  for (const auto &[k, x] : a) {
    const Miller h{k[0], k[1], k[2]};
    const auto it = b.find(pack(friedel(apply(t, h))));
    if (it == b.end())
      continue;
    const double y = it->second;
    sx += x;
    sy += y;
    sxx += x * x;
    syy += y * y;
    sxy += x * y;
    ++n;
  }
  if (n < 3)
    return {0.0, n};
  const double cov = sxy - sx * sy / static_cast<double>(n);
  const double vx = sxx - sx * sx / static_cast<double>(n);
  const double vy = syy - sy * sy / static_cast<double>(n);
  if (!(vx > 0.0) || !(vy > 0.0))
    return {0.0, n};
  return {cov / std::sqrt(vx * vy), n};
}

// A symmetric matrix's eigenvalues and eigenvectors (column q of `vectors`,
// row major, the q-th): Householder's reduction to tridiagonal form, then the
// implicit QL algorithm, as JAMA's (public domain, from EISPACK's tred2 and
// tql2). Once O(n^3), where Jacobi's rotations, as first written, were O(n^3)
// a sweep -- 864 x 864 for 36 crystals and 24 rotations, 83 s of it. The
// eigenvalues the same to rounding; an eigenvector's sign may differ, which no
// cosine sees.
void eigen(std::vector<double> a, std::size_t n, std::vector<double> *values,
           std::vector<double> *vectors) {
  values->assign(n, 0.0);
  vectors->assign(n * n, 0.0);
  if (n == 0)
    return;
  std::vector<double> &V = a; // V[i][j] = V[i * n + j]
  std::vector<double> d(n), e(n, 0.0);
  const auto at = [&](std::size_t i, std::size_t j) -> double & {
    return V[i * n + j];
  };
  // tred2
  for (std::size_t j = 0; j < n; ++j)
    d[j] = at(n - 1, j);
  for (std::size_t i = n - 1; i > 0; --i) {
    double scale = 0.0, h = 0.0;
    for (std::size_t k = 0; k < i; ++k)
      scale += std::fabs(d[k]);
    if (scale == 0.0) {
      e[i] = d[i - 1];
      for (std::size_t j = 0; j < i; ++j) {
        d[j] = at(i - 1, j);
        at(i, j) = 0.0;
        at(j, i) = 0.0;
      }
    } else {
      for (std::size_t k = 0; k < i; ++k) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      double f = d[i - 1];
      double g = std::sqrt(h);
      if (f > 0)
        g = -g;
      e[i] = scale * g;
      h -= f * g;
      d[i - 1] = f - g;
      for (std::size_t j = 0; j < i; ++j)
        e[j] = 0.0;
      for (std::size_t j = 0; j < i; ++j) {
        f = d[j];
        at(j, i) = f;
        g = e[j] + at(j, j) * f;
        for (std::size_t k = j + 1; k + 1 <= i; ++k) {
          g += at(k, j) * d[k];
          e[k] += at(k, j) * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (std::size_t j = 0; j < i; ++j) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      const double hh = f / (h + h);
      for (std::size_t j = 0; j < i; ++j)
        e[j] -= hh * d[j];
      for (std::size_t j = 0; j < i; ++j) {
        f = d[j];
        g = e[j];
        for (std::size_t k = j; k + 1 <= i; ++k)
          at(k, j) -= (f * e[k] + g * d[k]);
        d[j] = at(i - 1, j);
        at(i, j) = 0.0;
      }
    }
    d[i] = h;
  }
  for (std::size_t i = 0; i + 1 < n; ++i) {
    at(n - 1, i) = at(i, i);
    at(i, i) = 1.0;
    const double h = d[i + 1];
    if (h != 0.0) {
      for (std::size_t k = 0; k <= i; ++k)
        d[k] = at(k, i + 1) / h;
      for (std::size_t j = 0; j <= i; ++j) {
        double g = 0.0;
        for (std::size_t k = 0; k <= i; ++k)
          g += at(k, i + 1) * at(k, j);
        for (std::size_t k = 0; k <= i; ++k)
          at(k, j) -= g * d[k];
      }
    }
    for (std::size_t k = 0; k <= i; ++k)
      at(k, i + 1) = 0.0;
  }
  for (std::size_t j = 0; j < n; ++j) {
    d[j] = at(n - 1, j);
    at(n - 1, j) = 0.0;
  }
  at(n - 1, n - 1) = 1.0;
  e[0] = 0.0;
  // tql2, on W = V transposed: a rotation of columns i, i + 1 of V is of rows
  // i, i + 1 of W, contiguous.
  std::vector<double> W(n * n);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t k = 0; k < n; ++k)
      W[i * n + k] = at(k, i);
  for (std::size_t i = 1; i < n; ++i)
    e[i - 1] = e[i];
  e[n - 1] = 0.0;
  double f = 0.0, tst1 = 0.0;
  const double eps = std::pow(2.0, -52.0);
  for (std::size_t l = 0; l < n; ++l) {
    tst1 = std::max(tst1, std::fabs(d[l]) + std::fabs(e[l]));
    std::size_t m = l;
    while (m < n) {
      if (std::fabs(e[m]) <= eps * tst1)
        break;
      ++m;
    }
    if (m == n)
      m = n - 1; // e[n - 1] is zero: not reached
    if (m > l) {
      int iterations = 0;
      do {
        if (++iterations > 300)
          throw std::runtime_error("cosym: the eigen-decomposition did not "
                                   "converge");
        double g = d[l];
        double p = (d[l + 1] - g) / (2.0 * e[l]);
        double r = std::hypot(p, 1.0);
        if (p < 0)
          r = -r;
        d[l] = e[l] / (p + r);
        d[l + 1] = e[l] * (p + r);
        const double dl1 = d[l + 1];
        double h = g - d[l];
        for (std::size_t i = l + 2; i < n; ++i)
          d[i] -= h;
        f += h;
        p = d[m];
        double c = 1.0, c2 = c, c3 = c;
        const double el1 = e[l + 1];
        double s = 0.0, s2 = 0.0;
        for (std::size_t ii = m; ii-- > l;) {
          const std::size_t i = ii;
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = std::hypot(p, e[i]);
          e[i + 1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i + 1] = h + s * (c * g + s * d[i]);
          double *wi = &W[i * n], *wi1 = &W[(i + 1) * n];
          for (std::size_t k = 0; k < n; ++k) {
            h = wi1[k];
            wi1[k] = s * wi[k] + c * h;
            wi[k] = c * wi[k] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
      } while (std::fabs(e[l]) > eps * tst1);
    }
    d[l] += f;
    e[l] = 0.0;
  }
  *values = d;
  for (std::size_t q = 0; q < n; ++q)
    for (std::size_t k = 0; k < n; ++k)
      (*vectors)[k * n + q] = W[q * n + k];
}

} // namespace

void symmetric_eigen(const std::vector<double> &a, std::size_t n,
                     std::vector<double> *values,
                     std::vector<double> *vectors) {
  eigen(a, n, values, vectors);
}

CosymResult cosym(const std::vector<P1Intensities> &datasets,
                  const std::vector<Rotation> &operations) {
  CosymResult result;
  result.operations = operations;
  const std::size_t n = datasets.size(), m = operations.size();
  if (n == 0 || m == 0)
    return result;

  // The group's product table: which operation is k's inverse times l, the
  // relative rotation between data set i under k and j under l.
  std::vector<std::size_t> inverse(m, 0), relative(m * m, 0);
  const auto index_of = [&](const Rotation &r) {
    for (std::size_t q = 0; q < m; ++q)
      if (operations[q] == r)
        return q;
    throw std::runtime_error("cosym: the lattice's rotations are not a group");
  };
  const Rotation identity{1, 0, 0, 0, 1, 0, 0, 0, 1};
  for (std::size_t k = 0; k < m; ++k)
    for (std::size_t q = 0; q < m; ++q)
      if (multiply(operations[k], operations[q]) == identity)
        inverse[k] = q;
  for (std::size_t k = 0; k < m; ++k)
    for (std::size_t l = 0; l < m; ++l)
      relative[k * m + l] =
          index_of(multiply(operations[inverse[k]], operations[l]));

  // Each pair of data sets' correlation under each relative rotation.
  std::vector<std::map<Key, double>> maps;
  for (const P1Intensities &d : datasets)
    maps.push_back(by_reflection(d));
  const auto clock = [] {
    return std::chrono::duration<double>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
  };
  const double t_cc = clock();
  std::vector<Lookup> lookups(n);
  for_each_index(n, [&](std::size_t j) { lookups[j] = lookup_of(maps[j]); });
  std::vector<double> cc(n * n * m, 0.0);
  for_each_index(n * n, [&](std::size_t ij) {
    const std::size_t i = ij / n, j = ij % n;
    for (std::size_t t = 0; t < m; ++t)
      cc[ij * m + t] = correlation(maps[i], lookups[j], operations[t]).first;
  });

  // The matrix r_ik,jl, its diagonal one; its leading eigenvectors.
  result.seconds_correlations = clock() - t_cc;
  const std::size_t size = n * m;
  std::vector<double> r(size * size, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t k = 0; k < m; ++k)
      for (std::size_t j = 0; j < n; ++j)
        for (std::size_t l = 0; l < m; ++l) {
          const std::size_t a = i * m + k, b = j * m + l;
          r[a * size + b] =
              a == b ? 1.0 : cc[(i * n + j) * m + relative[k * m + l]];
        }
  std::vector<double> values, vectors;
  const double t_eigen = clock();
  eigen(r, size, &values, &vectors);
  result.seconds_eigen = clock() - t_eigen;
  std::vector<std::size_t> order(size);
  for (std::size_t q = 0; q < size; ++q)
    order[q] = q;
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return values[a] > values[b];
  });
  double positive = 0.0;
  for (double v : values)
    positive += std::max(v, 0.0);
  const double largest = values[order[0]];
  for (std::size_t q = 0; q < std::min<std::size_t>(size, 4); ++q)
    result.explained.push_back(
        positive > 0.0 ? std::max(values[order[q]], 0.0) / positive : 0.0);
  // Components holding a tenth of the largest's variance, or more.
  std::size_t d = 0;
  while (d < size && values[order[d]] >= 0.1 * largest)
    ++d;
  d = std::max<std::size_t>(d, 1);
  result.dimensions = d;
  const auto point = [&](std::size_t a) {
    std::vector<double> x(d);
    for (std::size_t q = 0; q < d; ++q)
      x[q] = vectors[a * size + order[q]] *
             std::sqrt(std::max(values[order[q]], 0.0));
    return x;
  };
  const auto cosine = [](const std::vector<double> &a,
                         const std::vector<double> &b) {
    double ab = 0, aa = 0, bb = 0;
    for (std::size_t q = 0; q < a.size(); ++q) {
      ab += a[q] * b[q];
      aa += a[q] * a[q];
      bb += b[q] * b[q];
    }
    return aa > 0 && bb > 0 ? ab / std::sqrt(aa * bb) : 0.0;
  };

  // Each data set: the rotation whose point is nearest the reference, the
  // first's identity; the identity kept where it is as near, within 0.1.
  const std::vector<double> reference = point(0);
  for (std::size_t j = 0; j < n; ++j) {
    std::size_t best = 0;
    double best_cos = -2.0;
    for (std::size_t l = 0; l < m; ++l) {
      const double c = cosine(reference, point(j * m + l));
      if (c > best_cos) {
        best_cos = c;
        best = l;
      }
    }
    const double at_identity = cosine(reference, point(j * m + 0));
    if (at_identity >= best_cos - 0.1) {
      best = 0;
      best_cos = at_identity;
    }
    result.chosen.push_back(best);
    result.cosine.push_back(best_cos);
    const auto c = correlation(maps[0], maps[j], operations[best]);
    result.cc.push_back(c.first);
    result.common.push_back(c.second);
  }
  return result;
}

} // namespace mxi
