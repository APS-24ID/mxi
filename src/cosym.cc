// cosym: see cosym.hh.

#include "cosym.hh"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
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

// A symmetric matrix's eigenvalues and eigenvectors, by Jacobi's rotations:
// the matrices here are a few hundred on a side at most.
void eigen(std::vector<double> a, std::size_t n, std::vector<double> *values,
           std::vector<double> *vectors) {
  std::vector<double> v(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    v[i * n + i] = 1.0;
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0;
    for (std::size_t p = 0; p < n; ++p)
      for (std::size_t q = p + 1; q < n; ++q)
        off += a[p * n + q] * a[p * n + q];
    if (off < 1e-20)
      break;
    for (std::size_t p = 0; p < n; ++p)
      for (std::size_t q = p + 1; q < n; ++q) {
        const double apq = a[p * n + q];
        if (std::fabs(apq) < 1e-15)
          continue;
        const double theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
        const double t = (theta >= 0 ? 1.0 : -1.0) /
                         (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
        for (std::size_t k = 0; k < n; ++k) {
          const double akp = a[k * n + p], akq = a[k * n + q];
          a[k * n + p] = c * akp - s * akq;
          a[k * n + q] = s * akp + c * akq;
        }
        for (std::size_t k = 0; k < n; ++k) {
          const double apk = a[p * n + k], aqk = a[q * n + k];
          a[p * n + k] = c * apk - s * aqk;
          a[q * n + k] = s * apk + c * aqk;
        }
        for (std::size_t k = 0; k < n; ++k) {
          const double vkp = v[k * n + p], vkq = v[k * n + q];
          v[k * n + p] = c * vkp - s * vkq;
          v[k * n + q] = s * vkp + c * vkq;
        }
      }
  }
  values->assign(n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    (*values)[i] = a[i * n + i];
  *vectors = std::move(v); // column j is the eigenvector of value j
}

} // namespace

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
  std::vector<double> cc(n * n * m, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t t = 0; t < m; ++t)
        cc[(i * n + j) * m + t] =
            correlation(maps[i], maps[j], operations[t]).first;

  // The matrix r_ik,jl, its diagonal one; its leading eigenvectors.
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
  eigen(r, size, &values, &vectors);
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
