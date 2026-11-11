#include "triage.hh"

#include "laue.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>

namespace mxi {

namespace {

Rotation multiply(const Rotation &a, const Rotation &b) {
  Rotation c{};
  for (std::int32_t i = 0; i < 3; ++i)
    for (std::int32_t j = 0; j < 3; ++j) {
      std::int32_t s = 0;
      for (std::int32_t k = 0; k < 3; ++k)
        s += a[i * 3 + k] * b[k * 3 + j];
      c[i * 3 + j] = s;
    }
  return c;
}

std::int32_t determinant(const Rotation &r) {
  return r[0] * (r[4] * r[8] - r[5] * r[7]) -
         r[1] * (r[3] * r[8] - r[5] * r[6]) +
         r[2] * (r[3] * r[7] - r[4] * r[6]);
}

std::array<double, 9> metric(const UnitCell &c) {
  const double deg = std::acos(-1.0) / 180.0;
  const double ab = c.a * c.b * std::cos(c.gamma * deg),
               ac = c.a * c.c * std::cos(c.beta * deg),
               bc = c.b * c.c * std::cos(c.alpha * deg);
  return {c.a * c.a, ab, ac, ab, c.b * c.b, bc, ac, bc, c.c * c.c};
}

std::vector<Rotation> point_group(const SpaceGroup &group) {
  std::set<Rotation> seen;
  for (const auto &r : group.rotations())
    seen.insert(r);
  return {seen.begin(), seen.end()};
}

Miller add(const Miller &a, const Miller &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

// Chebyshev nodes as cctbx's chebyshev_nodes(n, low, high, include_limits=
// True): -cos((2k - 1) pi / 2n), stretched to reach the limits.
std::vector<double> chebyshev_nodes(std::int32_t n, double low, double high) {
  std::vector<double> x(static_cast<std::size_t>(n));
  for (std::int32_t k = 0; k < n; ++k)
    x[k] = -std::cos((2.0 * (k + 1) - 1.0) / n * std::acos(-1.0) / 2.0);
  const double span = (x.back() - x.front()) / 2.0;
  for (double &v : x)
    v = 0.5 * (low + high) + 0.5 * (high - low) * (v / span);
  x.front() = low;
  x.back() = high;
  return x;
}

// A Chebyshev series on [low, high], fitted by least squares to (x, y).
struct Chebyshev {
  double low = 0, high = 1;
  std::vector<double> c;
  double operator()(double x) const {
    const double t = (2.0 * x - low - high) / (high - low);
    double t0 = 1, t1 = t, sum = c.empty() ? 0.0 : c[0];
    for (std::size_t k = 1; k < c.size(); ++k) {
      sum += c[k] * t1;
      const double t2 = 2 * t * t1 - t0;
      t0 = t1;
      t1 = t2;
    }
    return sum;
  }
};

Chebyshev chebyshev_fit(std::int32_t terms, const std::vector<double> &x,
                        const std::vector<double> &y, double low, double high) {
  Chebyshev f{low, high, {}};
  const std::size_t m = static_cast<std::size_t>(terms);
  std::vector<double> a(m * m, 0.0), b(m, 0.0), basis(m);
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double t = (2.0 * x[i] - low - high) / (high - low);
    basis[0] = 1;
    if (m > 1)
      basis[1] = t;
    for (std::size_t k = 2; k < m; ++k)
      basis[k] = 2 * t * basis[k - 1] - basis[k - 2];
    for (std::size_t p = 0; p < m; ++p) {
      b[p] += basis[p] * y[i];
      for (std::size_t q = 0; q < m; ++q)
        a[p * m + q] += basis[p] * basis[q];
    }
  }
  // The normal equations by Gaussian elimination with partial pivoting.
  for (std::size_t col = 0; col < m; ++col) {
    std::size_t pivot = col;
    for (std::size_t r = col + 1; r < m; ++r)
      if (std::abs(a[r * m + col]) > std::abs(a[pivot * m + col]))
        pivot = r;
    for (std::size_t q = 0; q < m; ++q)
      std::swap(a[col * m + q], a[pivot * m + q]);
    std::swap(b[col], b[pivot]);
    const double d = a[col * m + col];
    if (std::abs(d) < 1e-300)
      continue;
    for (std::size_t r = 0; r < m; ++r) {
      if (r == col)
        continue;
      const double factor = a[r * m + col] / d;
      for (std::size_t q = 0; q < m; ++q)
        a[r * m + q] -= factor * a[col * m + q];
      b[r] -= factor * b[col];
    }
  }
  f.c.resize(m);
  for (std::size_t k = 0; k < m; ++k)
    f.c[k] = std::abs(a[k * m + k]) < 1e-300 ? 0.0 : b[k] / a[k * m + k];
  return f;
}

// The reflections in [d_min, d_max], and each one's E^2, normalised as
// xtriage's kernel_normalisation (absolute_scaling.py, auto_kernel): <I/eps>
// by a Gaussian kernel -- its width the d*^2 of the 50th reflection from low
// resolution less the least -- at 23 Chebyshev nodes over the range, a
// 13-term Chebyshev series fitted to its log, and E^2 = I / that, eps NOT
// divided out (its normalised_miller; normalised_miller_dev_eps divides).
struct Used {
  std::vector<std::size_t> rows; //: in range, every one
  std::vector<double> e2;        //: each one's E^2
  std::vector<bool> centric;
};

Used normalised(const TriageData &data, const UnitCell &cell,
                const SpaceGroup &group, double d_min, double d_max) {
  Used u;
  std::vector<double> s2, scaled;
  for (std::size_t r = 0; r < data.hkl.size(); ++r) {
    const double d = d_spacing(cell, data.hkl[r]);
    if (!(d >= d_min && d <= d_max) || group.absent(data.hkl[r]))
      continue;
    u.rows.push_back(r);
    u.centric.push_back(group.centric(data.hkl[r]));
    s2.push_back(1.0 / (d * d));
    scaled.push_back(data.intensity[r] / epsilon(group, data.hkl[r]));
  }
  const std::size_t n = u.rows.size();
  u.e2.assign(n, 0.0);
  if (n < 2)
    return u;
  std::vector<double> sorted = s2;
  std::sort(sorted.begin(), sorted.end());
  const double low = sorted.front(), high = sorted.back();
  std::size_t number = std::min<std::size_t>(50, n - 1);
  double width = sorted[number] - low;
  while (width <= 0 && number + 50 < n / 2) {
    number += 50;
    width = sorted[number] - low;
  }
  if (!(width > 0) || !(high > low))
    return u;
  const std::vector<double> nodes = chebyshev_nodes(23, low, high);
  std::vector<double> x, y;
  for (double node : nodes) {
    double w = 0, sum = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const double dx = s2[i] - node;
      const double k = std::exp(-dx * dx / (2 * width * width));
      w += k;
      sum += k * scaled[i];
    }
    const double mean = sum / std::max(w, 1e-8);
    if (mean > 1e-16) {
      x.push_back(node);
      y.push_back(std::log(mean));
    }
  }
  if (x.size() < 2)
    return u;
  const Chebyshev fit = chebyshev_fit(13, x, y, low, high);
  for (std::size_t i = 0; i < n; ++i)
    // I / (eps <I/eps>): eps divided out, as E^2 is defined -- xtriage's
    // normalised_miller leaves it in, the few acentric reflections on an
    // axis eps times too strong (a perfect twin in P 4: 1.58, not 1.5).
    u.e2[i] = scaled[i] / std::exp(fit(s2[i]));
  return u;
}

Moments moments(const Used &u, bool centric) {
  Moments m;
  double e2 = 0, e4 = 0, e = 0, dev = 0;
  for (std::size_t i = 0; i < u.e2.size(); ++i) {
    // Every E^2, negative ones too: a weak reflection measured below zero is
    // still a weak reflection. xtriage keeps the positive only, which takes
    // the weakest out -- as twinning does -- and lowers N(z) at small z.
    if (u.centric[i] != centric)
      continue;
    const double x = u.e2[i];
    e2 += x;
    e4 += x * x;
    e += std::sqrt(std::max(x, 0.0)); // F from I: negatives as nothing
    dev += std::abs(x - 1.0);
    ++m.n;
  }
  if (m.n == 0)
    return m;
  const double n = static_cast<double>(m.n);
  m.i2_over_i_squared = (e4 / n) / ((e2 / n) * (e2 / n));
  m.f_squared_over_f2 = (e / n) * (e / n) / (e2 / n);
  m.e2_minus_one = dev / n;
  return m;
}

// The data as xtriage gives them to its twin analysis (data_statistics.py,
// wilson_scaling): over d*^2 0.008 to 0.690, put on an absolute scale with
// the anisotropic fall-off removed -- anisotropic_correction with its ML
// fit's scale and B; here the same model's fit, laue.hh's, and its f2(s) left
// in as there -- then possible_outliers removed: E^2 by kernel normalisation,
// divided by eps, an outlier where its Wilson probability is below 1e-6 or
// its extreme-value probability, the chance that the largest of n is as
// large, below 0.1. Sharpened, the H test's strongest half is chosen across
// resolution, not from the low resolution the raw intensities favour.
TriageData prepared(const TriageData &given, const UnitCell &cell,
                    const SpaceGroup &group, Triage *t) {
  const double d_lo = std::sqrt(1.0 / 0.008), d_hi = std::sqrt(1.0 / 0.690);
  TriageData window;
  P1Intensities p1;
  for (std::size_t r = 0; r < given.hkl.size(); ++r) {
    const double d = d_spacing(cell, given.hkl[r]);
    if (d > d_lo || d < d_hi)
      continue;
    window.hkl.push_back(given.hkl[r]);
    window.intensity.push_back(given.intensity[r]);
    window.sigma.push_back(given.sigma[r]);
    p1.hkl.push_back(given.hkl[r]);
    p1.i.push_back(given.intensity[r]);
    p1.sigma.push_back(given.sigma[r]);
    p1.d.push_back(d);
  }
  WilsonFit fit;
  P1Intensities scratch = p1;
  if (!p1.hkl.empty())
    normalise(scratch, 200, &fit);
  t->sharpened = fit.fitted;
  if (fit.fitted)
    for (std::size_t r = 0; r < window.hkl.size(); ++r) {
      const Miller &m = window.hkl[r];
      const double q[6] = {double(m[0]) * m[0], double(m[1]) * m[1],
                           double(m[2]) * m[2], 2.0 * m[0] * m[1],
                           2.0 * m[0] * m[2],   2.0 * m[1] * m[2]};
      double qb = 0.0;
      for (std::int32_t j = 0; j < 6; ++j)
        qb += q[j] * fit.b[j];
      const double factor = std::exp(qb - fit.log_scale);
      window.intensity[r] *= factor;
      window.sigma[r] *= factor;
    }
  // possible_outliers: over the whole window, E^2 / eps.
  const Used all = normalised(window, cell, group, 0.0, 1e9);
  std::size_t n_acentric = 0, n_centric = 0;
  for (bool c : all.centric)
    (c ? n_centric : n_acentric) += 1;
  std::vector<bool> outlier(window.hkl.size(), false);
  for (std::size_t i = 0; i < all.rows.size(); ++i) {
    const std::size_t r = all.rows[i];
    const double e2 = std::max(all.e2[i], 0.0); // eps divided out already
    double single, extreme;
    if (all.centric[i]) {
      const double cdf = std::erf(std::sqrt(e2 / 2.0));
      single = 1.0 - cdf;
      extreme = 1.0 - std::pow(cdf, static_cast<double>(n_centric));
    } else {
      single = std::exp(-e2);
      extreme =
          1.0 - std::pow(1.0 - std::exp(-e2), static_cast<double>(n_acentric));
    }
    outlier[r] = extreme < 0.1 || single < 1e-6;
  }
  TriageData out;
  for (std::size_t r = 0; r < window.hkl.size(); ++r) {
    if (outlier[r]) {
      ++t->outliers;
      continue;
    }
    out.hkl.push_back(window.hkl[r]);
    out.intensity.push_back(window.intensity[r]);
    out.sigma.push_back(window.sigma[r]);
  }
  return out;
}

} // namespace

std::int32_t epsilon(const SpaceGroup &group, const Miller &hkl) {
  std::int32_t n = 0;
  for (const Rotation &r : point_group(group))
    if (mxi::apply(r, hkl) == hkl)
      ++n;
  return std::max(n, 1);
}

std::string twin_law_text(const Rotation &r) {
  // h' = R^T h: component i of h' is sum_j R[j][i] h_j.
  const char *names = "hkl";
  std::string out;
  for (std::int32_t i = 0; i < 3; ++i) {
    std::string term;
    for (std::int32_t j = 0; j < 3; ++j) {
      const std::int32_t c = r[j * 3 + i];
      if (c == 0)
        continue;
      if (c < 0)
        term += "-";
      else if (!term.empty())
        term += "+";
      if (std::abs(c) != 1)
        term += std::to_string(std::abs(c));
      term += names[j];
    }
    out += (i ? "," : "") + (term.empty() ? std::string("0") : term);
  }
  return out;
}

std::vector<TwinLaw> twin_laws(const UnitCell &cell, const SpaceGroup &group,
                               double max_delta) {
  const std::vector<Rotation> points = point_group(group);
  const std::set<Rotation> in_group(points.begin(), points.end());
  const auto g = metric(cell);
  double norm = 0.0;
  for (double v : g)
    norm += v * v;
  std::set<std::set<Rotation>> cosets;
  std::vector<TwinLaw> out;
  for (const Rotation &l : lattice_symmetry(cell, max_delta)) {
    if (in_group.count(l) || determinant(l) != 1)
      continue;
    std::set<Rotation> coset;
    for (const Rotation &p : points)
      coset.insert(multiply(l, p));
    if (!cosets.insert(coset).second)
      continue;
    // Of the coset, a two-fold where there is one -- as xtriage names a twin
    // law -- else its least.
    Rotation pick = *coset.begin();
    for (const Rotation &c : coset)
      if (c[0] + c[4] + c[8] == -1) {
        pick = c;
        break;
      }
    TwinLaw law;
    law.op = pick;
    law.text = twin_law_text(pick);
    // G' = R^T G R: the metric the twin law would make of the cell.
    std::array<double, 9> moved{};
    for (std::int32_t i = 0; i < 3; ++i)
      for (std::int32_t j = 0; j < 3; ++j) {
        double s = 0.0;
        for (std::int32_t a = 0; a < 3; ++a)
          for (std::int32_t b = 0; b < 3; ++b)
            s += pick[a * 3 + i] * g[a * 3 + b] * pick[b * 3 + j];
        moved[i * 3 + j] = s;
      }
    double diff = 0.0;
    for (std::int32_t k = 0; k < 9; ++k)
      diff += (moved[k] - g[k]) * (moved[k] - g[k]);
    law.r_metric = norm > 0 ? diff / norm : 0.0;
    law.merohedral = law.r_metric < 1e-8;
    out.push_back(law);
  }
  return out;
}

double twin_test_d_min(const TriageData &data, const UnitCell &cell,
                       const SpaceGroup &group, const TriageOptions &o) {
  if (data.hkl.empty())
    return o.resolution_at_least;
  double d_lo = 0.0, d_hi = 1e9;
  for (const Miller &m : data.hkl) {
    const double d = d_spacing(cell, m);
    d_lo = std::max(d_lo, d);
    d_hi = std::min(d_hi, d);
  }
  // Shells of equal volume over the data's range, as cctbx's binner sets
  // them up for n bins: equal steps in d*^3.
  const std::int32_t bins = std::max(o.strength_bins, 2);
  const double lo = 1.0 / (d_lo * d_lo * d_lo), hi = 1.0 / (d_hi * d_hi * d_hi);
  auto bin_of = [&](double d) {
    const double v = 1.0 / (d * d * d);
    return std::clamp(static_cast<std::int32_t>((v - lo) / (hi - lo) * bins), 0,
                      bins - 1);
  };
  std::vector<double> possible(bins, 0.0), strong(bins, 0.0);
  const std::int32_t hm = static_cast<std::int32_t>(cell.a / d_hi) + 1,
                     km = static_cast<std::int32_t>(cell.b / d_hi) + 1,
                     lm = static_cast<std::int32_t>(cell.c / d_hi) + 1;
  std::set<Miller> counted;
  for (std::int32_t h = -hm; h <= hm; ++h)
    for (std::int32_t k = -km; k <= km; ++k)
      for (std::int32_t l = -lm; l <= lm; ++l) {
        const Miller m{h, k, l};
        if (h == 0 && k == 0 && l == 0)
          continue;
        const double d = d_spacing(cell, m);
        if (d < d_hi || d > d_lo || group.absent(m))
          continue;
        if (counted.insert(group.unique(m)).second)
          possible[bin_of(d)] += 1.0;
      }
  for (std::size_t r = 0; r < data.hkl.size(); ++r)
    if (data.intensity[r] > o.isigi_cut * data.sigma[r])
      strong[bin_of(d_spacing(cell, data.hkl[r]))] += 1.0;
  // From low resolution, of the shells finer than resolution_at_least, to
  // the last complete enough before the first that is not. xtriage takes the
  // finest complete enough wherever it lies, a shell beyond a failing one by
  // chance, and its loop leaves the last shell out; here each counts.
  double cut = 4.0;
  for (std::int32_t b = 0; b < bins; ++b) {
    const double d_edge = std::cbrt(1.0 / (lo + (hi - lo) * (b + 1) / bins));
    if (d_edge >= o.resolution_at_least)
      continue;
    if (!(possible[b] > 0 && strong[b] / possible[b] > o.completeness))
      break;
    cut = d_edge;
  }
  return cut > o.resolution_at_least ? o.resolution_at_least : cut;
}

Triage triage(const TriageData &given, const UnitCell &cell,
              const SpaceGroup &group, const TriageOptions &options) {
  Triage t;
  t.d_max = options.d_max;
  // The range from the data as given, as xtriage chooses it; every test on
  // the data as xtriage prepares them.
  t.d_min = options.d_min > 0 ? options.d_min
                              : twin_test_d_min(given, cell, group, options);
  const TriageData data = prepared(given, cell, group, &t);
  const Used u = normalised(data, cell, group, t.d_min, t.d_max);
  for (std::size_t i = 0; i < u.rows.size(); ++i)
    (u.centric[i] ? t.centric : t.acentric) += 1;
  t.used = t.acentric + t.centric;
  t.acentric_moments = moments(u, false);
  t.centric_moments = moments(u, true);

  // N(z), z = 0.1 .. 1.0.
  for (std::int32_t k = 1; k <= 10; ++k) {
    const double z = 0.1 * k;
    double na = 0, nc = 0;
    for (std::size_t i = 0; i < u.e2.size(); ++i)
      if (u.e2[i] <= z)
        (u.centric[i] ? nc : na) += 1;
    const double fa = t.acentric ? na / t.acentric : 0.0,
                 fc = t.centric ? nc / t.centric : 0.0;
    const double ta = 1.0 - std::exp(-z), tc = std::erf(std::sqrt(z / 2.0));
    t.nz.z.push_back(z);
    t.nz.acentric.push_back(fa);
    t.nz.acentric_theory.push_back(ta);
    t.nz.centric.push_back(fc);
    t.nz.centric_theory.push_back(tc);
    if (t.acentric)
      t.nz.max_deviation_acentric =
          std::max(t.nz.max_deviation_acentric, std::abs(fa - ta));
    if (t.centric)
      t.nz.max_deviation_centric =
          std::max(t.nz.max_deviation_centric, std::abs(fc - tc));
  }

  // The L test as xtriage's (twinning.h, l_test): over the positive acentric
  // E^2, each reflection one neighbour drawn from the offsets (dh, dk, dl),
  // each even, 2 <= |dh|+|dk|+|dl| < 8; none where it was not measured.
  std::map<Miller, std::size_t> acentric_at; //: unique index -> position in u
  for (std::size_t i = 0; i < u.rows.size(); ++i)
    if (!u.centric[i] && u.e2[i] > 0)
      acentric_at[data.hkl[u.rows[i]]] = i;
  std::vector<Miller> offsets;
  for (std::int32_t a = -8; a <= 8; ++a)
    for (std::int32_t b = -8; b <= 8; ++b)
      for (std::int32_t c = -8; c <= 8; ++c) {
        const std::int32_t s = std::abs(a) + std::abs(b) + std::abs(c);
        if (a % 2 == 0 && b % 2 == 0 && c % 2 == 0 && s >= 2 && s < 8)
          offsets.push_back({a, b, c});
      }
  // Each reflection's neighbours tried in an order drawn afresh, to the first
  // measured -- xtriage draws one and gives up where it is missing -- and
  // none the lattice's own symmetry relates to it: such pairs are not
  // independent, equal where a symmetry was missed (L 0, as if twinned).
  const std::vector<Rotation> lattice =
      lattice_symmetry(cell, options.max_delta);
  std::mt19937 random(options.seed);
  double sum_abs = 0, sum_sq = 0;
  for (const auto &[h, i] : acentric_at) {
    std::shuffle(offsets.begin(), offsets.end(), random);
    auto found = acentric_at.end();
    for (const Miller &dh : offsets) {
      const Miller mate = group.unique(add(h, dh));
      bool related = mate == h;
      for (const Rotation &r : lattice)
        related = related || group.unique(mxi::apply(r, h)) == mate;
      if (related)
        continue;
      found = acentric_at.find(mate);
      if (found != acentric_at.end())
        break;
    }
    if (found == acentric_at.end())
      continue;
    const double p = u.e2[i], q = u.e2[found->second];
    const double l = (p - q) / (p + q);
    sum_abs += std::abs(l);
    sum_sq += l * l;
    ++t.l_test.pairs;
  }
  if (t.l_test.pairs) {
    const double n = static_cast<double>(t.l_test.pairs);
    t.l_test.mean_abs_l = sum_abs / n;
    t.l_test.mean_l2 = sum_sq / n;
    // xtriage's (twin_results_interpretation.compute_maha_l): from the mean
    // and inverse covariance of some 1200 PDB data sets.
    const double x = t.l_test.mean_abs_l - 0.48513455414,
                 y = t.l_test.mean_l2 - 0.316418789809;
    t.l_test.z_score = std::sqrt(std::max(
        0.0, x * x * 679728.0 + y * y * 506232.0 - 2.0 * 583582.0 * x * y));
  }
  t.abnormal = t.l_test.pairs > 0 && t.l_test.z_score >= 3.5 &&
               t.l_test.mean_abs_l < 0.5;
  t.wide = t.l_test.pairs > 0 && t.l_test.z_score >= 3.5 &&
           t.l_test.mean_abs_l >= 0.5;

  // Each twin law on the intensities as merged, in range, as xtriage's
  // (twinning.h: twin_r, h_test): R over every reflection with a mate other
  // than itself, centric too, each from both sides; the H test over the
  // acentric, the strongest half of its pairs by I1 + I2.
  std::map<Miller, std::size_t> in_range;
  for (std::size_t i = 0; i < u.rows.size(); ++i)
    in_range[data.hkl[u.rows[i]]] = i;
  for (TwinLaw law : twin_laws(cell, group, options.max_delta)) {
    double top = 0, bottom = 0;
    std::vector<std::pair<double, double>> pairs;
    for (const auto &[h, i] : in_range) {
      const Miller mate = group.unique(mxi::apply(law.op, h));
      const auto found = in_range.find(mate);
      if (mate == h || found == in_range.end())
        continue;
      const double p = data.intensity[u.rows[i]],
                   q = data.intensity[u.rows[found->second]];
      top += std::abs(p - q);
      bottom += std::abs(p + q);
      if (!u.centric[i] && !u.centric[found->second])
        pairs.emplace_back(p, q);
    }
    law.r_obs = bottom > 0 ? top / bottom : 0.0;
    law.pairs = pairs.size() / 2; // each from both sides
    std::sort(pairs.begin(), pairs.end(), [](const auto &a, const auto &b) {
      return a.first + a.second > b.first + b.second;
    });
    const std::size_t half = pairs.size() / 2;
    double h_abs = 0;
    std::size_t used = 0;
    for (std::size_t k = 0; k < half; ++k) {
      const double s = pairs[k].first + pairs[k].second;
      if (!(s > 0))
        continue;
      h_abs += std::abs(pairs[k].first - pairs[k].second) / s;
      ++used;
    }
    law.h_alpha = used ? std::clamp(0.5 - h_abs / used, 0.0, 0.5) : 0.0;
    t.twin_laws.push_back(law);
  }

  // The summary, as xtriage's: twinned only where the statistics independent
  // of twin laws are abnormal and a twin law could explain it; a twin law's
  // pairs nearly equal with the statistics normal is a symmetry missing.
  std::ostringstream v;
  const TwinLaw *highest = nullptr;
  for (const TwinLaw &law : t.twin_laws)
    if (!highest || law.h_alpha > highest->h_alpha)
      highest = &law;
  char z[32];
  std::snprintf(z, sizeof z, "%.2f", t.l_test.z_score);
  auto round_to = [](double x, double unit) {
    return std::round(x / unit) * unit;
  };
  const bool missed = highest && highest->h_alpha > 0.4 && !t.abnormal;
  if (missed) {
    v << "The intensities related by " << highest->text
      << " are nearly equal (H-test fraction "
      << round_to(highest->h_alpha, 0.01) << ", R "
      << round_to(highest->r_obs, 0.001)
      << ") while the L test shows no twinning: more likely a symmetry element "
         "missed than a twin. The space group may be too low.";
    if (t.wide)
      v << " (The L test's spread is wider than untwinned data's, Z " << z
        << ", as the zones a higher symmetry makes centric would make it.)";
  } else if (t.abnormal && highest)
    v << "The L test's statistics differ from untwinned data's as twinning "
         "makes them (Z "
      << z << ", <|L|> " << round_to(t.l_test.mean_abs_l, 0.001)
      << "). As twin laws are possible, twinning could be why: the largest "
         "H-test fraction is "
      << round_to(highest->h_alpha, 0.01) << ", for " << highest->text
      << ". Refining against a twin target may be worth trying.";
  else if (t.abnormal)
    v << "The L test's statistics differ from untwinned data's as twinning "
         "makes them (Z "
      << z
      << "), yet no twin law is possible here: something else -- the "
         "data themselves, anisotropy -- is worth looking for.";
  else if (t.wide)
    v << "The L test's spread is wider than untwinned data's (Z " << z
      << ", <|L|> " << round_to(t.l_test.mean_abs_l, 0.001)
      << "): not twinning's signature. Translational NCS can do it.";
  else if (highest)
    v << "The L test's statistics are as expected (Z " << z
      << "): no twinning is suspected. Twin laws are possible, though, and "
         "a small twin fraction is hard to tell from none.";
  else
    v << "The L test's statistics are as expected (Z " << z
      << ") and no twin law is possible: no twinning.";
  t.verdict = v.str();
  return t;
}

} // namespace mxi
