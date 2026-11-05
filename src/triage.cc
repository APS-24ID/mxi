#include "triage.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>

namespace mxi {

namespace {

Rotation multiply(const Rotation &a, const Rotation &b) {
  Rotation c{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      int s = 0;
      for (int k = 0; k < 3; ++k)
        s += a[i * 3 + k] * b[k * 3 + j];
      c[i * 3 + j] = s;
    }
  return c;
}

int determinant(const Rotation &r) {
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

// The reflections used: in [d_min, d_max], each with its E^2 and whether it is
// centric, by index.
struct Used {
  std::vector<std::size_t> rows;
  std::map<Miller, std::size_t> at; //: unique index -> position in rows
  std::vector<double> e2;
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
    u.at[data.hkl[r]] = u.rows.size();
    u.rows.push_back(r);
    u.centric.push_back(group.centric(data.hkl[r]));
    s2.push_back(1.0 / (d * d));
    scaled.push_back(data.intensity[r] / epsilon(group, data.hkl[r]));
  }
  // <I/eps> in shells of equal count by resolution, about 400 a shell, 1 to
  // 20 of them; E^2 each one's I/eps over its shell's mean.
  const std::size_t n = u.rows.size();
  u.e2.assign(n, 0.0);
  if (n == 0)
    return u;
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](std::size_t a, std::size_t b) { return s2[a] < s2[b]; });
  const std::size_t shells = std::clamp<std::size_t>(n / 400, 1, 20);
  for (std::size_t k = 0; k < shells; ++k) {
    const std::size_t lo = k * n / shells, hi = (k + 1) * n / shells;
    double mean = 0.0;
    for (std::size_t i = lo; i < hi; ++i)
      mean += scaled[order[i]];
    mean /= static_cast<double>(std::max<std::size_t>(hi - lo, 1));
    for (std::size_t i = lo; i < hi; ++i)
      u.e2[order[i]] = mean > 0.0 ? scaled[order[i]] / mean : 0.0;
  }
  return u;
}

Moments moments(const Used &u, bool centric) {
  Moments m;
  double e2 = 0, e4 = 0, e = 0, dev = 0;
  for (std::size_t i = 0; i < u.e2.size(); ++i) {
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

} // namespace

int epsilon(const SpaceGroup &group, const Miller &hkl) {
  int n = 0;
  for (const Rotation &r : point_group(group))
    if (mxi::apply(r, hkl) == hkl)
      ++n;
  return std::max(n, 1);
}

std::string twin_law_text(const Rotation &r) {
  // h' = R^T h: component i of h' is sum_j R[j][i] h_j.
  const char *names = "hkl";
  std::string out;
  for (int i = 0; i < 3; ++i) {
    std::string term;
    for (int j = 0; j < 3; ++j) {
      const int c = r[j * 3 + i];
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
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        double s = 0.0;
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            s += pick[a * 3 + i] * g[a * 3 + b] * pick[b * 3 + j];
        moved[i * 3 + j] = s;
      }
    double diff = 0.0;
    for (int k = 0; k < 9; ++k)
      diff += (moved[k] - g[k]) * (moved[k] - g[k]);
    law.r_metric = norm > 0 ? diff / norm : 0.0;
    law.merohedral = law.r_metric < 1e-8;
    out.push_back(law);
  }
  return out;
}

double twin_test_d_min(const TriageData &data, const UnitCell &cell,
                       const SpaceGroup &group, const TriageOptions &o) {
  double data_d_min = 1e9;
  for (const Miller &m : data.hkl)
    data_d_min = std::min(data_d_min, d_spacing(cell, m));
  if (data.hkl.empty() || data_d_min >= o.d_max)
    return data_d_min;
  // The reflections possible in the asymmetric unit, from d_max to the data's
  // limit, and of those, how many have I/sigma at least isigi_cut: shells of
  // equal volume in 1/d^3, ten.
  const int shells = 10;
  const double lo = 1.0 / (o.d_max * o.d_max * o.d_max),
               hi = 1.0 / (data_d_min * data_d_min * data_d_min);
  auto shell_of = [&](double d) {
    const double x = 1.0 / (d * d * d);
    return std::clamp(static_cast<int>((x - lo) / (hi - lo) * shells), 0,
                      shells - 1);
  };
  std::vector<double> possible(shells, 0.0), strong(shells, 0.0);
  const int hm = static_cast<int>(cell.a / data_d_min) + 1,
            km = static_cast<int>(cell.b / data_d_min) + 1,
            lm = static_cast<int>(cell.c / data_d_min) + 1;
  std::set<Miller> counted;
  for (int h = -hm; h <= hm; ++h)
    for (int k = -km; k <= km; ++k)
      for (int l = -lm; l <= lm; ++l) {
        const Miller m{h, k, l};
        if (h == 0 && k == 0 && l == 0)
          continue;
        const double d = d_spacing(cell, m);
        if (d < data_d_min || d > o.d_max || group.absent(m))
          continue;
        if (counted.insert(group.unique(m)).second)
          possible[shell_of(d)] += 1.0;
      }
  for (std::size_t r = 0; r < data.hkl.size(); ++r) {
    const double d = d_spacing(cell, data.hkl[r]);
    if (d > o.d_max || d < data_d_min || !(data.sigma[r] > 0))
      continue;
    if (data.intensity[r] / data.sigma[r] >= o.isigi_cut)
      strong[shell_of(d)] += 1.0;
  }
  double d_min = o.d_max;
  for (int s = 0; s < shells; ++s) {
    if (possible[s] == 0 || strong[s] / possible[s] < o.completeness)
      break;
    d_min = std::cbrt(1.0 / (lo + (hi - lo) * (s + 1) / shells));
  }
  return std::max(std::min(d_min, o.coarsest_d_min), data_d_min);
}

Triage triage(const TriageData &data, const UnitCell &cell,
              const SpaceGroup &group, const TriageOptions &options) {
  Triage t;
  t.d_max = options.d_max;
  t.d_min = options.d_min > 0 ? options.d_min
                              : twin_test_d_min(data, cell, group, options);
  const Used u = normalised(data, cell, group, t.d_min, t.d_max);
  t.used = u.rows.size();
  for (bool c : u.centric)
    (c ? t.centric : t.acentric) += 1;
  t.acentric_moments = moments(u, false);
  t.centric_moments = moments(u, true);

  // N(z), z = 0.1 .. 1.0.
  for (int k = 1; k <= 10; ++k) {
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

  // The L test: each acentric reflection and a neighbour (2h', 2k', 2l'),
  // 1 <= |h'|+|k'|+|l'| <= 4, tried in an order drawn afresh for each.
  std::vector<Miller> offsets;
  for (int a = -4; a <= 4; ++a)
    for (int b = -4; b <= 4; ++b)
      for (int c = -4; c <= 4; ++c) {
        const int s = std::abs(a) + std::abs(b) + std::abs(c);
        if (s >= 1 && s <= 4)
          offsets.push_back({2 * a, 2 * b, 2 * c});
      }
  std::mt19937 random(options.seed);
  double sum_abs = 0, sum_sq = 0;
  // Not a neighbour the lattice's own symmetry relates to it: equal where a
  // symmetry was missed, or near it where twinned, L = 0 then, and the test
  // had called P 4 2 2 data taken as P 4 abnormal (Z 7.8) -- L is to compare
  // unrelated reflections, the mixture twinning makes of each still showing.
  const std::vector<Rotation> lattice =
      lattice_symmetry(cell, options.max_delta);
  for (std::size_t i = 0; i < u.rows.size(); ++i) {
    if (u.centric[i])
      continue;
    const Miller &h = data.hkl[u.rows[i]];
    std::shuffle(offsets.begin(), offsets.end(), random);
    for (const Miller &dh : offsets) {
      const Miller mate = group.unique(add(h, dh));
      const auto found = u.at.find(mate);
      if (mate == h || found == u.at.end() || u.centric[found->second])
        continue;
      bool related = false;
      for (const Rotation &r : lattice)
        related = related || group.unique(mxi::apply(r, h)) == mate;
      if (related)
        continue;
      // E^2, not I: a neighbour two indices off is a little further out, and
      // the fall-off alone had made the two differ -- <|L|> up from 1/2,
      // untwinned data abnormal to the Z score over 30 000 pairs.
      const double i1 = u.e2[i], i2 = u.e2[found->second];
      if (!(i1 + i2 > 0))
        break;
      const double l = std::clamp((i1 - i2) / (i1 + i2), -1.0, 1.0);
      sum_abs += std::abs(l);
      sum_sq += l * l;
      ++t.l_test.pairs;
      break;
    }
  }
  if (t.l_test.pairs) {
    const double n = static_cast<double>(t.l_test.pairs);
    t.l_test.mean_abs_l = sum_abs / n;
    t.l_test.mean_l2 = sum_sq / n;
    // Mahalanobis distance from (1/2, 1/3), L uniform on [-1, 1]: var |L| =
    // 1/12, var L^2 = 4/45, cov = 1/12; of the means, over n.
    const double a = 1.0 / 12.0 / n, b = 1.0 / 12.0 / n, c = 4.0 / 45.0 / n;
    const double det = a * c - b * b;
    const double x = t.l_test.mean_abs_l - 0.5,
                 y = t.l_test.mean_l2 - 1.0 / 3.0;
    t.l_test.z_score =
        std::sqrt(std::max(0.0, (c * x * x - 2 * b * x * y + a * y * y) / det));
  }
  t.abnormal = t.l_test.pairs > 0 && t.l_test.z_score > 3.5;
  t.twin_like = t.abnormal && t.l_test.mean_abs_l < 0.5;

  // Each twin law: its pairs of acentric reflections in range, each once.
  for (TwinLaw law : twin_laws(cell, group, options.max_delta)) {
    std::vector<std::pair<double, double>> pairs;
    for (std::size_t i = 0; i < u.rows.size(); ++i) {
      if (u.centric[i])
        continue;
      const Miller &h = data.hkl[u.rows[i]];
      const Miller mate = group.unique(mxi::apply(law.op, h));
      const auto found = u.at.find(mate);
      if (!(h < mate) || found == u.at.end() || u.centric[found->second])
        continue;
      pairs.emplace_back(data.intensity[u.rows[i]],
                         data.intensity[u.rows[found->second]]);
    }
    law.pairs = pairs.size();
    double num = 0, den = 0;
    for (const auto &[a, b] : pairs) {
      num += std::abs(a - b);
      den += a + b;
    }
    law.r_obs = den > 0 ? num / den : 0.0;
    std::sort(pairs.begin(), pairs.end(), [](const auto &p, const auto &q) {
      return p.first + p.second > q.first + q.second;
    });
    double h_abs = 0;
    std::size_t used = 0;
    for (std::size_t k = 0; k < (pairs.size() + 1) / 2; ++k) {
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
  const bool missed = highest && highest->h_alpha > 0.4 && !t.twin_like;
  if (missed) {
    v << "The intensities related by " << highest->text
      << " are nearly equal (H-test fraction "
      << round_to(highest->h_alpha, 0.01) << ", R "
      << round_to(highest->r_obs, 0.001)
      << ") while the L test shows no twinning: more likely a symmetry element "
         "missed than a twin. The space group may be too low.";
    if (t.abnormal)
      v << " (The L test's spread is wider than untwinned data's, Z " << z
        << ", as the zones a higher symmetry makes centric would make it.)";
  } else if (t.twin_like && highest)
    v << "The L test's statistics differ from untwinned data's as twinning "
         "makes them (Z "
      << z << ", <|L|> " << round_to(t.l_test.mean_abs_l, 0.001)
      << "). As twin laws are possible, twinning could be why: the largest "
         "H-test fraction is "
      << round_to(highest->h_alpha, 0.01) << ", for " << highest->text
      << ". Refining against a twin target may be worth trying.";
  else if (t.twin_like)
    v << "The L test's statistics differ from untwinned data's as twinning "
         "makes them (Z "
      << z
      << "), yet no twin law is possible here: something else -- the "
         "data themselves, anisotropy -- is worth looking for.";
  else if (t.abnormal)
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
