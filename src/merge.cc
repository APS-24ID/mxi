#include "merge.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <unordered_map>

#include "french_wilson_tables.hh"

namespace mxi {

MergedSet merge_equivalents(const std::vector<Miller> &hkl,
                            const std::vector<double> &intensity,
                            const std::vector<double> &variance,
                            const gemmi::SpaceGroup &group, bool anomalous) {
  const gemmi::GroupOps ops = group.operations();
  const gemmi::ReciprocalAsu asu(&group);
  struct Sum {
    Miller h;
    std::int8_t hemisphere = 0;
    bool centric = false;
    double w = 0.0, wi = 0.0;
    std::int32_t n = 0;
  };
  std::vector<Sum> sums;
  std::unordered_map<std::uint64_t, std::size_t> slot;
  for (std::size_t i = 0; i < hkl.size(); ++i) {
    if (!(variance[i] > 0.0))
      continue;
    const gemmi::Op::Miller m{hkl[i][0], hkl[i][1], hkl[i][2]};
    const auto [h, isym] = asu.to_asu(m, ops);
    const bool centric = ops.is_reflection_centric(m);
    const std::int8_t hemisphere =
        anomalous && !centric ? (isym % 2 == 1 ? 1 : -1) : 0;
    const std::uint64_t key = (static_cast<std::uint64_t>(h[0] + 1024) << 43) |
                              (static_cast<std::uint64_t>(h[1] + 1024) << 22) |
                              (static_cast<std::uint64_t>(h[2] + 1024) << 1) |
                              (hemisphere < 0 ? 1u : 0u);
    auto it = slot.find(key);
    if (it == slot.end()) {
      it = slot.emplace(key, sums.size()).first;
      Sum s;
      s.h = Miller{h[0], h[1], h[2]};
      s.hemisphere = hemisphere;
      s.centric = centric;
      sums.push_back(s);
    }
    Sum &s = sums[it->second];
    const double w = 1.0 / variance[i];
    s.w += w;
    s.wi += w * intensity[i];
    ++s.n;
  }
  std::sort(sums.begin(), sums.end(), [](const Sum &a, const Sum &b) {
    return a.h != b.h ? a.h < b.h : a.hemisphere > b.hemisphere;
  });
  MergedSet out;
  for (const Sum &s : sums) {
    out.hkl.push_back(s.h);
    out.value.push_back(s.wi / s.w);
    out.sigma.push_back(1.0 / std::sqrt(s.w));
    out.n.push_back(s.n);
    out.centric.push_back(s.centric);
    out.hemisphere.push_back(s.hemisphere);
  }
  return out;
}

namespace {

// np.interp's sense: x within xs, linear between; outside, the end values.
double interpolate(const std::vector<double> &xs, const std::vector<double> &ys,
                   double x) {
  if (x <= xs.front())
    return ys.front();
  if (x >= xs.back())
    return ys.back();
  const auto it = std::upper_bound(xs.begin(), xs.end(), x);
  const std::size_t j = static_cast<std::size_t>(it - xs.begin());
  const double x0 = xs[j - 1], x1 = xs[j];
  if (x1 == x0)
    return ys[j - 1];
  return ys[j - 1] + (ys[j] - ys[j - 1]) * (x - x0) / (x1 - x0);
}

// A table on x = i / 10 - 4, as scipy's interp1d reads it.
template <std::size_t N> double table(const double (&t)[N], double h) {
  const double u = (h + 4.0) * 10.0;
  std::size_t i = static_cast<std::size_t>(std::floor(u));
  if (i >= N - 1)
    i = N - 2;
  const double f = u - static_cast<double>(i);
  return t[i] + (t[i + 1] - t[i]) * f;
}

double median(std::vector<double> v) {
  if (v.empty())
    return std::numeric_limits<double>::quiet_NaN();
  const std::size_t k = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(k),
                   v.end());
  const double upper = v[k];
  if (v.size() % 2 == 1)
    return upper;
  const double lower =
      *std::max_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(k));
  return 0.5 * (lower + upper);
}

} // namespace

FrenchWilson french_wilson(const std::vector<double> &intensity,
                           const std::vector<double> &sigma,
                           const std::vector<bool> &centric,
                           const std::vector<double> &d_star_cubed,
                           std::size_t max_bins) {
  namespace t = french_wilson_tables;
  const std::size_t n = intensity.size();
  FrenchWilson out;
  out.F.assign(n, 0.0);
  out.sigF.assign(n, 0.0);
  out.J.assign(n, 0.0);
  out.sigJ.assign(n, 0.0);
  out.valid.assign(n, false);
  if (n == 0)
    return out;
  // determine_binning_counting_sorted: edges at np.interp(linspace(0, n,
  // n_bins + 1), arange(n), sort(x)).
  std::vector<double> sorted = d_star_cubed;
  std::sort(sorted.begin(), sorted.end());
  const std::size_t per_bin = std::max<std::size_t>(n / max_bins, 20);
  const std::size_t n_bins = std::max<std::size_t>(n / per_bin, 1);
  std::vector<double> edges(n_bins + 1);
  for (std::size_t b = 0; b <= n_bins; ++b) {
    const double at = static_cast<double>(n) * static_cast<double>(b) /
                      static_cast<double>(n_bins);
    const std::size_t lo =
        std::min<std::size_t>(static_cast<std::size_t>(std::floor(at)), n - 1);
    const std::size_t hi = std::min<std::size_t>(lo + 1, n - 1);
    const double f = at - static_cast<double>(lo);
    edges[b] = lo == n - 1 ? sorted[n - 1]
                           : sorted[lo] + (sorted[hi] - sorted[lo]) * f;
  }
  // scipy's binned_statistic: [e_i, e_i+1), the last bin closed.
  std::vector<std::vector<double>> in_bin(n_bins);
  for (std::size_t i = 0; i < n; ++i) {
    const double x = d_star_cubed[i];
    if (x < edges.front() || x > edges.back())
      continue;
    std::size_t b = static_cast<std::size_t>(
        std::upper_bound(edges.begin(), edges.end(), x) - edges.begin());
    b = b == 0 ? 0 : b - 1;
    if (b >= n_bins)
      b = n_bins - 1;
    in_bin[b].push_back(intensity[i]);
  }
  std::vector<double> centres, statistic;
  for (std::size_t b = 0; b < n_bins; ++b) {
    centres.push_back(0.5 * (edges[b] + edges[b + 1]));
    statistic.push_back(median(in_bin[b]) / std::numbers::ln2);
  }
  for (std::size_t i = 0; i < n; ++i) {
    const double I = intensity[i], s = sigma[i];
    if (!(s > 0.0))
      continue;
    const double expected = interpolate(centres, statistic, d_star_cubed[i]);
    if (expected == 0.0 || !std::isfinite(expected))
      continue;
    const double i_over_s = I / s;
    if (!centric[i]) {
      const double h = i_over_s - std::fabs(s / expected);
      if (!(i_over_s >= -3.7 && h >= -4.0))
        continue;
      if (h < 3.0) {
        out.J[i] = table(t::ac_zj, h) * s;
        out.sigJ[i] = table(t::ac_zj_sd, h) * s;
        out.F[i] = table(t::ac_zf, h) * std::sqrt(s);
        out.sigF[i] = table(t::ac_zf_sd, h) * std::sqrt(s);
      } else {
        out.J[i] = h * s;
        out.sigJ[i] = s;
        out.F[i] = std::sqrt(out.J[i]);
        out.sigF[i] = 0.5 * s / out.F[i];
      }
    } else {
      const double h = i_over_s - std::fabs(s / (2.0 * expected));
      if (!(i_over_s >= -3.7 && h >= -4.0))
        continue;
      if (h < 4.0) {
        out.J[i] = table(t::c_zj, h) * s;
        out.sigJ[i] = table(t::c_zj_sd, h) * s;
        out.F[i] = table(t::c_zf, h) * std::sqrt(s);
        out.sigF[i] = table(t::c_zf_sd, h) * std::sqrt(s);
      } else {
        // French and Wilson's expansion with the x^6 term, as dials'.
        const double h2 = 1.0 / (h * h), h4 = h2 * h2, h6 = h2 * h4;
        const double f =
            std::sqrt(h) * (1.0 - (3.0 / 8.0) * h2 - (87.0 / 128.0) * h4 -
                            (2889.0 / 1024.0) * h6);
        const double sf = std::sqrt(
            h * ((1.0 / 4.0) * h2 + (15.0 / 32.0) * h4 + (273.0 / 128.0) * h6));
        out.J[i] = h * s * (1.0 - 0.5 * h2 - 0.75 * h4 - 3.0 * h6);
        out.sigJ[i] = 2.0 * s * f * sf;
        out.F[i] = f * std::sqrt(s);
        out.sigF[i] = sf * std::sqrt(s);
      }
    }
    out.valid[i] = true;
  }
  return out;
}

std::int32_t free_flag(const Miller &hkl, const std::vector<Rotation> &lattice,
                       std::int32_t groups) {
  Miller best = hkl;
  bool first = true;
  for (const Rotation &r : lattice) {
    const Miller m = apply(r, hkl);
    for (const Miller &c : {m, Miller{-m[0], -m[1], -m[2]}})
      if (first || c < best) {
        best = c;
        first = false;
      }
  }
  // SplitMix64 of the canonical index: the same flag wherever it is found.
  std::uint64_t z = (static_cast<std::uint64_t>(best[0] + 1024) << 42) |
                    (static_cast<std::uint64_t>(best[1] + 1024) << 21) |
                    static_cast<std::uint64_t>(best[2] + 1024);
  z += 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  z ^= z >> 31;
  return static_cast<std::int32_t>(z % static_cast<std::uint64_t>(groups));
}

} // namespace mxi
