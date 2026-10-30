#include "verify.hh"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <tuple>

#include "predict.hh"

namespace mxi {

namespace {

double d_of(const Mat3 &A, int h, int k, int l) {
  const Vec3 q = A * Vec3{static_cast<double>(h), static_cast<double>(k),
                          static_cast<double>(l)};
  const double n = q.norm();
  return n > 0.0 ? 1.0 / n : 0.0;
}

} // namespace

VerifyResult verify_lattice(const ExperimentList &experiments, std::size_t id,
                            const Table &table, const VerifyOptions &options,
                            const std::vector<std::uint8_t> *mask,
                            std::size_t mask_width) {
  const Experiment &e = experiments[id];
  if (!e.crystal)
    throw std::runtime_error("experiment " + std::to_string(id) +
                             " has no crystal");
  const Mat3 A = e.crystal->A;
  bool invertible = false;
  const Mat3 A_inverse = A.inverse(&invertible);
  if (!invertible)
    throw std::runtime_error("experiment " + std::to_string(id) +
                             ": its crystal's matrix is singular");
  VerifyResult result;
  result.masked = mask != nullptr;

  // The sweep's strong spots: any experiment of the same sweep, indexed or
  // not; rows with no experiment (DIALS's -1) where there is one sweep only.
  const Column &ids = table.at("id");
  const Column &xyz = table.at("xyzobs.px.value");
  const Column *hkl =
      table.has("miller_index") ? &table.at("miller_index") : nullptr;
  const Column *flags = table.has("flags") ? &table.at("flags") : nullptr;
  const Column *panels = table.has("panel") ? &table.at("panel") : nullptr;
  const std::size_t sweep = experiments.sweep(id);
  bool one_sweep = true;
  for (std::size_t i = 0; i < experiments.size(); ++i)
    one_sweep = one_sweep && experiments.sweep(i) == sweep;

  struct Spot {
    int h, k, l;
    double z;
    bool indexed_here; // indexed to this lattice, with these indices
  };
  std::vector<Spot> spots;
  std::vector<double> inverse_d;
  for (std::size_t r = 0; r < table.nrows; ++r) {
    const std::int64_t rid = ids.ints[r];
    const bool in_sweep =
        rid >= 0 ? static_cast<std::size_t>(rid) < experiments.size() &&
                       experiments.sweep(static_cast<std::size_t>(rid)) == sweep
                 : one_sweep;
    if (!in_sweep)
      continue;
    if (flags && (flags->ints[r] & (flag::kIndexed | flag::kStrong)) == 0)
      continue;
    const double x = xyz.reals[r * 3], y = xyz.reals[r * 3 + 1],
                 z = xyz.reals[r * 3 + 2];
    if (!std::isfinite(z) || (x == 0.0 && y == 0.0))
      continue;
    ++result.spots;
    const std::size_t panel =
        panels ? static_cast<std::size_t>(panels->ints[r]) : 0;
    const Vec3 f = A_inverse * reciprocal_lattice_point(e, panel, x, y, z);
    const int h = static_cast<int>(std::lround(f.x));
    const int k = static_cast<int>(std::lround(f.y));
    const int l = static_cast<int>(std::lround(f.z));
    if ((h == 0 && k == 0 && l == 0) ||
        std::max({std::fabs(f.x - h), std::fabs(f.y - k),
                  std::fabs(f.z - l)}) >= options.tolerance)
      continue;
    bool here = false;
    if (rid == static_cast<std::int64_t>(id) && hkl) {
      here = hkl->ints[r * 3] == h && hkl->ints[r * 3 + 1] == k &&
             hkl->ints[r * 3 + 2] == l;
      if (here) {
        ++result.indexed;
        const double d = d_of(A, h, k, l);
        if (d > 0.0)
          inverse_d.push_back(1.0 / d);
      }
    }
    spots.push_back({h, k, l, z, here});
  }
  if (inverse_d.empty() && options.d_min <= 0.0)
    return result;

  // The limit: given, or the one 99.5 per cent of the indexed spots reach.
  if (options.d_min > 0.0) {
    result.d_min = options.d_min;
    result.d_min_from_spots = false;
  } else {
    std::sort(inverse_d.begin(), inverse_d.end());
    const std::size_t at = static_cast<std::size_t>(
        std::floor(0.995 * static_cast<double>(inverse_d.size() - 1)));
    result.d_min = 1.0 / inverse_d[at];
  }

  PredictOptions po;
  po.d_min = result.d_min;
  const std::vector<Prediction> predicted = predict(e, po);

  // Shells of equal volume in reciprocal space, to the limit from the lowest
  // resolution predicted.
  double d_max = 0.0;
  for (const Prediction &p : predicted)
    d_max = std::max(d_max, d_of(A, p.h, p.k, p.l));
  const std::size_t n = std::max<std::size_t>(1, options.shells);
  const double lo = std::pow(1.0 / std::max(d_max, result.d_min), 3);
  const double hi = std::pow(1.0 / result.d_min, 3);
  result.shells.resize(n);
  for (std::size_t s = 0; s < n; ++s) {
    const double a =
        lo + (hi - lo) * static_cast<double>(s) / static_cast<double>(n);
    const double b =
        lo + (hi - lo) * static_cast<double>(s + 1) / static_cast<double>(n);
    result.shells[s].d_max = a > 0.0 ? 1.0 / std::cbrt(a) : d_max;
    result.shells[s].d_min = 1.0 / std::cbrt(b);
  }
  result.all.d_max = result.shells.front().d_max;
  result.all.d_min = result.d_min;

  // Each spot to the nearest in rotation of its reflection's crossings.
  std::map<std::tuple<int, int, int>, std::vector<std::size_t>> crossings;
  for (std::size_t i = 0; i < predicted.size(); ++i)
    crossings[{predicted[i].h, predicted[i].k, predicted[i].l}].push_back(i);
  std::vector<char> seen(predicted.size(), 0);
  std::vector<char> seen_indexed(predicted.size(), 0);
  for (const Spot &s : spots) {
    const auto it = crossings.find({s.h, s.k, s.l});
    if (it == crossings.end())
      continue; // beyond the limit, or no crossing in the scan
    std::size_t best = it->second.front();
    for (std::size_t i : it->second)
      if (std::fabs(predicted[i].z - s.z) < std::fabs(predicted[best].z - s.z))
        best = i;
    seen[best] = 1;
    if (s.indexed_here)
      seen_indexed[best] = 1;
  }

  for (std::size_t i = 0; i < predicted.size(); ++i) {
    const Prediction &p = predicted[i];
    const double d = d_of(A, p.h, p.k, p.l);
    if (!(d >= result.d_min))
      continue;
    const double v = std::pow(1.0 / d, 3);
    std::size_t s =
        hi > lo ? static_cast<std::size_t>(
                      std::floor((v - lo) / (hi - lo) * static_cast<double>(n)))
                : 0;
    s = std::min(s, n - 1);
    VerifyShell &shell = result.shells[s];
    ++shell.predicted;
    ++result.all.predicted;
    if (mask && p.panel == 0) {
      const long x = static_cast<long>(std::floor(p.px_fast));
      const long y = static_cast<long>(std::floor(p.px_slow));
      const std::size_t at = static_cast<std::size_t>(y) * mask_width +
                             static_cast<std::size_t>(x);
      if (x < 0 || y < 0 || static_cast<std::size_t>(x) >= mask_width ||
          at >= mask->size() || (*mask)[at] != 0) {
        ++shell.unobservable;
        ++result.all.unobservable;
        continue;
      }
    }
    if (seen[i]) {
      ++shell.observed;
      ++result.all.observed;
    }
    if (seen_indexed[i]) {
      ++shell.indexed;
      ++result.all.indexed;
    }
  }
  return result;
}

} // namespace mxi
