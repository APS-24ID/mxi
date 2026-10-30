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

VerifyResult verify_lattice(const Experiment &e, std::size_t id,
                            const Table &table, const VerifyOptions &options,
                            const std::vector<std::uint8_t> *mask,
                            std::size_t mask_width) {
  if (!e.crystal)
    throw std::runtime_error("experiment " + std::to_string(id) +
                             " has no crystal");
  const Mat3 A = e.crystal->A;
  VerifyResult result;
  result.masked = mask != nullptr;

  // The lattice's indexed spots, by their indices: where in the scan each was.
  const Column &ids = table.at("id");
  const Column &hkl = table.at("miller_index");
  const Column &xyz = table.at("xyzobs.px.value");
  const Column *flags = table.has("flags") ? &table.at("flags") : nullptr;
  std::map<std::tuple<int, int, int>, std::vector<double>> seen;
  std::vector<double> inverse_d;
  for (std::size_t r = 0; r < table.nrows; ++r) {
    if (ids.ints[r] != static_cast<std::int64_t>(id))
      continue;
    const int h = static_cast<int>(hkl.ints[r * 3]);
    const int k = static_cast<int>(hkl.ints[r * 3 + 1]);
    const int l = static_cast<int>(hkl.ints[r * 3 + 2]);
    if (h == 0 && k == 0 && l == 0)
      continue;
    if (flags && (flags->ints[r] & (flag::kIndexed | flag::kStrong)) == 0)
      continue;
    const double z = xyz.reals[r * 3 + 2];
    if (!std::isfinite(z) ||
        (xyz.reals[r * 3] == 0.0 && xyz.reals[r * 3 + 1] == 0.0))
      continue;
    seen[{h, k, l}].push_back(z);
    const double d = d_of(A, h, k, l);
    if (d > 0.0)
      inverse_d.push_back(1.0 / d);
  }
  result.indexed = inverse_d.size();
  if (inverse_d.empty())
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

  // Within the tolerance in rotation, in images.
  const double width = std::fabs(e.scan.osc_width);
  const double tolerance =
      width > 0.0 ? options.tolerance_degrees / width : 1.0e300;

  for (const Prediction &p : predicted) {
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
    const auto it = seen.find({p.h, p.k, p.l});
    if (it == seen.end())
      continue;
    for (double z : it->second)
      if (std::fabs(z - p.z) <= tolerance) {
        ++shell.observed;
        ++result.all.observed;
        break;
      }
  }
  return result;
}

} // namespace mxi
