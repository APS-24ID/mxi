// Intensity-based clustering of isomorphous data sets: see correlation.hh.

#include "correlation.hh"

#include <algorithm>
#include <cmath>
#include <map>

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

} // namespace mxi
