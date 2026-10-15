#include "scale.hh"

#include "parallel.hh"
#include "processing.hh"
#include "timing.hh"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>

namespace mxi {

ScaleData build_scale_data(const ExperimentList &experiments,
                           const Table &reflections, const SpaceGroup &group,
                           const ScaleModelShape &shape,
                           const ScaleDataOptions &options) {
  if (experiments.size() == 0)
    throw std::runtime_error("scaling needs an experiment");
  // Each observation from its own sweep: its id names the experiment whose
  // scan gives its place in the rotation and whose goniometer and beam give
  // its crystal frame. With one sweep every id is 0 and all is as it was.
  const std::size_t sweeps = experiments.size();
  const bool profile = options.use == IntensityChoice::profile;
  const Column &miller = reflections.at("miller_index");
  const Column &flags = reflections.at("flags");
  const Column &dcol = reflections.at("d");
  const Column &cal = reflections.at("xyzcal.px");
  const bool has_lp = reflections.has("lp"), has_qe = reflections.has("qe"),
             has_part = reflections.has("partiality"),
             has_s1 = reflections.has("s1"), has_id = reflections.has("id");
  const std::int64_t wanted =
      profile ? flag::kIntegratedPrf : flag::kIntegratedSum;
  // The columns looked up once, not once a row: eight lookups by name a
  // reflection were most of gathering, which was serial.
  const Column *id_col = has_id ? &reflections.at("id") : nullptr;
  const Column *part_col = has_part ? &reflections.at("partiality") : nullptr;
  const Column *lp_col = has_lp ? &reflections.at("lp") : nullptr;
  const Column *qe_col = has_qe ? &reflections.at("qe") : nullptr;
  const Column &pv = reflections.at("intensity.prf.value");
  const Column &pvar = reflections.at("intensity.prf.variance");
  const Column &sv = reflections.at("intensity.sum.value");
  const Column &svar = reflections.at("intensity.sum.variance");
  // Each variance with no signal, if integration wrote them.
  const bool has_bg = reflections.has("intensity.prf.variance_background") &&
                      reflections.has("intensity.sum.variance_background");
  const Column *pvar_bg =
      has_bg ? &reflections.at("intensity.prf.variance_background") : nullptr;
  const Column *svar_bg =
      has_bg ? &reflections.at("intensity.sum.variance_background") : nullptr;

  // First each row on its own, in parallel: whether it is taken, and its
  // symmetry -- absent, its unique index, which Friedel mate -- which loop over
  // the group's operators. Then, in row order, the groups numbered as they
  // first appear and the observations appended, so that everything comes out as
  // it did in one serial loop.
  const std::size_t rows = reflections.nrows;
  // A reflection split by integration across blocks of images, its pieces
  // sharing a partial_id, made one observation first, as dials.scale does
  // (sum_partial_reflections): each piece's profile-fitted intensity already
  // estimates the whole, so their mean weighted by I^2 / sigma^2; the summed
  // intensities and their variances added; the partialities added. Pieces of
  // partiality under 0.99 with every variance positive, within an experiment;
  // the first piece in the table carries the result, the others are merged
  // away. Scaled apart, each piece divided by its own partiality, they were
  // 8 per cent of the strongest observations on 36 insulin crystals.
  std::vector<double> prf_value(pv.reals.begin(), pv.reals.end()),
      prf_var(pvar.reals.begin(), pvar.reals.end()),
      sum_value(sv.reals.begin(), sv.reals.end()),
      sum_var(svar.reals.begin(), svar.reals.end());
  std::vector<double> partiality(rows, 1.0);
  if (part_col)
    partiality.assign(part_col->reals.begin(), part_col->reals.end());
  std::vector<std::uint8_t> merged_away(rows, 0);
  std::size_t combined = 0;
  if (part_col && id_col && reflections.has("partial_id")) {
    const Column &pid = reflections.at("partial_id");
    std::map<std::pair<std::int64_t, std::int64_t>, std::vector<std::size_t>>
        pieces;
    for (std::size_t i = 0; i < rows; ++i) {
      if ((flags.ints[i] & wanted) == 0 || id_col->ints[i] < 0)
        continue;
      if (!(partiality[i] < 0.99) || !(pvar.reals[i] > 0.0) ||
          !(svar.reals[i] > 0.0))
        continue;
      pieces[{id_col->ints[i], pid.ints[i]}].push_back(i);
    }
    for (const auto &[key, j] : pieces) {
      if (j.size() < 2)
        continue;
      double value = 0.0, variance = 0.0, weight = 0.0;
      double svalue = 0.0, svariance = 0.0, total = 0.0;
      for (std::size_t i : j) {
        const double w = pv.reals[i] * pv.reals[i] / pvar.reals[i];
        value += w * pv.reals[i];
        variance += w * pvar.reals[i];
        weight += w;
        svalue += sv.reals[i];
        svariance += svar.reals[i];
        total += partiality[i];
      }
      const std::size_t first = j[0];
      if (weight > 0.0) {
        prf_value[first] = value / weight;
        prf_var[first] = variance / weight;
      } else {
        double v = 0.0;
        for (std::size_t i : j)
          v += pvar.reals[i];
        prf_value[first] = 0.0;
        prf_var[first] = v;
      }
      sum_value[first] = svalue;
      sum_var[first] = svariance;
      partiality[first] = total;
      for (std::size_t k = 1; k < j.size(); ++k)
        merged_away[j[k]] = 1;
      combined += j.size() - 1;
    }
  }
  const std::vector<double> &chosen_value = profile ? prf_value : sum_value;
  const std::vector<double> &chosen_var = profile ? prf_var : sum_var;
  std::vector<std::uint8_t> taken(rows, 0), plus_of(rows, 0);
  std::vector<Miller> unique_of(rows);
  // Two factors a row: lp / qe for the profile-fitted intensity, which
  // estimates the whole reflection whatever part of it was recorded, and lp /
  // (qe partiality) for the summed one, which counts only the recorded part --
  // as dials.scale (combine_intensities). Dividing both by the partiality made
  // every partial profile-fitted intensity too large by 1 / p.
  std::vector<double> prf_factor_of(rows, 0.0), sum_factor_of(rows, 0.0);
  std::vector<std::size_t> sweep_of(rows, 0);
  // An id naming no experiment is not a reflection to drop quietly: the
  // table and the list do not belong together.
  if (id_col) {
    for (std::size_t i = 0; i < rows; ++i) {
      const std::int64_t id = id_col->ints[i];
      if ((flags.ints[i] & wanted) != 0 &&
          (id < 0 || static_cast<std::size_t>(id) >= sweeps))
        throw std::runtime_error(
            "a reflection's id is " + std::to_string(id) + " and there " +
            (sweeps == 1 ? std::string("is one experiment")
                         : "are " + std::to_string(sweeps) + " experiments") +
            ": the reflections are not these experiments'");
    }
  }
  for_each_index(rows, [&](std::size_t i) {
    if ((flags.ints[i] & wanted) == 0 || merged_away[i])
      return;
    if (id_col)
      sweep_of[i] = static_cast<std::size_t>(id_col->ints[i]);
    const double v = chosen_var[i];
    const double d = dcol.reals[i];
    if (!(v > 0.0) || !std::isfinite(chosen_value[i]) || !(d > 0.0) ||
        !std::isfinite(d) || d < options.d_min ||
        (options.d_max > 0.0 && d > options.d_max))
      return;
    const double part = partiality[i];
    if (!(part >= options.partiality_cutoff))
      return;
    const Miller h{static_cast<int>(miller.ints[i * 3]),
                   static_cast<int>(miller.ints[i * 3 + 1]),
                   static_cast<int>(miller.ints[i * 3 + 2])};
    if ((h[0] == 0 && h[1] == 0 && h[2] == 0) || group.absent(h))
      return;
    double factor = 1.0;
    if (lp_col)
      factor *= lp_col->reals[i];
    if (qe_col && qe_col->reals[i] > 0.0)
      factor /= qe_col->reals[i];
    taken[i] = 1;
    prf_factor_of[i] = factor;
    sum_factor_of[i] = factor / part;
    unique_of[i] = group.unique(h);
    plus_of[i] = group.friedel_plus(h) ? 1 : 0;
  });
  std::size_t count = 0;
  for (std::uint8_t t : taken)
    count += t;

  ScaleData data;
  data.partials_combined = combined;
  data.intensity.reserve(count);
  data.variance.reserve(count);
  data.variance_before.reserve(count);
  data.scale_term.reserve(count);
  data.prf.reserve(count);
  data.prf_variance.reserve(count);
  data.sum.reserve(count);
  data.sum_variance.reserve(count);
  data.observation.reserve(count);
  data.group.reserve(count);
  data.row.reserve(count);
  data.d.reserve(count);
  data.outlier.reserve(count);
  // A group's number is the order it first appears in; a hash from the index
  // finds it, where an ordered map searched a tree.
  const auto key = [](const Miller &m) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(m[0]) &
                                       0x1FFFFF)
            << 42) |
           (static_cast<std::uint64_t>(static_cast<std::uint32_t>(m[1]) &
                                       0x1FFFFF)
            << 21) |
           (static_cast<std::uint64_t>(static_cast<std::uint32_t>(m[2]) &
                                       0x1FFFFF));
  };
  std::unordered_map<std::uint64_t, std::size_t> groups, pairs;
  groups.reserve(count / 2 + 16);
  if (options.anomalous) {
    pairs.reserve(count / 2 + 16);
    data.pair.reserve(count);
  }
  for (std::size_t i = 0; i < rows; ++i) {
    if (!taken[i])
      continue;
    const Miller &u = unique_of[i];
    // With anomalous an acentric reflection's I(-) is a group of its own: the
    // top bit, which the index's 63 do not use, says which mate. A centric
    // reflection's mates are never apart, since friedel_plus calls every
    // centric reflection I(+): they are equal by symmetry.
    std::uint64_t group_key = key(u);
    if (options.anomalous && !plus_of[i])
      group_key |= std::uint64_t{1} << 63;
    const auto [it, fresh] = groups.emplace(group_key, data.unique.size());
    if (fresh) {
      data.unique.push_back(u);
      data.centric.push_back(group.centric(u));
    }
    if (options.anomalous) {
      const auto [pt, pair_fresh] =
          pairs.emplace(key(u), data.pair_unique.size());
      if (pair_fresh) {
        data.pair_unique.push_back(u);
        data.pair_centric.push_back(group.centric(u));
      }
      data.pair.push_back(pt->second);
    }
    const double pf = prf_factor_of[i], sf = sum_factor_of[i];
    const double factor = profile ? pf : sf;
    const double v = chosen_var[i];
    const double d = dcol.reals[i];
    ScaleObservation o;
    const double z = cal.reals[i * 3 + 2];
    const Experiment &e = experiments[sweep_of[i]];
    o.rotation = e.scan.num_images() > 0 ? e.scan.fraction(z) : 0.0;
    o.sweep = sweep_of[i];
    o.time = o.rotation;
    o.inv_2d2 = 1.0 / (2.0 * d * d);
    data.intensity.push_back(chosen_value[i] * factor);
    data.variance.push_back(v * factor * factor);
    data.variance_before.push_back(v * factor * factor);
    data.scale_term.push_back(0.0);
    const bool summed =
        (flags.ints[i] & flag::kIntegratedSum) != 0 && sum_var[i] > 0.0;
    data.prf.push_back(prf_value[i] * pf);
    data.prf_variance.push_back(prf_var[i] * pf * pf);
    data.sum.push_back(sum_value[i] * sf);
    data.sum_variance.push_back(sum_var[i] * sf * sf);
    data.factor.push_back(factor);
    if (has_bg) {
      const double pb = pvar_bg->reals[i] * pf * pf;
      const double sb = svar_bg->reals[i] * sf * sf;
      data.prf_variance_bg.push_back(pb);
      data.sum_variance_bg.push_back(sb);
      data.variance_background.push_back(profile ? pb : sb);
    }
    data.has_sum.push_back(summed);
    data.plus.push_back(plus_of[i] != 0);
    data.observation.push_back(std::move(o));
    data.group.push_back(it->second);
    data.row.push_back(i);
    data.d.push_back(d);
    data.outlier.push_back(false);
  }
  data.has_background = has_bg;
  if (shape.lmax > 0 && has_s1) {
    const Column &s1 = reflections.at("s1");
    for_each_index(data.size(), [&](std::size_t k) {
      const std::size_t i = data.row[k];
      const double z = cal.reals[i * 3 + 2];
      const Experiment &e = experiments[data.observation[k].sweep];
      // The crystal frame: undo the goniometer's rotation at this angle.
      const Mat3 r = e.goniometer.setting *
                     rotation(e.goniometer.axis, e.scan.phi_from_z(z)) *
                     e.goniometer.fixed;
      const Mat3 back = r.transpose();
      const Vec3 s1c = (back * Vec3{s1.reals[i * 3], s1.reals[i * 3 + 1],
                                    s1.reals[i * 3 + 2]})
                           .normalized();
      const Vec3 s0c =
          (back * e.beam.direction).normalized(); // the reverse incident beam
      std::vector<double> y1, y0;
      real_spherical_harmonics(shape.lmax, s1c, &y1);
      real_spherical_harmonics(shape.lmax, s0c, &y0);
      std::vector<double> &out = data.observation[k].absorption;
      out.resize(y1.size());
      for (std::size_t q = 0; q < y1.size(); ++q)
        out[q] = 0.5 * (y1[q] + y0[q]);
    });
  }
  return data;
}

std::vector<double> inverse_scales(const ScaleModel &model,
                                   const ScaleData &data) {
  std::vector<double> g(data.size());
  for_each_index(data.size(), [&](std::size_t i) {
    g[i] = model.inverse_scale(data.observation[i]);
  });
  return g;
}

std::vector<double> merged_intensities(const ScaleData &data,
                                       const std::vector<double> &g) {
  std::vector<double> num(data.unique.size(), 0.0),
      den(data.unique.size(), 0.0);
  for (std::size_t i = 0; i < data.size(); ++i) {
    if (data.outlier[i])
      continue;
    const double w = 1.0 / data.variance[i];
    num[data.group[i]] += w * g[i] * data.intensity[i];
    den[data.group[i]] += w * g[i] * g[i];
  }
  for (std::size_t h = 0; h < num.size(); ++h)
    num[h] = den[h] > 0.0 ? num[h] / den[h] : 0.0;
  return num;
}

std::vector<std::size_t> select_for_fitting(const ScaleData &data,
                                            std::size_t min_groups,
                                            std::size_t min_observations,
                                            unsigned seed) {
  std::vector<std::vector<std::size_t>> members(data.unique.size());
  for (std::size_t i = 0; i < data.size(); ++i)
    if (!data.outlier[i])
      members[data.group[i]].push_back(i);
  std::vector<std::size_t> order(members.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::mt19937 rng(seed);
  std::shuffle(order.begin(), order.end(), rng);
  std::vector<std::size_t> out;
  std::size_t groups = 0;
  for (std::size_t h : order) {
    if (groups >= min_groups && out.size() >= min_observations)
      break;
    if (members[h].size() < 2)
      continue; // a single observation says nothing about scale
    out.insert(out.end(), members[h].begin(), members[h].end());
    ++groups;
  }
  std::sort(out.begin(), out.end());
  return out;
}

namespace {

//: The target over `use`, with <I_h> from the same observations, plus the
//: restraints.
// The restraints' share of the target: every sweep's decay and harmonic terms,
// each pulled towards zero -- added to `r` one by one, in the order they always
// were, so that one sweep's target is the same to the last bit.
double restraint_of(const ScaleModel &model, const ScaleFitOptions &options,
                    double r = 0.0) {
  for (std::size_t w = 0; w < model.sweeps(); ++w) {
    const ScaleModelShape &s = model.shape(w);
    for (std::size_t k = 0; k < s.decay_points; ++k)
      r += options.decay_restraint *
           model.parameters[model.first_decay(w) + k] *
           model.parameters[model.first_decay(w) + k];
    for (std::size_t k = 0; k < harmonic_count(s.lmax); ++k)
      r += options.absorption_restraint *
           model.parameters[model.first_absorption(w) + k] *
           model.parameters[model.first_absorption(w) + k];
  }
  return r;
}

double target(const ScaleModel &model, const ScaleData &data,
              const std::vector<std::vector<std::size_t>> &members,
              const std::vector<std::size_t> &use,
              const ScaleFitOptions &options, std::vector<double> *g_out,
              std::vector<double> *merged_out) {
  // In parallel over groups, each group's sums in its own members' order --
  // the order of `use` within it -- and the total summed in `use` order: the
  // same arithmetic as one thread, on any number.
  std::vector<double> g(data.size(), 0.0), merged(data.unique.size(), 0.0),
      term(data.size(), 0.0);
  const std::size_t groups = members.size();
  const std::size_t blocks =
      std::min<std::size_t>(std::max<std::size_t>(groups, 1), kParallelBlocks);
  for_each_block(blocks, [&](std::size_t blk) {
    for (std::size_t hh = block_begin(groups, blocks, blk);
         hh < block_begin(groups, blocks, blk + 1); ++hh) {
      const std::vector<std::size_t> &obs = members[hh];
      if (obs.empty())
        continue;
      double num = 0.0, den = 0.0;
      for (std::size_t i : obs) {
        g[i] = model.inverse_scale(data.observation[i]);
        const double w = 1.0 / data.variance[i];
        num += w * g[i] * data.intensity[i];
        den += w * g[i] * g[i];
      }
      const double m = den > 0.0 ? num / den : 0.0;
      merged[hh] = m;
      for (std::size_t i : obs) {
        const double r = data.intensity[i] - g[i] * m;
        term[i] = r * r / data.variance[i];
      }
    }
  });
  double phi = 0.0;
  for (std::size_t i : use)
    phi += term[i];
  phi = restraint_of(model, options, phi);
  if (g_out)
    *g_out = std::move(g);
  if (merged_out)
    *merged_out = std::move(merged);
  return phi;
}

} // namespace

namespace {

//: The normal equations of the fit at the model's current parameters, over the
//: groups in `members`: J^T J and -J^T r, with J the variable-projection
//: Jacobian of the residuals sqrt(w)(I - g <I>), and the restraints. `g` and
//: `merged` are the model's scales and merged intensities over the same
//: observations. The fit steps with these; the covariance inverts them.
void normal_equations(const ScaleModel &model, const ScaleData &data,
                      const std::vector<std::vector<std::size_t>> &members,
                      const std::vector<double> &g,
                      const std::vector<double> &merged,
                      const ScaleFitOptions &options, std::vector<double> &N,
                      std::vector<double> &b) {
  const std::size_t n = model.size();
  N.assign(n * n, 0.0);
  b.assign(n, 0.0);
  // In blocks of groups, each with its own N and b, combined in block order:
  // the same sums on any number of threads. (The only step of scaling whose
  // order of summation differs from one thread's -- by rounding, the same on
  // every count.)
  const std::size_t groups = members.size();
  const std::size_t blocks =
      std::min<std::size_t>(std::max<std::size_t>(groups, 1), kParallelBlocks);
  std::vector<std::vector<double>> Nb(blocks), bb(blocks);
  for_each_block(blocks, [&](std::size_t blk) {
    // N and b here are this block's.
    std::vector<double> &N = Nb[blk];
    std::vector<double> &b = bb[blk];
    N.assign(n * n, 0.0);
    b.assign(n, 0.0);
    std::vector<std::pair<std::size_t, double>> grad;
    std::vector<double> G(n, 0.0);
    std::vector<std::vector<std::pair<std::size_t, double>>> grads;
    std::vector<std::size_t> touched;
    for (std::size_t hh = block_begin(groups, blocks, blk);
         hh < block_begin(groups, blocks, blk + 1); ++hh) {
      const std::vector<std::size_t> &obs = members[hh];
      if (obs.empty())
        continue;
      const std::size_t h = data.group[obs.front()];
      const double m = merged[h];
      grads.resize(obs.size());
      std::fill(G.begin(), G.end(), 0.0);
      touched.clear();
      double den = 0.0;
      for (std::size_t k = 0; k < obs.size(); ++k) {
        const std::size_t i = obs[k];
        const double w = 1.0 / data.variance[i];
        model.inverse_scale(data.observation[i], &grads[k]);
        den += w * g[i] * g[i];
        for (const auto &[a, da] : grads[k]) {
          if (G[a] == 0.0)
            touched.push_back(a);
          G[a] += w * g[i] * da;
        }
      }
      if (!(den > 0.0))
        continue;
      std::sort(touched.begin(), touched.end());
      touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
      for (std::size_t a : touched)
        G[a] /= den;
      // The group's part, by its rows d_i - g_i G: sum w m^2 (d_i - g_i G)
      // (d_i - g_i G)^T is m^2 (sum w d_i d_i^T - den G G^T), since sum w g_i
      // d_i is den G; and sum w m r_i (d_i - g_i G) is m (sum w r_i d_i - (sum
      // w g_i r_i) G). Each observation's own few parameters, then the dense G
      // G^T once a group, not once an observation: G touches every parameter
      // any of the group's observations does, at multiplicity 374 over ten
      // turns all 448, and N had cost 448^2 an observation -- three fits, 70 s
      // of mxi_scale's 95 on Graeme's 36000 images. The same sum, rearranged.
      double wgr = 0.0;
      for (std::size_t k = 0; k < obs.size(); ++k) {
        const std::size_t i = obs[k];
        const double w = 1.0 / data.variance[i];
        const double r = data.intensity[i] - g[i] * m;
        wgr += w * g[i] * r;
        for (const auto &[a, da] : grads[k]) {
          b[a] += w * m * r * da;
          double *row_a = &N[a * n];
          for (const auto &[c, dc] : grads[k])
            row_a[c] += w * m * m * da * dc;
        }
      }
      for (std::size_t a : touched) {
        b[a] -= m * wgr * G[a];
        const double ga = m * m * den * G[a];
        double *row_a = &N[a * n];
        for (std::size_t c : touched)
          row_a[c] -= ga * G[c];
      }
      for (std::size_t a : touched)
        G[a] = 0.0;
    }
  });
  for (std::size_t blk = 0; blk < blocks; ++blk) {
    for (std::size_t k = 0; k < n * n; ++k)
      N[k] += Nb[blk][k];
    for (std::size_t k = 0; k < n; ++k)
      b[k] += bb[blk][k];
  }

  for (std::size_t w = 0; w < model.sweeps(); ++w) {
    const ScaleModelShape &s = model.shape(w);
    for (std::size_t k = 0; k < s.decay_points; ++k) {
      const std::size_t a = model.first_decay(w) + k;
      N[a * n + a] += options.decay_restraint;
      b[a] -= options.decay_restraint * model.parameters[a];
    }
    for (std::size_t k = 0; k < harmonic_count(s.lmax); ++k) {
      const std::size_t a = model.first_absorption(w) + k;
      N[a * n + a] += options.absorption_restraint;
      b[a] -= options.absorption_restraint * model.parameters[a];
    }
  }
}

} // namespace

ScaleFitResult fit_scale_model(ScaleModel &model, const ScaleData &data,
                               const ScaleFitOptions &options,
                               const std::vector<std::size_t> &given) {
  std::vector<std::size_t> use = given;
  if (use.empty())
    for (std::size_t i = 0; i < data.size(); ++i)
      if (!data.outlier[i])
        use.push_back(i);
  ScaleFitResult result;
  result.observations = use.size();
  const std::size_t n = model.size();
  if (use.empty() || n == 0)
    return result;

  // The observations of each group being fitted.
  std::vector<std::vector<std::size_t>> members(data.unique.size());
  for (std::size_t i : use)
    members[data.group[i]].push_back(i);
  std::vector<double> g, merged;
  double phi = target(model, data, members, use, options, &g, &merged);
  result.target_start = phi;
  double mu = 1e-3;
  for (int it = 0; it < options.max_iterations; ++it) {
    result.iterations = it + 1;
    // Normal equations of the residuals sqrt(w)(I - g <I>), differentiated
    // with <I> as a function of the parameters (variable projection, in
    // Kaufman's form): d<I>/dp = -<I> G_h, G_h = sum w g dg/dp / sum w g^2.
    // Holding <I> instead gives the right gradient but overstates the
    // curvature in every direction <I> can partly follow, so the steps come
    // out short: three fits of 50 steps on a 300 image sweep ended still
    // creeping. The row for observation i is dg_i/dp - g_i G_h.
    std::vector<double> N, b;
    normal_equations(model, data, members, g, merged, options, N, b);
    bool improved = false;
    for (int tries = 0; tries < 12 && !improved; ++tries) {
      std::vector<double> A = N, step = b;
      for (std::size_t a = 0; a < n; ++a)
        A[a * n + a] += mu * std::fmax(N[a * n + a], 1e-12);
      if (!solve_spd(A.data(), step.data(), n)) {
        mu *= 10.0;
        continue;
      }
      ScaleModel trial = model;
      for (std::size_t a = 0; a < n; ++a)
        trial.parameters[a] += step[a];
      trial.normalise();
      std::vector<double> g2, m2;
      const double phi2 = target(trial, data, members, use, options, &g2, &m2);
      if (phi2 < phi) {
        const double gain = phi - phi2;
        model = trial;
        g = std::move(g2);
        merged = std::move(m2);
        phi = phi2;
        mu = std::fmax(mu / 10.0, 1e-9);
        improved = true;
        if (gain < 1e-8 * std::fmax(1.0, phi)) {
          result.converged = true;
          result.target_end = phi;
          return result;
        }
      } else {
        mu *= 10.0;
      }
    }
    if (!improved) {
      result.converged = true; // no step lowers the target: a minimum
      break;
    }
  }
  result.target_end = phi;
  return result;
}

namespace {

//: Weighted sums over a group's observations that are not outliers.
struct GroupSums {
  std::vector<double> num, den;
  std::vector<std::size_t> count;
};

GroupSums group_sums(const ScaleData &data, const std::vector<double> &g) {
  GroupSums s;
  s.num.assign(data.unique.size(), 0.0);
  s.den.assign(data.unique.size(), 0.0);
  s.count.assign(data.unique.size(), 0);
  for (std::size_t i = 0; i < data.size(); ++i) {
    if (data.outlier[i])
      continue;
    const double w = 1.0 / data.variance[i];
    s.num[data.group[i]] += w * g[i] * data.intensity[i];
    s.den[data.group[i]] += w * g[i] * g[i];
    ++s.count[data.group[i]];
  }
  return s;
}

//: The standard normal quantile, by Acklam's rational approximation refined
//: with a Newton step on erfc: better than 1e-12.
double normal_quantile(double p) {
  const double a[6] = {-3.969683028665376e+01, 2.209460984245205e+02,
                       -2.759285104469687e+02, 1.383577518672690e+02,
                       -3.066479806614716e+01, 2.506628277459239e+00};
  const double b[5] = {-5.447609879822406e+01, 1.615858368580409e+02,
                       -1.556989798598866e+02, 6.680131188771972e+01,
                       -1.328068155288572e+01};
  const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
                       -2.400758277161838e+00, -2.549732539343734e+00,
                       4.374664141464968e+00,  2.938163982698783e+00};
  const double d[4] = {7.784695709041462e-03, 3.224671290700398e-01,
                       2.445134137142996e+00, 3.754408661907416e+00};
  double x;
  if (p < 0.02425) {
    const double q = std::sqrt(-2.0 * std::log(p));
    x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  } else if (p > 1.0 - 0.02425) {
    const double q = std::sqrt(-2.0 * std::log(1.0 - p));
    x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  } else {
    const double q = p - 0.5, r = q * q;
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) *
        q /
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
  }
  const double e = 0.5 * std::erfc(-x / std::sqrt(2.0)) - p;
  const double u = e * std::sqrt(2.0 * std::acos(-1.0)) * std::exp(x * x / 2.0);
  return x - u / (1.0 + x * u / 2.0);
}

} // namespace

std::size_t reject_outliers(ScaleData &data, const std::vector<double> &g,
                            double zmax) {
  std::fill(data.outlier.begin(), data.outlier.end(), false);
  std::vector<std::vector<std::size_t>> members(data.unique.size());
  for (std::size_t i = 0; i < data.size(); ++i)
    members[data.group[i]].push_back(i);
  // Each group judged apart from every other: in parallel, in blocks of groups,
  // each block counting what it flags.
  const std::size_t groups = members.size();
  const std::size_t blocks =
      std::min<std::size_t>(std::max<std::size_t>(groups, 1), kParallelBlocks);
  std::vector<std::size_t> count(blocks, 0);
  for_each_block(blocks, [&](std::size_t blk) {
    for (std::size_t hh = block_begin(groups, blocks, blk);
         hh < block_begin(groups, blocks, blk + 1); ++hh) {
      const std::vector<std::size_t> &m = members[hh];
      while (true) {
        double sn = 0.0, sd = 0.0;
        std::size_t live = 0;
        for (std::size_t i : m) {
          if (data.outlier[i])
            continue;
          const double w = 1.0 / data.variance[i];
          sn += w * g[i] * data.intensity[i];
          sd += w * g[i] * g[i];
          ++live;
        }
        if (live < 3)
          break; // with two, which of them is wrong cannot be said
        double worst = 0.0;
        std::size_t at = m.front();
        for (std::size_t i : m) {
          if (data.outlier[i])
            continue;
          const double w = 1.0 / data.variance[i];
          const double num = sn - w * g[i] * data.intensity[i];
          const double den = sd - w * g[i] * g[i];
          const double z = std::abs(data.intensity[i] - g[i] * num / den) /
                           std::sqrt(data.variance[i] + g[i] * g[i] / den);
          if (z > worst) {
            worst = z;
            at = i;
          }
        }
        if (worst <= zmax)
          break;
        data.outlier[at] = true;
        ++count[blk];
      }
    }
  });
  std::size_t flagged = 0;
  for (std::size_t k : count)
    flagged += k;
  return flagged;
}

namespace {

//: The minimum of f on [lo, hi] by Brent's method (Brent 1973, "localmin"):
//: golden-section steps, and parabolic interpolation where it can be trusted,
//: to within tol.
template <typename F>
double brent_minimum(F f, double lo, double hi, double tol) {
  const double golden = 0.5 * (3.0 - std::sqrt(5.0));
  double a = lo, b = hi;
  double x = a + golden * (b - a), w = x, v = x;
  double fx = f(x), fw = fx, fv = fx;
  double d = 0.0, e = 0.0;
  for (int it = 0; it < 100; ++it) {
    const double m = 0.5 * (a + b);
    const double t1 = 1e-10 * std::abs(x) + tol / 3.0, t2 = 2.0 * t1;
    if (std::abs(x - m) <= t2 - 0.5 * (b - a))
      break;
    double step;
    bool parabolic = false;
    if (std::abs(e) > t1) {
      double r = (x - w) * (fx - fv), q = (x - v) * (fx - fw),
             pp = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0)
        pp = -pp;
      else
        q = -q;
      if (std::abs(pp) < std::abs(0.5 * q * e) && pp > q * (a - x) &&
          pp < q * (b - x)) {
        e = d;
        step = pp / q;
        const double u = x + step;
        if (u - a < t2 || b - u < t2)
          step = x < m ? t1 : -t1;
        parabolic = true;
      }
    }
    if (!parabolic) {
      e = (x < m ? b : a) - x;
      step = golden * e;
    }
    d = step;
    const double u = x + (std::abs(d) >= t1 ? d : (d > 0.0 ? t1 : -t1));
    const double fu = f(u);
    if (fu <= fx) {
      (u < x ? b : a) = x;
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    } else {
      (u < x ? a : b) = u;
      if (fu <= fw || w == x) {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }
  return x;
}

} // namespace

ErrorModel refine_error_model(const ScaleData &data,
                              const std::vector<double> &g) {
  const GroupSums s = group_sums(data, g);
  // Groups fit to judge an error model by, as dials.scale: <I> above 25 and
  // <I / sigma^2> above 0.85, with at least two observations.
  std::vector<double> ratio(data.unique.size(), 0.0);
  for (std::size_t i = 0; i < data.size(); ++i)
    if (!data.outlier[i])
      ratio[data.group[i]] +=
          data.intensity[i] / (data.variance_before[i] + data.scale_term_at(i));
  std::vector<std::size_t> use;
  for (std::size_t i = 0; i < data.size(); ++i) {
    const std::size_t h = data.group[i];
    if (data.outlier[i] || s.count[h] < 2 || !(s.den[h] > 0.0))
      continue;
    const double mean = s.num[h] / s.den[h];
    if (mean > 25.0 && ratio[h] / static_cast<double>(s.count[h]) > 0.85)
      use.push_back(i);
  }
  ErrorModel em;
  em.used = use.size();
  if (use.size() < 50)
    return em;
  // Normalised deviations of x_i = I_i / g_i from the merged <I_h>, with the
  // EXACT variance of the difference rather than eqn 12's prefactor. With
  // <I_h> = sum c_j x_j, weights held as dials.scale holds them, and v_j the
  // variance of x_j under the candidate error model,
  //     Var(x_i - <I_h>) = v_i (1 - 2 c_i) + sum_j c_j^2 v_j,
  // which is v (n - 1) / n when every weight is equal. dials.scale multiplies
  // by sqrt((n - 1) / n) instead (calc_deltahl), which leaves the deviations'
  // spread at (n - 1) / n and a low by that factor: planted a = 1.3 in groups
  // of eight came back 1.158.
  std::vector<double> c(data.size(), 0.0);
  for (std::size_t i : use) {
    const std::size_t h = data.group[i];
    c[i] = (g[i] * g[i] / data.variance[i]) / s.den[h];
  }
  // Everything in a deviation but a and b, once a fit rather than once an
  // evaluation. With A = (sigma^2 + the scale's term) / g^2 and B = I^2 / g^2,
  // a variance under the error model is v = a^2 (A + b^2 B), and the variance
  // of a deviation v (1 - 2c) + sum_h c^2 v = a^2 (P + b^2 Q), P and Q summed a
  // group at a time here -- so an evaluation is a square root an observation.
  // It was every group's spread summed again, for each of some sixty values of
  // b a round.
  std::vector<double> A_of(data.size(), 0.0), B_of(data.size(), 0.0);
  std::vector<double> SA(data.unique.size(), 0.0), SB(data.unique.size(), 0.0);
  for (std::size_t i : use) {
    const double g2 = g[i] * g[i];
    A_of[i] = (data.variance_before[i] + data.scale_term_at(i)) / g2;
    // The fractional term's intensity the model's, or else the reflection's
    // merged one -- not the observation's own, whose square overstates the
    // truth's by its own variance.
    const double Ih = s.num[data.group[i]] / s.den[data.group[i]];
    B_of[i] = data.model_intensity.empty()
                  ? Ih * Ih
                  : data.b_intensity(i) * data.b_intensity(i) / g2;
    SA[data.group[i]] += c[i] * c[i] * A_of[i];
    SB[data.group[i]] += c[i] * c[i] * B_of[i];
  }
  std::vector<double> P(use.size()), Q(use.size()), R(use.size());
  for (std::size_t k = 0; k < use.size(); ++k) {
    const std::size_t i = use[k], h = data.group[i];
    P[k] = A_of[i] * (1.0 - 2.0 * c[i]) + SA[h];
    Q[k] = B_of[i] * (1.0 - 2.0 * c[i]) + SB[h];
    R[k] = data.intensity[i] / g[i] - s.num[h] / s.den[h];
  }
  // Intensity bins for b: logarithmically spaced over the used intensities.
  double lo = HUGE_VAL, hi = 0.0;
  for (std::size_t i : use) {
    const double v =
        std::fmax(s.num[data.group[i]] / s.den[data.group[i]], 1e-3);
    lo = std::fmin(lo, v);
    hi = std::fmax(hi, v);
  }
  const int nbins = 20;
  std::vector<int> bin(use.size());
  for (std::size_t k = 0; k < use.size(); ++k) {
    // Binned by the reflection's merged intensity, as dials.scale: by the
    // observation's own, one that scattered high joined a higher bin, and
    // every bin's deviations read wide -- a = 1.18 for 1.0, b = 0.15 planted.
    const std::size_t h = data.group[use[k]];
    const double v = std::fmax(s.num[h] / s.den[h], 1e-3);
    bin[k] = hi > lo ? std::min(nbins - 1,
                                static_cast<int>(nbins * std::log(v / lo) /
                                                 std::log(hi / lo)))
                     : 0;
  }
  // a and b together, by regression, as dials.scale's default: in each
  // intensity bin the deviations' mean square should be a^2 (P + b^2 Q),
  // linear in a^2 and a^2 b^2 -- fitted across the bins by weighted least
  // squares, each bin's weight its count over its mean square squared. a had
  // been the slope of the normal probability plot over every observation,
  // b then fitted about it: on Graeme's 6S0R set, its strong reflections out by
  // a tenth of their intensity, that gave a = 4.81 and b = 0.023 -- every sigma
  // five times over, the weak ones' three times too large (the anomalous
  // slope 0.33, its control 0.34) -- where the data wanted b near 0.1.
  std::vector<double> sy(nbins, 0.0), sp(nbins, 0.0), sq(nbins, 0.0);
  std::vector<std::size_t> count(nbins, 0);
  for (std::size_t k = 0; k < use.size(); ++k) {
    sy[bin[k]] += R[k] * R[k];
    sp[bin[k]] += P[k];
    sq[bin[k]] += Q[k];
    ++count[bin[k]];
  }
  double Spp = 0.0, Spq = 0.0, Sqq = 0.0, Spy = 0.0, Sqy = 0.0;
  for (int j = 0; j < nbins; ++j) {
    if (count[j] < 10)
      continue;
    const double c = static_cast<double>(count[j]);
    const double y = sy[j] / c, pp = sp[j] / c, qq = sq[j] / c;
    if (!(y > 0.0))
      continue;
    const double w = c / (y * y);
    Spp += w * pp * pp;
    Spq += w * pp * qq;
    Sqq += w * qq * qq;
    Spy += w * pp * y;
    Sqy += w * qq * y;
  }
  double alpha = 0.0, beta = 0.0;
  const double det = Spp * Sqq - Spq * Spq;
  // The fit, bin by bin, under --verbose: what drove a and b.
  const auto show = [&](double al, double be) {
    if (!verbose())
      return;
    std::vector<double> mean_ih(nbins, 0.0);
    for (std::size_t k = 0; k < use.size(); ++k) {
      const std::size_t h = data.group[use[k]];
      mean_ih[bin[k]] += s.num[h] / s.den[h];
    }
    detail("  the error model's fit, a bin at a time: y the deviations' mean "
           "square, p the variance as given, q the intensity squared\n");
    detail("  %4s %8s %10s %12s %12s %12s %12s %8s\n", "bin", "n", "<I>", "y",
           "p", "q", "fit", "y/fit");
    for (int j = 0; j < nbins; ++j) {
      if (count[j] == 0)
        continue;
      const double c = static_cast<double>(count[j]);
      const double fit = al * sp[j] / c + be * sq[j] / c;
      detail("  %4d %8zu %10.1f %12.4g %12.4g %12.4g %12.4g %8.3f%s\n", j,
             count[j], mean_ih[j] / c, sy[j] / c, sp[j] / c, sq[j] / c, fit,
             fit > 0.0 ? (sy[j] / c) / fit : 0.0,
             count[j] < 10 ? "  (too few: not fitted)" : "");
    }
    detail("  a^2 = %.5g, a^2 b^2 = %.5g\n", al, be);
  };
  if (det > 0.0) {
    alpha = (Spy * Sqq - Sqy * Spq) / det;
    beta = (Sqy * Spp - Spy * Spq) / det;
  }
  if (!(beta > 0.0) || !(alpha > 0.0)) { // b at its bound, a alone
    beta = 0.0;
    alpha = Spp > 0.0 ? Spy / Spp : 1.0;
  }
  show(alpha, beta);
  if (alpha > 0.0) {
    em.a = std::sqrt(alpha);
    em.b = std::min(std::sqrt(beta / alpha), 0.5);
  }
  em.a_bg = em.a;
  // Split, where integration wrote each variance with no signal: one a for a
  // background's variance, one for a signal's. One a for both had taken the
  // strong reflections' excess -- 6S0R's, a tenth of their intensity -- into
  // the weak ones' sigmas, whose variances were near right: final calibration
  // 0.05 to 0.08 (outstanding item 64). Every group of two or more, the weak
  // included, in twenty bins of merged intensity, equal counts.
  if (data.has_background && data.variance_background.size() == data.size()) {
    std::vector<std::size_t> all;
    for (std::size_t i = 0; i < data.size(); ++i) {
      const std::size_t h = data.group[i];
      if (!data.outlier[i] && s.count[h] >= 2 && s.den[h] > 0.0 && g[i] > 0.0)
        all.push_back(i);
    }
    if (all.size() >= 400) {
      std::vector<double> cc(data.size(), 0.0), Vb(data.size(), 0.0),
          Vs(data.size(), 0.0), Iq(data.size(), 0.0);
      std::vector<double> SVb(data.unique.size(), 0.0),
          SVs(data.unique.size(), 0.0), SIq(data.unique.size(), 0.0);
      for (std::size_t i : all) {
        const std::size_t h = data.group[i];
        const double g2 = g[i] * g[i];
        cc[i] = (g2 / data.variance[i]) / s.den[h];
        const double bg = std::fmax(data.variance_background[i], 0.0);
        Vb[i] = bg / g2;
        Vs[i] = std::fmax(data.variance_before[i] + data.scale_term_at(i) - bg,
                          0.0) /
                g2;
        const double Ih = s.num[h] / s.den[h];
        Iq[i] = data.model_intensity.empty()
                    ? Ih * Ih
                    : data.b_intensity(i) * data.b_intensity(i) / g2;
        SVb[h] += cc[i] * cc[i] * Vb[i];
        SVs[h] += cc[i] * cc[i] * Vs[i];
        SIq[h] += cc[i] * cc[i] * Iq[i];
      }
      std::sort(all.begin(), all.end(), [&](std::size_t x, std::size_t y) {
        return s.num[data.group[x]] / s.den[data.group[x]] <
               s.num[data.group[y]] / s.den[data.group[y]];
      });
      const int nb = 20;
      double M[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, v[3] = {0, 0, 0};
      struct Row {
        std::size_t n;
        double ih, y, pb, ps, q;
      };
      std::vector<Row> rows;
      for (int j = 0; j < nb; ++j) {
        const std::size_t from = all.size() * static_cast<std::size_t>(j) / nb,
                          to =
                              all.size() * static_cast<std::size_t>(j + 1) / nb;
        Row r{to - from, 0, 0, 0, 0, 0};
        for (std::size_t k = from; k < to; ++k) {
          const std::size_t i = all[k], h = data.group[i];
          const double dev = data.intensity[i] / g[i] - s.num[h] / s.den[h];
          const double f = 1.0 - 2.0 * cc[i];
          r.ih += s.num[h] / s.den[h];
          r.y += dev * dev;
          r.pb += Vb[i] * f + SVb[h];
          r.ps += Vs[i] * f + SVs[h];
          r.q += Iq[i] * f + SIq[h];
        }
        const double n = static_cast<double>(r.n);
        r.ih /= n;
        r.y /= n;
        r.pb /= n;
        r.ps /= n;
        r.q /= n;
        rows.push_back(r);
        if (!(r.y > 0.0))
          continue;
        const double w = n / (r.y * r.y), x[3] = {r.pb, r.ps, r.q};
        for (int a = 0; a < 3; ++a) {
          v[a] += w * x[a] * r.y;
          for (int b = 0; b < 3; ++b)
            M[a][b] += w * x[a] * x[b];
        }
      }
      // Solved by Cholesky; b at its bound refitted without it.
      const auto solve = [](double A[3][3], double rhs[3], int n, double *out) {
        std::vector<double> m(static_cast<std::size_t>(n * n)), r(rhs, rhs + n);
        for (int a = 0; a < n; ++a)
          for (int b = 0; b < n; ++b)
            m[static_cast<std::size_t>(a * n + b)] = A[a][b];
        if (!solve_spd(m.data(), r.data(), static_cast<std::size_t>(n)))
          return false;
        for (int a = 0; a < n; ++a)
          out[a] = r[static_cast<std::size_t>(a)];
        return true;
      };
      double x3[3] = {0, 0, 0};
      bool ok = solve(M, v, 3, x3);
      if (ok && !(x3[2] > 0.0)) {
        x3[2] = 0.0;
        ok = solve(M, v, 2, x3);
      }
      if (verbose()) {
        detail("  the error model split, a bin at a time: y the deviations' "
               "mean square, background's and signal's variances, intensity "
               "squared\n");
        detail("  %4s %8s %10s %12s %12s %12s %12s %12s %8s\n", "bin", "n",
               "<I>", "y", "bg", "signal", "q", "fit", "y/fit");
        for (std::size_t j = 0; j < rows.size(); ++j) {
          const Row &r = rows[j];
          const double fit = x3[0] * r.pb + x3[1] * r.ps + x3[2] * r.q;
          detail(
              "  %4zu %8zu %10.1f %12.4g %12.4g %12.4g %12.4g %12.4g %8.3f\n",
              j, r.n, r.ih, r.y, r.pb, r.ps, r.q, fit,
              fit > 0.0 ? r.y / fit : 0.0);
        }
        detail("  a_bg^2 = %.5g, a^2 = %.5g, a^2 b^2 = %.5g%s\n", x3[0], x3[1],
               x3[2], ok && x3[0] > 0.0 && x3[1] > 0.0 ? "" : ": not kept");
      }
      if (ok && x3[0] > 0.0 && x3[1] > 0.0) {
        em.split = true;
        em.a_bg = std::sqrt(x3[0]);
        em.a = std::sqrt(x3[1]);
        em.b = std::min(std::sqrt(std::fmax(x3[2], 0.0) / x3[1]), 0.5);
        em.used = all.size();
      }
    }
  }
  return em;
}

void apply_error_model(ScaleData &data, const ErrorModel &model) {
  const bool split = model.split && data.has_background &&
                     data.variance_background.size() == data.size();
  for (std::size_t i = 0; i < data.size(); ++i) {
    const double I = data.b_intensity(i);
    if (split) {
      const double bg = std::fmax(data.variance_background[i], 0.0);
      const double rest =
          std::fmax(data.variance_before[i] + data.scale_term_at(i) - bg, 0.0);
      data.variance[i] = model.a_bg * model.a_bg * bg +
                         model.a * model.a * (rest + model.b * model.b * I * I);
    } else {
      data.variance[i] = model.a * model.a *
                         (data.variance_before[i] + data.scale_term_at(i) +
                          model.b * model.b * I * I);
    }
  }
}

std::vector<CalibrationRow> calibration(const ScaleData &data,
                                        const std::vector<double> &g,
                                        bool by_resolution, int bins) {
  const std::size_t groups = data.unique.size();
  std::vector<std::vector<std::size_t>> members(groups);
  for (std::size_t i = 0; i < data.size(); ++i)
    if (!data.outlier[i] && g[i] > 0.0)
      members[data.group[i]].push_back(i);
  struct Dev {
    double key, raw, final;
  };
  std::vector<Dev> devs;
  for (std::size_t h = 0; h < groups; ++h) {
    const std::vector<std::size_t> &m = members[h];
    if (m.size() < 2)
      continue;
    double W[2] = {0, 0}, S[2] = {0, 0};
    for (std::size_t i : m) {
      const double x = data.intensity[i] / g[i];
      const double v[2] = {data.variance_before[i] + data.scale_term_at(i),
                           data.variance[i]};
      for (int k = 0; k < 2; ++k)
        if (v[k] > 0.0) {
          const double w = g[i] * g[i] / v[k];
          W[k] += w;
          S[k] += w * x;
        }
    }
    const double mean_final = W[1] > 0.0 ? S[1] / W[1] : 0.0;
    for (std::size_t i : m) {
      const double x = data.intensity[i] / g[i];
      const double v[2] = {data.variance_before[i] + data.scale_term_at(i),
                           data.variance[i]};
      double d[2] = {0, 0};
      bool ok = true;
      for (int k = 0; k < 2; ++k) {
        const double w = v[k] > 0.0 ? g[i] * g[i] / v[k] : 0.0;
        const double rest = W[k] - w;
        if (!(v[k] > 0.0) || !(rest > 0.0)) {
          ok = false;
          break;
        }
        const double others = (S[k] - w * x) / rest;
        d[k] = (x - others) / std::sqrt(v[k] / (g[i] * g[i]) + 1.0 / rest);
      }
      if (ok)
        devs.push_back({by_resolution ? data.d[i] : mean_final, d[0], d[1]});
    }
  }
  std::vector<CalibrationRow> out;
  if (devs.size() < static_cast<std::size_t>(bins) * 10)
    return out;
  // Equal counts a bin: resolution from low to high, intensity from weak up.
  std::sort(devs.begin(), devs.end(), [&](const Dev &a, const Dev &b) {
    return by_resolution ? a.key > b.key : a.key < b.key;
  });
  const auto robust = [](std::vector<double> x) {
    std::nth_element(x.begin(), x.begin() + x.size() / 2, x.end());
    const double median = x[x.size() / 2];
    for (double &v : x)
      v = std::abs(v - median);
    std::nth_element(x.begin(), x.begin() + x.size() / 2, x.end());
    const double mad = 1.4826 * x[x.size() / 2];
    return mad * mad;
  };
  for (int b = 0; b < bins; ++b) {
    const std::size_t from = devs.size() * static_cast<std::size_t>(b) / bins,
                      to = devs.size() * static_cast<std::size_t>(b + 1) / bins;
    CalibrationRow r;
    r.lo = devs[from].key;
    r.hi = devs[to - 1].key;
    r.n = to - from;
    std::vector<double> raw, fin;
    for (std::size_t k = from; k < to; ++k) {
      r.raw_ms += devs[k].raw * devs[k].raw;
      r.final_ms += devs[k].final * devs[k].final;
      raw.push_back(devs[k].raw);
      fin.push_back(devs[k].final);
    }
    r.raw_ms /= static_cast<double>(r.n);
    r.final_ms /= static_cast<double>(r.n);
    r.raw_robust = robust(raw);
    r.final_robust = robust(fin);
    out.push_back(r);
  }
  return out;
}

bool model_variances(ScaleData &data, const std::vector<double> &g) {
  if (!data.has_background || data.variance_background.size() != data.size())
    return false;
  const GroupSums s = group_sums(data, g);
  data.model_intensity.assign(data.size(), 0.0);
  for (std::size_t i = 0; i < data.size(); ++i) {
    const std::size_t h = data.group[i];
    const double mean = s.den[h] > 0.0 ? s.num[h] / s.den[h] : 0.0;
    const double m = g[i] * mean;
    data.model_intensity[i] = m;
    data.variance_before[i] =
        data.variance_background[i] + data.factor[i] * std::fmax(m, 0.0);
    data.variance[i] = data.variance_before[i] + data.scale_term_at(i);
  }
  return true;
}

void combine_intensities(ScaleData &data, double i_mid) {
  for_each_index(data.size(), [&](std::size_t i) {
    double w = 1.0;
    if (!data.has_sum[i])
      w = 1.0;
    else if (std::isinf(i_mid))
      w = 0.0;
    else if (i_mid > 0.0)
      w = 1.0 / (1.0 + std::pow(std::fmax(data.sum[i], 0.0) / i_mid, 3.0));
    data.intensity[i] = w * data.prf[i] + (1.0 - w) * data.sum[i];
    const double sp = std::sqrt(data.prf_variance[i]),
                 ss = std::sqrt(data.sum_variance[i]);
    const double sigma = w * sp + (1.0 - w) * (data.has_sum[i] ? ss : sp);
    data.variance[i] = sigma * sigma;
    data.variance_before[i] = sigma * sigma;
    if (data.has_background) {
      // The variances with no signal, mixed as the variances are.
      const double bp = std::sqrt(std::fmax(data.prf_variance_bg[i], 0.0)),
                   bs = std::sqrt(std::fmax(data.sum_variance_bg[i], 0.0));
      const double b = w * bp + (1.0 - w) * (data.has_sum[i] ? bs : bp);
      data.variance_background[i] = b * b;
    }
  });
}

double rmeas(const ScaleData &data, const std::vector<double> &g) {
  const GroupSums s = group_sums(data, g);
  double top = 0.0, bottom = 0.0;
  for (std::size_t i = 0; i < data.size(); ++i) {
    const std::size_t h = data.group[i];
    if (data.outlier[i] || s.count[h] < 2 || !(s.den[h] > 0.0))
      continue;
    const double n = static_cast<double>(s.count[h]);
    const double mean = s.num[h] / s.den[h];
    top += std::sqrt(n / (n - 1.0)) * std::abs(data.intensity[i] / g[i] - mean);
    bottom += data.intensity[i] / g[i];
  }
  return bottom > 0.0 ? top / bottom : 0.0;
}

double choose_intensity_combination(ScaleData &data,
                                    const std::vector<double> &g) {
  double biggest = 1.0;
  for (std::size_t i = 0; i < data.size(); ++i)
    if (data.has_sum[i])
      biggest = std::fmax(biggest, data.sum[i]);
  std::vector<double> candidates = {0.0, HUGE_VAL};
  for (double m = 10.0; m < biggest; m *= 10.0)
    candidates.push_back(m);
  double best = 0.0, best_r = HUGE_VAL;
  for (double m : candidates) {
    combine_intensities(data, m);
    const double r = rmeas(data, g);
    if (r < best_r) {
      best_r = r;
      best = m;
    }
  }
  combine_intensities(data, best);
  return best;
}

namespace {

//: Plain and weighted means of scaled observations x = I/g.
struct Merged {
  std::size_t n = 0;
  double sum = 0.0;           //: sum of x
  double sw = 0.0, swx = 0.0; //: sum of w = g^2 / var, and of w x
  double mean() const { return n ? sum / static_cast<double>(n) : 0.0; }
  double weighted() const { return sw > 0.0 ? swx / sw : 0.0; }
  double sigma() const { return sw > 0.0 ? 1.0 / std::sqrt(sw) : 0.0; }
};

double pearson_of(const std::vector<std::pair<double, double>> &pairs) {
  if (pairs.size() < 3)
    return 0.0;
  const double n = static_cast<double>(pairs.size());
  double ma = 0.0, mb = 0.0;
  for (const auto &[a, b] : pairs) {
    ma += a / n;
    mb += b / n;
  }
  double sab = 0.0, saa = 0.0, sbb = 0.0;
  for (const auto &[a, b] : pairs) {
    sab += (a - ma) * (b - mb);
    saa += (a - ma) * (a - ma);
    sbb += (b - mb) * (b - mb);
  }
  return saa > 0.0 && sbb > 0.0 ? sab / std::sqrt(saa * sbb) : 0.0;
}

//: R factor sums over one set of observations of a group, plain mean.
void r_sums(const ScaleData &data, const std::vector<double> &g,
            const std::vector<std::size_t> &o, double *merge, double *meas,
            double *pim, double *denominator) {
  if (o.size() < 2)
    return;
  const double n = static_cast<double>(o.size());
  double mean = 0.0;
  for (std::size_t i : o)
    mean += data.intensity[i] / g[i] / n;
  double dev = 0.0, sum = 0.0;
  for (std::size_t i : o) {
    dev += std::abs(data.intensity[i] / g[i] - mean);
    sum += data.intensity[i] / g[i];
  }
  *merge += dev;
  *meas += std::sqrt(n / (n - 1.0)) * dev;
  *pim += std::sqrt(1.0 / (n - 1.0)) * dev;
  *denominator += sum;
}

//: The mean of a random half of o, and of the rest.
std::pair<double, double> halves(const ScaleData &data,
                                 const std::vector<double> &g,
                                 std::vector<std::size_t> o,
                                 std::mt19937 &rng) {
  std::shuffle(o.begin(), o.end(), rng);
  const std::size_t split = o.size() / 2;
  double a = 0.0, b = 0.0;
  for (std::size_t k = 0; k < o.size(); ++k)
    (k < split ? a : b) += data.intensity[o[k]] / g[o[k]];
  return {a / static_cast<double>(split),
          b / static_cast<double>(o.size() - split)};
}

} // namespace

MergingShell merge_groups(const ScaleData &data, const std::vector<double> &g,
                          const std::vector<std::size_t> &groups) {
  MergingShell m;
  std::vector<std::vector<std::size_t>> all(data.stats_unique().size()),
      plus(data.stats_unique().size()), minus(data.stats_unique().size());
  std::vector<bool> wanted(data.stats_unique().size(), false);
  for (std::size_t h : groups)
    wanted[h] = true;
  for (std::size_t i = 0; i < data.size(); ++i) {
    const std::size_t h = data.stats_group()[i];
    if (data.outlier[i] || !wanted[h])
      continue;
    all[h].push_back(i);
    const bool centric =
        h < data.stats_centric().size() && data.stats_centric()[h];
    const bool p = i < data.plus.size() ? data.plus[i] : true;
    (centric || p ? plus : minus)[h].push_back(i);
  }
  std::mt19937 rng(20);
  double denominator = 0.0, denominator_anom = 0.0;
  std::size_t anomalous_groups = 0;
  std::vector<std::pair<double, double>> cc_pairs, anom_pairs;
  std::vector<double> deltas;        // dI / sigma(dI)
  std::vector<double> deltas_random; // the same, random halves of each pair
  std::mt19937 dealer(21);
  double sum_abs_di = 0.0, sum_sig_di = 0.0, sum_df2 = 0.0, sum_f2 = 0.0;
  std::size_t n_df = 0;
  const auto merged_of = [&](const std::vector<std::size_t> &o) {
    Merged r;
    for (std::size_t i : o) {
      const double x = data.intensity[i] / g[i];
      const double w = g[i] * g[i] / data.variance[i];
      ++r.n;
      r.sum += x;
      r.sw += w;
      r.swx += w * x;
    }
    return r;
  };
  for (std::size_t h : groups) {
    if (all[h].empty())
      continue;
    const Merged mg = merged_of(all[h]);
    m.observations += mg.n;
    ++m.unique;
    m.mean_i += mg.weighted();
    m.i_over_sigma += mg.weighted() / mg.sigma();
    r_sums(data, g, all[h], &m.rmerge, &m.rmeas, &m.rpim, &denominator);
    if (all[h].size() >= 2)
      cc_pairs.push_back(halves(data, g, all[h], rng));
    // With Friedel mates apart: a centric reflection is one group, an acentric
    // one two.
    for (const std::vector<std::size_t> *side : {&plus[h], &minus[h]}) {
      if (side->empty())
        continue;
      ++anomalous_groups;
      r_sums(data, g, *side, &m.rmerge_anom, &m.rmeas_anom, &m.rpim_anom,
             &denominator_anom);
    }
    const bool centric =
        h < data.stats_centric().size() && data.stats_centric()[h];
    if (centric || plus[h].empty() || minus[h].empty())
      continue;
    ++m.anomalous_pairs;
    const Merged p = merged_of(plus[h]), q = merged_of(minus[h]);
    const double di = p.weighted() - q.weighted();
    const double sdi = std::sqrt(p.sigma() * p.sigma() + q.sigma() * q.sigma());
    sum_abs_di += std::abs(di);
    sum_sig_di += sdi;
    if (sdi > 0.0)
      deltas.push_back(di / sdi);
    {
      // The control: this pair's observations pooled and dealt into two
      // groups as large as I+'s and I-'s.
      std::vector<std::size_t> pooled = plus[h];
      pooled.insert(pooled.end(), minus[h].begin(), minus[h].end());
      std::shuffle(pooled.begin(), pooled.end(), dealer);
      const std::vector<std::size_t> one(pooled.begin(),
                                         pooled.begin() + plus[h].size()),
          two(pooled.begin() + plus[h].size(), pooled.end());
      const Merged a = merged_of(one), b = merged_of(two);
      const double s2 =
          std::sqrt(a.sigma() * a.sigma() + b.sigma() * b.sigma());
      if (s2 > 0.0)
        deltas_random.push_back((a.weighted() - b.weighted()) / s2);
    }
    if (p.weighted() > 0.0 && q.weighted() > 0.0) {
      const double fp = std::sqrt(p.weighted()), fm = std::sqrt(q.weighted());
      sum_df2 += (fp - fm) * (fp - fm);
      sum_f2 += fp * fp + fm * fm;
      ++n_df;
    }
    if (plus[h].size() >= 2 && minus[h].size() >= 2) {
      const auto a = halves(data, g, plus[h], rng),
                 b = halves(data, g, minus[h], rng);
      anom_pairs.emplace_back(a.first - b.first, a.second - b.second);
    }
  }
  if (m.unique) {
    m.multiplicity = static_cast<double>(m.observations) / m.unique;
    m.mean_i /= m.unique;
    m.i_over_sigma /= m.unique;
  }
  if (denominator > 0.0) {
    m.rmerge /= denominator;
    m.rmeas /= denominator;
    m.rpim /= denominator;
  }
  if (denominator_anom > 0.0) {
    m.rmerge_anom /= denominator_anom;
    m.rmeas_anom /= denominator_anom;
    m.rpim_anom /= denominator_anom;
  }
  m.cc_half = pearson_of(cc_pairs);
  m.cc_anom = pearson_of(anom_pairs);
  m.anom_multiplicity =
      anomalous_groups ? static_cast<double>(m.observations) / anomalous_groups
                       : 0.0;
  m.di_over_sig_di = sum_sig_di > 0.0 ? sum_abs_di / sum_sig_di : 0.0;
  m.df_over_f = n_df && sum_f2 > 0.0 ? std::sqrt(2.0 * sum_df2 / sum_f2) : 0.0;
  // The slope of the normal probability plot over |x| < 0.9, as dials.scale:
  // the Friedel pairs', and the control's.
  const auto slope_of = [](std::vector<double> d) {
    if (d.size() < 10)
      return 0.0;
    std::sort(d.begin(), d.end());
    const double n = static_cast<double>(d.size());
    double sxx = 0.0, sxy = 0.0, sx = 0.0, sy = 0.0, k = 0.0;
    for (std::size_t i = 0; i < d.size(); ++i) {
      const double x = normal_quantile((static_cast<double>(i) + 0.5) / n);
      if (std::abs(x) < 0.9) {
        sx += x;
        sy += d[i];
        sxx += x * x;
        sxy += x * d[i];
        k += 1.0;
      }
    }
    const double den = sxx - sx * sx / k;
    return k > 2.0 && den > 0.0 ? (sxy - sx * sy / k) / den : 0.0;
  };
  m.anom_slope = slope_of(deltas);
  m.anom_slope_random = slope_of(deltas_random);
  return m;
}

std::vector<MergingShell> merging_statistics(const ScaleData &data,
                                             const std::vector<double> &g,
                                             const SpaceGroup &group,
                                             const Crystal &crystal, int shells,
                                             MergingShell *overall) {
  shells = std::max(shells, 1);
  // One d for everything: the crystal's, at each group's unique index. The d
  // column is the scan-varying crystal's at each observation, and mixing the
  // two put a shell at 100.2 per cent.
  std::vector<bool> seen(data.stats_unique().size(), false);
  for (std::size_t i = 0; i < data.size(); ++i)
    if (!data.outlier[i])
      seen[data.stats_group()[i]] = true;
  std::vector<double> d_of(data.stats_unique().size(), 0.0);
  double lo = HUGE_VAL, hi = 0.0; // 1/d^3
  for (std::size_t h = 0; h < data.stats_unique().size(); ++h) {
    const Miller &u = data.stats_unique()[h];
    d_of[h] = crystal.d_spacing(u[0], u[1], u[2]);
    if (!seen[h])
      continue;
    const double v = 1.0 / (d_of[h] * d_of[h] * d_of[h]);
    lo = std::fmin(lo, v);
    hi = std::fmax(hi, v);
  }
  std::vector<MergingShell> out(static_cast<std::size_t>(shells));
  if (!(hi > 0.0))
    return out;
  const auto shell_of = [&](double d) {
    const double v = 1.0 / (d * d * d);
    const auto s =
        hi > lo ? static_cast<int>((v - lo) / (hi - lo) * shells) : 0;
    return static_cast<std::size_t>(std::min(std::max(s, 0), shells - 1));
  };
  std::vector<std::vector<std::size_t>> in_shell(out.size());
  std::vector<std::size_t> everything;
  for (std::size_t h = 0; h < data.stats_unique().size(); ++h) {
    if (!seen[h])
      continue;
    in_shell[shell_of(d_of[h])].push_back(h);
    everything.push_back(h);
  }
  for (std::size_t s = 0; s < out.size(); ++s) {
    out[s] = merge_groups(data, g, in_shell[s]);
    out[s].d_max =
        std::cbrt(1.0 / (lo + (hi - lo) * static_cast<double>(s) / shells));
    out[s].d_min =
        std::cbrt(1.0 / (lo + (hi - lo) * static_cast<double>(s + 1) / shells));
  }
  MergingShell all = merge_groups(data, g, everything);
  all.d_max = out.front().d_max;
  all.d_min = out.back().d_min;
  // Possible reflections, and the acentric ones among them.
  const UnitCell cell = crystal.cell();
  const int hmax = static_cast<int>(std::ceil(cell.a / all.d_min)),
            kmax = static_cast<int>(std::ceil(cell.b / all.d_min)),
            lmax = static_cast<int>(std::ceil(cell.c / all.d_min));
  for (int hh = -hmax; hh <= hmax; ++hh)
    for (int kk = -kmax; kk <= kmax; ++kk)
      for (int ll = -lmax; ll <= lmax; ++ll) {
        if (hh == 0 && kk == 0 && ll == 0)
          continue;
        const double d = crystal.d_spacing(hh, kk, ll);
        if (d < all.d_min || d > all.d_max)
          continue;
        const Miller m{hh, kk, ll};
        if (group.absent(m) || group.unique(m) != m)
          continue;
        const std::size_t s = shell_of(d);
        const bool acentric = !group.centric(m);
        ++out[s].possible;
        ++all.possible;
        if (acentric) {
          ++out[s].possible_acentric;
          ++all.possible_acentric;
        }
      }
  for (MergingShell *m : [&] {
         std::vector<MergingShell *> v;
         for (MergingShell &s : out)
           v.push_back(&s);
         v.push_back(&all);
         return v;
       }()) {
    m->completeness =
        m->possible ? static_cast<double>(m->unique) / m->possible : 0.0;
    m->anom_completeness =
        m->possible_acentric
            ? std::fmin(1.0, static_cast<double>(m->anomalous_pairs) /
                                 m->possible_acentric)
            : 0.0;
  }
  if (overall)
    *overall = all;
  return out;
}

ScaleRun scale_sweep(const ExperimentList &experiments,
                     const Table &reflections, const SpaceGroup &group,
                     const ScaleRunOptions &options) {
  ScaleRun run;
  double mark = Timing::now();
  const auto step = [&](const char *name) {
    const double t = Timing::now();
    run.timing.emplace_back(name, t - mark);
    mark = t;
  };
  // A shape a sweep, from its own width of rotation; the harmonics gathered
  // to the largest any sweep has, each sweep reading as many as its own.
  std::vector<ScaleModelShape> shapes;
  ScaleModelShape shape;
  for (const Experiment &e : experiments) {
    const double degrees = std::abs(e.scan.osc_width) * e.scan.num_images();
    ScaleModelShape own = default_shape(degrees);
    if (options.lmax >= 0)
      own.lmax = options.lmax; // asked for, whatever the sweep's width
    else if (own.lmax > 0)
      own.lmax = options.level_lmax; // the absorption level's, wide sweeps only
    if (!options.absorption)
      own.lmax = 0;
    shapes.push_back(own);
    if (shapes.size() == 1 || own.lmax > shape.lmax)
      shape.lmax = own.lmax;
  }
  shape.scale_points = shapes[0].scale_points;
  shape.decay_points = shapes[0].decay_points;
  ScaleDataOptions data_options;
  data_options.d_min = options.d_min;
  data_options.d_max = options.d_max;
  data_options.anomalous = options.anomalous;
  run.data =
      build_scale_data(experiments, reflections, group, shape, data_options);
  step("gathering the observations");
  ScaleData &data = run.data;
  run.model = ScaleModel(shapes);
  if (data.size() == 0)
    return run;
  // A sweep with nothing to scale has parameters nothing determines.
  {
    std::vector<std::size_t> per(experiments.size(), 0);
    for (const ScaleObservation &o : data.observation)
      ++per[o.sweep];
    for (std::size_t w = 0; w < per.size(); ++w)
      if (per[w] == 0)
        throw std::runtime_error("sweep " + std::to_string(w) +
                                 " has no observations to scale");
  }

  const auto fit = [&]() {
    const std::vector<std::size_t> use = select_for_fitting(data);
    run.fitted_on = use.size();
    run.fits.push_back(fit_scale_model(run.model, data, options.fit, use));
    run.g = inverse_scales(run.model, data);
  };
  run.g.assign(data.size(), 1.0);
  reject_outliers(data, run.g); // on the unscaled intensities
  step("outliers, unscaled");
  fit();
  step("fit 1");
  reject_outliers(data, run.g);
  if (options.combine) {
    run.i_mid = choose_intensity_combination(data, run.g);
    reject_outliers(data, run.g);
  }
  step("outliers and profile against summation");
  // From here each variance from the model, not the observation's own noise.
  if (options.model_variances)
    model_variances(data, run.g);
  run.error_model = refine_error_model(data, run.g);
  apply_error_model(data, run.error_model);
  step("error model 1");
  fit();
  step("fit 2");
  reject_outliers(data, run.g);
  if (options.model_variances)
    model_variances(data, run.g);
  run.error_model = refine_error_model(data, run.g);
  apply_error_model(data, run.error_model);
  step("outliers and error model 2");
  fit();
  step("fit 3");
  // The scale's uncertainty into each observation's variance BEFORE the final
  // error model, so that a and b correct what remains after it rather than
  // absorbing it. (dials.scale applies its error model first and then
  // multiplies each variance by 1 + sigma_g / g, a factor linear in the
  // fractional error and the same at every intensity; propagated, the term is
  // I^2 var(g) / g^2.)
  run.covariance = parameter_covariance(run.model, data, options.fit,
                                        select_for_fitting(data));
  run.g_variance = inverse_scale_variances(run.model, data, run.covariance);
  propagate_scale_variances(data, run.g, run.g_variance);
  step("the scale's uncertainty");
  run.outliers = reject_outliers(data, run.g);
  if (options.model_variances)
    model_variances(data, run.g);
  run.error_model = refine_error_model(data, run.g);
  apply_error_model(data, run.error_model);
  step("final outliers and error model");
  return run;
}

void write_scaling(Table &reflections, const ScaleData &data,
                   const std::vector<double> &g,
                   const std::vector<double> &g_variance) {
  const std::size_t n = reflections.nrows;
  Column &isf = reflections.real_column("inverse_scale_factor", "double", 1);
  Column &isf_var =
      reflections.real_column("inverse_scale_factor_variance", "double", 1);
  Column &value = reflections.real_column("intensity.scale.value", "double", 1);
  Column &variance =
      reflections.real_column("intensity.scale.variance", "double", 1);
  // A copy of the flags, changed and set back. int_column() makes a new column
  // of zeroes: used here, it wiped every flag integration had set.
  Column flags = reflections.at("flags");
  std::vector<bool> scaled(n, false);
  for (std::size_t r = 0; r < n; ++r) {
    isf.reals[r] = 1.0;
    isf_var.reals[r] = 0.0;
    value.reals[r] = 0.0;
    variance.reals[r] = 0.0;
    flags.ints[r] &=
        ~(flag::kOutlierInScaling | flag::kExcludedForScaling | flag::kScaled);
  }
  for (std::size_t i = 0; i < data.size(); ++i) {
    const std::size_t r = data.row[i];
    scaled[r] = true;
    isf.reals[r] = g[i];
    isf_var.reals[r] = i < g_variance.size() ? g_variance[i] : 0.0;
    value.reals[r] = data.intensity[i];
    variance.reals[r] = data.variance[i];
    flags.ints[r] |= data.outlier[i] ? flag::kOutlierInScaling : flag::kScaled;
  }
  for (std::size_t r = 0; r < n; ++r)
    if (!scaled[r])
      flags.ints[r] |= flag::kExcludedForScaling;
  reflections.set("flags", std::move(flags));
}

ParameterCovariance
parameter_covariance(const ScaleModel &model, const ScaleData &data,
                     const ScaleFitOptions &options,
                     const std::vector<std::size_t> &given) {
  ParameterCovariance out;
  std::vector<std::size_t> use = given;
  if (use.empty())
    for (std::size_t i = 0; i < data.size(); ++i)
      if (!data.outlier[i])
        use.push_back(i);
  const std::size_t n = model.size();
  if (use.empty() || n == 0)
    return out;
  std::vector<std::vector<std::size_t>> members(data.unique.size());
  for (std::size_t i : use)
    members[data.group[i]].push_back(i);
  std::vector<double> g, merged;
  const double phi = target(model, data, members, use, options, &g, &merged);
  std::vector<double> N, b;
  normal_equations(model, data, members, g, merged, options, N, b);

  // Z: the directions keeping the scale's sum and the relative B's sum fixed,
  // e_i - e_last over the indices of each -- every sweep's scale points
  // together, and every sweep's decay points: one overall scale and one
  // overall B, not one a sweep; the absorption terms are free.
  std::vector<std::vector<double>> Z;
  const auto block = [&](const std::vector<std::size_t> &indices) {
    for (std::size_t i = 0; i + 1 < indices.size(); ++i) {
      std::vector<double> z(n, 0.0);
      z[indices[i]] = 1.0;
      z[indices.back()] = -1.0;
      Z.push_back(std::move(z));
    }
  };
  std::vector<std::size_t> scales, decays;
  for (std::size_t w = 0; w < model.sweeps(); ++w) {
    for (std::size_t i = 0; i < model.shape(w).scale_points; ++i)
      scales.push_back(model.first_scale(w) + i);
    for (std::size_t i = 0; i < model.shape(w).decay_points; ++i)
      decays.push_back(model.first_decay(w) + i);
  }
  block(scales);
  block(decays);
  for (std::size_t w = 0; w < model.sweeps(); ++w)
    for (std::size_t k = 0; k < harmonic_count(model.shape(w).lmax); ++k) {
      std::vector<double> z(n, 0.0);
      z[model.first_absorption(w) + k] = 1.0;
      Z.push_back(std::move(z));
    }
  const std::size_t m = Z.size();
  if (m == 0)
    return out;
  std::vector<double> NZ(n * m, 0.0), M(m * m, 0.0);
  for (std::size_t a = 0; a < n; ++a)
    for (std::size_t j = 0; j < m; ++j)
      for (std::size_t c = 0; c < n; ++c)
        NZ[a * m + j] += N[a * n + c] * Z[j][c];
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < m; ++j)
      for (std::size_t a = 0; a < n; ++a)
        M[i * m + j] += Z[i][a] * NZ[a * m + j];
  // M's inverse from one Cholesky factorisation, a column at a time from the
  // factor: a factorisation a column, as it was, is m^4 / 3 -- 446
  // constrained parameters on Graeme's 36000 images, 10^10 operations, the
  // 3 s of 'the scale's uncertainty'.
  std::vector<double> L = M;
  for (std::size_t j = 0; j < m; ++j) {
    double d = L[j * m + j];
    for (std::size_t k = 0; k < j; ++k)
      d -= L[j * m + k] * L[j * m + k];
    if (!(d > 0.0))
      return out;
    const double root = std::sqrt(d);
    L[j * m + j] = root;
    for (std::size_t i = j + 1; i < m; ++i) {
      double v = L[i * m + j];
      for (std::size_t k = 0; k < j; ++k)
        v -= L[i * m + k] * L[j * m + k];
      L[i * m + j] = v / root;
    }
  }
  std::vector<double> Minv(m * m, 0.0);
  for_each_index(m, [&](std::size_t j) {
    std::vector<double> x(m, 0.0);
    x[j] = 1.0;
    for (std::size_t i = 0; i < m; ++i) { // L y = e_j
      double v = x[i];
      for (std::size_t k = 0; k < i; ++k)
        v -= L[i * m + k] * x[k];
      x[i] = v / L[i * m + i];
    }
    for (std::size_t i = m; i-- > 0;) { // L^T x = y
      double v = x[i];
      for (std::size_t k = i + 1; k < m; ++k)
        v -= L[k * m + i] * x[k];
      x[i] = v / L[i * m + i];
    }
    for (std::size_t i = 0; i < m; ++i)
      Minv[i * m + j] = x[i];
  });
  // The degrees of freedom: observations, less the merged intensities, less
  // the free parameters.
  std::size_t groups = 0;
  for (const auto &group_members : members)
    if (!group_members.empty())
      ++groups;
  const std::size_t fitted = groups + m;
  out.degrees_of_freedom = use.size() > fitted ? use.size() - fitted : 1;
  const double restraint = restraint_of(model, options);
  out.goodness_of_fit =
      (phi - restraint) / static_cast<double>(out.degrees_of_freedom);
  out.matrix.assign(n * n, 0.0);
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < m; ++j) {
      const double v = out.goodness_of_fit * Minv[i * m + j];
      if (v == 0.0)
        continue;
      for (std::size_t a = 0; a < n; ++a) {
        if (Z[i][a] == 0.0)
          continue;
        for (std::size_t c = 0; c < n; ++c)
          out.matrix[a * n + c] += Z[i][a] * v * Z[j][c];
      }
    }
  out.ok = true;
  return out;
}

std::vector<double>
inverse_scale_variances(const ScaleModel &model, const ScaleData &data,
                        const ParameterCovariance &covariance) {
  std::vector<double> out(data.size(), 0.0);
  if (!covariance.ok)
    return out;
  const std::size_t n = model.size();
  // A block of observations a thread, one gradient buffer each: a vector made
  // and freed for every observation, 6.7 million on Graeme's 36000 images,
  // was most of the 3 s this took.
  const std::size_t count = data.size();
  const std::size_t blocks =
      std::min<std::size_t>(std::max<std::size_t>(count, 1), kParallelBlocks);
  for_each_block(blocks, [&](std::size_t blk) {
    std::vector<std::pair<std::size_t, double>> grad;
    for (std::size_t i = block_begin(count, blocks, blk);
         i < block_begin(count, blocks, blk + 1); ++i) {
      grad.clear();
      model.inverse_scale(data.observation[i], &grad);
      double v = 0.0;
      for (const auto &[a, da] : grad)
        for (const auto &[c, dc] : grad)
          v += da * covariance.matrix[a * n + c] * dc;
      out[i] = std::fmax(v, 0.0);
    }
  });
  return out;
}

void propagate_scale_variances(ScaleData &data, const std::vector<double> &g,
                               const std::vector<double> &g_variance) {
  data.scale_term.resize(data.size(), 0.0);
  for (std::size_t i = 0; i < data.size(); ++i)
    data.scale_term[i] = g[i] > 0.0 ? data.intensity[i] * data.intensity[i] *
                                          g_variance[i] / (g[i] * g[i])
                                    : 0.0;
}

} // namespace mxi
