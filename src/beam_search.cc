// The beam centre found from the strong spots: see beam_search.hh.
#include "beam_search.hh"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <random>
#include <vector>

#include "index.hh"
#include "parallel.hh"
#include "portable_random.hh"

namespace mxi {

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;

// One periodicity measured on the spots: the vector t refitted to them, and
// S = mean exp(2 pi i r.t) there. Refitted at every trial, as DIALS finds
// each direction's period again by a one-dimensional FFT -- but the whole
// vector, not its length alone: the FFT's vectors are only as fine as its
// grid, and a direction slightly off scatters the phases, which a small move
// of the detector partly mends, so that even on an exact lattice the best
// score lay 0.1 mm off the truth. Three steps of Gauss-Newton: the spots'
// phases about their mean phase, those within a quarter turn of it, the
// change in t that flattens them by least squares. On an exact lattice at the
// true centre t comes to the lattice's vector and the mean phase to zero.
std::complex<double> measure(const std::vector<Vec3> &points, Vec3 t) {
  if (points.empty() || !(t.norm() > 0.0))
    return 0.0;
  const double n = static_cast<double>(points.size());
  std::complex<double> sum = 0.0;
  for (std::int32_t step = 0; step < 4; ++step) {
    sum = 0.0;
    for (const Vec3 &r : points)
      sum += std::polar(1.0, kTwoPi * r.dot(t));
    if (step == 3 || std::abs(sum) < 1e-9 * n)
      break;
    const double mean = std::arg(sum);
    // Normal equations for dt: sum r r^T dt = -sum r psi / (2 pi).
    double a[6] = {0, 0, 0, 0, 0, 0}, b[3] = {0, 0, 0};
    std::size_t used = 0;
    for (const Vec3 &r : points) {
      double psi = kTwoPi * r.dot(t) - mean;
      psi -= kTwoPi * std::round(psi / kTwoPi);
      if (std::abs(psi) > 0.5 * 3.14159265358979323846)
        continue;
      const double y = -psi / kTwoPi;
      a[0] += r.x * r.x, a[1] += r.x * r.y, a[2] += r.x * r.z;
      a[3] += r.y * r.y, a[4] += r.y * r.z, a[5] += r.z * r.z;
      b[0] += r.x * y, b[1] += r.y * y, b[2] += r.z * y;
      ++used;
    }
    if (used < 10)
      break;
    const Mat3 m{a[0], a[1], a[2], a[1], a[3], a[4], a[2], a[4], a[5]};
    bool ok = false;
    const Mat3 inverse = m.inverse(&ok);
    if (!ok)
      break;
    t = t + inverse * Vec3{b[0], b[1], b[2]};
  }
  return sum / n;
}

// The score: over the periodicities, Re S, |S| cos of its phase -- with the
// centre right every spot projects to an integer, every phase zero, and the
// lattice undistorted, every |S| at its largest. DIALS sums cos of the phase
// alone; with |S| too, a centre that brings the phases back by shifting the
// lattice one of its own vectors -- an alias, distorted by the Ewald sphere's
// curvature -- scores lower: on 30 degrees of insulin started 1.8 mm out,
// cos alone found the alias, this the truth.
double score(const std::vector<Vec3> &points, const std::vector<Vec3> &periods) {
  double total = 0.0;
  for (const Vec3 &t : periods) {
    const std::complex<double> s = measure(points, t);
    total += s.real();
  }
  return total;
}

} // namespace

bool beam_centre(const Experiment &e, std::size_t *panel, double mm[2],
                 double px[2]) {
  const Vec3 s0 = e.beam.s0().normalized();
  for (std::size_t p = 0; p < e.detector.panels.size(); ++p) {
    const Panel &q = e.detector.panels[p];
    const Vec3 normal = q.fast.cross(q.slow);
    const double along = s0.dot(normal);
    if (std::abs(along) < 1e-12)
      continue;
    const double t = q.origin.dot(normal) / along;
    if (!(t > 0.0))
      continue;
    const Vec3 on = s0 * t - q.origin;
    const double x = on.dot(q.fast.normalized()),
                 y = on.dot(q.slow.normalized());
    const double fx = x / q.pixel_size[0], fy = y / q.pixel_size[1];
    if (fx < 0.0 || fy < 0.0 || fx > static_cast<double>(q.image_size[0]) ||
        fy > static_cast<double>(q.image_size[1]))
      continue;
    *panel = p;
    mm[0] = x;
    mm[1] = y;
    px[0] = fx;
    px[1] = fy;
    return true;
  }
  return false;
}

void move_detectors(ExperimentList &experiments, const Vec3 &offset) {
  for (Experiment &e : experiments.experiments)
    for (Panel &p : e.detector.panels)
      p.origin = p.origin + offset;
}

BeamSearchResult search_beam(ExperimentList &experiments, const Table &strong,
                             const BeamSearchOptions &options) {
  BeamSearchResult out;
  if (experiments.experiments.empty()) {
    out.why = "no experiments";
    return out;
  }
  const Experiment &first = experiments.experiments[0];
  if (!beam_centre(first, &out.panel_before, out.mm_before, out.px_before)) {
    out.why = "the beam does not meet the detector";
    return out;
  }

  // The spots: those coarser than --d-min, then a random subset of at most
  // max_reflections, as DIALS takes them.
  std::vector<Vec3> all = reciprocal_lattice_points(experiments, strong);
  std::vector<std::size_t> rows;
  for (std::size_t r = 0; r < all.size(); ++r) {
    const double n = all[r].norm();
    if (n > 0.0 && (options.d_min <= 0.0 || 1.0 / n > options.d_min))
      rows.push_back(r);
  }
  if (rows.size() > options.max_reflections) {
    std::mt19937 random(options.seed);
    portable_shuffle(rows.begin(), rows.end(), random);
    rows.resize(options.max_reflections);
    std::sort(rows.begin(), rows.end());
  }
  if (rows.size() < 50) {
    out.why = "too few spots (" + std::to_string(rows.size()) + ")";
    return out;
  }
  const Table spots = select_rows(strong, rows);
  std::vector<Vec3> points;
  for (std::size_t r : rows)
    points.push_back(all[r]);
  out.spots = rows.size();

  // The periodicities: the three-dimensional FFT's strongest candidate
  // vectors, as indexing finds them, from the centre as given.
  IndexOptions index_options;
  index_options.max_cell = options.max_cell;
  SearchScale scale;
  if (!search_scale(experiments, spots, points, index_options, &scale)) {
    out.why = "no largest cell or resolution could be had from the spots";
    return out;
  }
  out.max_cell = scale.max_cell;
  out.d_min = scale.d_min;
  // The FFT's candidates are padded with peaks indexing would throw away --
  // on 30 degrees of cubic insulin 4 of the strongest 20 were the lattice's,
  // |S|/N 0.95 at the true centre, the rest 0.01, and two short ones at 0.07
  // and 0.13 pointing the wrong way, which left the true centre scoring 2. So
  // only the coherent: |S| at least 40 per cent of the best's, from the centre
  // as given (where the lattice's own stayed 0.64 to 0.91 with the centre 4
  // mm out, the rest 0.12 at most); with their sums and differences, lattice
  // vectors too, for more directions -- up to options.candidates, the most
  // coherent first, none within 5 degrees of a better.
  const std::vector<Vec3> found = find_candidate_vectors(
      points, index_options, scale.d_min, scale.max_cell, scale.grid);
  std::vector<std::pair<double, Vec3>> scored;
  double strongest = 0.0;
  for (const Vec3 &t : found) {
    const double c = std::abs(measure(points, t));
    scored.emplace_back(c, t);
    strongest = std::max(strongest, c);
  }
  std::vector<Vec3> kept;
  for (const auto &[c, t] : scored)
    if (c >= 0.4 * strongest)
      kept.push_back(t);
  const std::size_t base = kept.size();
  for (std::size_t i = 0; i < base; ++i)
    for (std::size_t j = i + 1; j < base; ++j)
      for (const Vec3 &t : {kept[i] + kept[j], kept[i] - kept[j]})
        if (t.norm() <= 1.5 * scale.max_cell && t.norm() > 2.0 * scale.d_min) {
          const double c = std::abs(measure(points, t));
          if (c >= 0.4 * strongest)
            scored.emplace_back(c, t);
        }
  std::sort(scored.begin(), scored.end(),
            [](const auto &a, const auto &b) { return a.first > b.first; });
  std::vector<Vec3> periods;
  const double parallel = std::cos(Scan::radians(5.0));
  for (const auto &[c, t] : scored) {
    if (c < 0.4 * strongest || periods.size() >= options.candidates)
      continue;
    bool near = false;
    for (const Vec3 &q : periods)
      near = near ||
             std::abs(q.dot(t)) >= parallel * q.norm() * t.norm();
    if (!near)
      periods.push_back(t);
  }
  out.candidates = periods.size();
  if (periods.size() < 3) {
    out.why = "too few periodicities in the spots (" +
              std::to_string(periods.size()) + ")";
    return out;
  }

  // A trial: the detector moved across the beam by (u, v) mm, the spots
  // mapped again, scored.
  const Vec3 s0 = first.beam.s0().normalized();
  const Vec3 across0 = s0.cross(Vec3{1.0, 0.0, 0.0}).normalized();
  const Vec3 across1 = across0.cross(s0).normalized();
  const Vec3 across2 = across1.cross(s0).normalized();
  const auto trial = [&](double u, double v) {
    ExperimentList moved = experiments;
    move_detectors(moved, across1 * u + across2 * v);
    return score(reciprocal_lattice_points(moved, spots), periods);
  };
  out.score_start = trial(0.0, 0.0);

  // The grid: the scope either way, in steps of the pixel times the binning.
  const double step =
      first.detector.panels[out.panel_before].pixel_size[0] * options.binning;
  const std::int32_t reach = std::max<std::int32_t>(
      1, static_cast<std::int32_t>(options.scope_mm / step));
  const std::int32_t side = 2 * reach + 1;
  std::vector<double> scores(static_cast<std::size_t>(side) *
                             static_cast<std::size_t>(side));
  for_each_index(scores.size(), [&](std::size_t i) {
    const std::int32_t x = static_cast<std::int32_t>(i) % side - reach;
    const std::int32_t y = static_cast<std::int32_t>(i) / side - reach;
    scores[i] = trial(x * step, y * step);
  });
  out.trials = scores.size();
  const double best = *std::max_element(scores.begin(), scores.end());
  if (!(best > 0.0)) {
    out.why = "no trial centre scored above zero";
    return out;
  }
  // Of the scores within 90 per cent of the best, DIALS takes the nearest
  // the centre as given. But a detector moved by the right amount shifts the
  // reciprocal lattice by one of its own vectors, every spot on a lattice
  // point again: on 30 degrees of insulin started 24 px out, the grid's
  // best was such an alias 3.1 mm from the truth -- the shortest reciprocal
  // vector, 0.018 A^-1, at 170 mm -- nearer the start than the truth was.
  // The Ewald sphere's curvature makes the shift uneven, so the truth scores
  // higher once refined. So every local maximum within 90 per cent of the
  // best, at most five, refined; the highest refined score taken, the
  // nearest the start only among those within half a per cent of it.
  std::vector<std::pair<double, std::size_t>> peaks;
  for (std::size_t i = 0; i < scores.size(); ++i) {
    if (scores[i] < 0.9 * best)
      continue;
    const std::int32_t x = static_cast<std::int32_t>(i) % side,
                       y = static_cast<std::int32_t>(i) / side;
    bool peak = true;
    for (std::int32_t dy = -1; dy <= 1 && peak; ++dy)
      for (std::int32_t dx = -1; dx <= 1 && peak; ++dx) {
        const std::int32_t xx = x + dx, yy = y + dy;
        if ((dx || dy) && xx >= 0 && yy >= 0 && xx < side && yy < side)
          peak = scores[static_cast<std::size_t>(yy * side + xx)] <= scores[i];
      }
    if (peak)
      peaks.emplace_back(scores[i], i);
  }
  std::sort(peaks.begin(), peaks.end(),
            [](const auto &a, const auto &b) { return a.first > b.first; });
  if (peaks.size() > 5)
    peaks.resize(5);
  // A pattern search from a point: its step halved from the grid's down to
  // a thousandth of a millimetre, moving to the best of the eight neighbours
  // while one is better.
  const auto refine = [&](double &pu, double &pv) {
    double at = trial(pu, pv);
    for (double h = step / 2.0; h >= 0.001; h /= 2.0) {
      for (std::int32_t moves = 0; moves < 50; ++moves) {
        double bu = pu, bv = pv, bs = at;
        std::vector<double> around(8);
        const std::int32_t dx[8] = {1, -1, 0, 0, 1, 1, -1, -1};
        const std::int32_t dy[8] = {0, 0, 1, -1, 1, -1, 1, -1};
        for_each_index(8, [&](std::size_t k) {
          around[k] = trial(pu + dx[k] * h, pv + dy[k] * h);
        });
        out.trials += 8;
        for (std::size_t k = 0; k < 8; ++k)
          if (around[k] > bs)
            bs = around[k], bu = pu + dx[k] * h, bv = pv + dy[k] * h;
        if (bs <= at)
          break;
        pu = bu;
        pv = bv;
        at = bs;
      }
    }
    return at;
  };
  struct Refined {
    double score, u, v;
  };
  std::vector<Refined> refined;
  for (const auto &[s0_unused, i] : peaks) {
    double pu = (static_cast<std::int32_t>(i) % side - reach) * step;
    double pv = (static_cast<std::int32_t>(i) / side - reach) * step;
    const double at = refine(pu, pv);
    refined.push_back({at, pu, pv});
  }
  double top = -1e300;
  for (const Refined &r : refined)
    top = std::max(top, r.score);
  double u = 0.0, v = 0.0, here = top, nearest = 1e300;
  for (const Refined &r : refined)
    if (r.score >= top - 0.005 * std::abs(top) &&
        r.u * r.u + r.v * r.v < nearest)
      nearest = r.u * r.u + r.v * r.v, u = r.u, v = r.v, here = r.score;
  out.peaks = refined.size();
  out.score_best = here;
  out.offset = across1 * u + across2 * v;
  move_detectors(experiments, out.offset);
  if (!beam_centre(experiments.experiments[0], &out.panel_after, out.mm_after,
                   out.px_after)) {
    move_detectors(experiments, out.offset * -1.0);
    out.why = "the centre found is off the detector";
    return out;
  }
  out.found = true;
  return out;
}

} // namespace mxi
