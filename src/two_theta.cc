// The unit cell from two-theta: see two_theta.hh.

#include "two_theta.hh"

#include <algorithm>
#include <cmath>

#include "linalg.hh"

namespace mxi {

namespace {

constexpr double kDeg = M_PI / 180.0;

// 1/d^2 for hkl in a cell: h^T G* h, G* the inverse of the real-space metric.
double inverse_d2(const UnitCell &c, const Miller &h) {
  const double ca = std::cos(c.alpha * kDeg), cb = std::cos(c.beta * kDeg),
               cg = std::cos(c.gamma * kDeg);
  const Mat3 G{c.a * c.a,      c.a * c.b * cg, c.a * c.c * cb,
               c.a * c.b * cg, c.b * c.b,      c.b * c.c * ca,
               c.a * c.c * cb, c.b * c.c * ca, c.c * c.c};
  const Mat3 S = G.inverse();
  const double x[3] = {static_cast<double>(h[0]), static_cast<double>(h[1]),
                       static_cast<double>(h[2])};
  double out = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out += x[i] * S(i, j) * x[j];
  return out;
}

// The parameters a crystal system refines, and the cell they make. A
// monoclinic cell's unique angle and a trigonal cell's setting are read from
// the cell it starts from.
struct Parametrisation {
  CrystalSystem system;
  int unique = 1;            // monoclinic: 0 alpha, 1 beta, 2 gamma
  bool rhombohedral = false; // trigonal on rhombohedral axes
  // For each cell parameter, which refined parameter it is (-1: fixed).
  std::array<int, 6> source{};
  std::array<double, 6> fixed{};

  std::vector<double> start(const UnitCell &c) const {
    const std::array<double, 6> v = {c.a, c.b, c.c, c.alpha, c.beta, c.gamma};
    int n = 0;
    for (int s : source)
      n = std::max(n, s + 1);
    std::vector<double> p(static_cast<std::size_t>(n), 0.0);
    std::vector<int> seen(static_cast<std::size_t>(n), 0);
    for (int k = 0; k < 6; ++k)
      if (source[static_cast<std::size_t>(k)] >= 0) {
        const auto s =
            static_cast<std::size_t>(source[static_cast<std::size_t>(k)]);
        p[s] += v[static_cast<std::size_t>(k)];
        ++seen[s];
      }
    for (std::size_t s = 0; s < p.size(); ++s)
      p[s] /= std::max(seen[s], 1);
    return p;
  }

  UnitCell cell(const std::vector<double> &p) const {
    std::array<double, 6> v{};
    for (std::size_t k = 0; k < 6; ++k)
      v[k] = source[k] >= 0 ? p[static_cast<std::size_t>(source[k])] : fixed[k];
    return {v[0], v[1], v[2], v[3], v[4], v[5]};
  }
};

Parametrisation parametrise(CrystalSystem system, const UnitCell &c) {
  Parametrisation out;
  out.system = system;
  out.fixed = {0, 0, 0, 90, 90, 90};
  switch (system) {
  case CrystalSystem::cubic:
    out.source = {0, 0, 0, -1, -1, -1};
    break;
  case CrystalSystem::tetragonal:
    out.source = {0, 0, 1, -1, -1, -1};
    break;
  case CrystalSystem::hexagonal:
    out.source = {0, 0, 1, -1, -1, -1};
    out.fixed[5] = 120.0;
    break;
  case CrystalSystem::trigonal:
    if (std::fabs(c.gamma - 120.0) < 1.0) {
      out.source = {0, 0, 1, -1, -1, -1};
      out.fixed[5] = 120.0;
    } else {
      out.rhombohedral = true;
      out.source = {0, 0, 0, 1, 1, 1};
    }
    break;
  case CrystalSystem::orthorhombic:
    out.source = {0, 1, 2, -1, -1, -1};
    break;
  case CrystalSystem::monoclinic: {
    const double d[3] = {std::fabs(c.alpha - 90.0), std::fabs(c.beta - 90.0),
                         std::fabs(c.gamma - 90.0)};
    out.unique = static_cast<int>(std::max_element(d, d + 3) - d);
    out.source = {0, 1, 2, -1, -1, -1};
    out.source[static_cast<std::size_t>(3 + out.unique)] = 3;
    break;
  }
  case CrystalSystem::triclinic:
    out.source = {0, 1, 2, 3, 4, 5};
    break;
  }
  return out;
}

double two_theta_of(const UnitCell &c, const TwoThetaObservation &o) {
  const double s =
      o.wavelength * std::sqrt(std::fmax(inverse_d2(c, o.hkl), 0.0)) / 2.0;
  return s < 1.0 ? 2.0 * std::asin(s) : M_PI;
}

} // namespace

CrystalSystem crystal_system_of(int n) {
  if (n <= 2)
    return CrystalSystem::triclinic;
  if (n <= 15)
    return CrystalSystem::monoclinic;
  if (n <= 74)
    return CrystalSystem::orthorhombic;
  if (n <= 142)
    return CrystalSystem::tetragonal;
  if (n <= 167)
    return CrystalSystem::trigonal;
  if (n <= 194)
    return CrystalSystem::hexagonal;
  return CrystalSystem::cubic;
}

const char *name_of(CrystalSystem system) {
  switch (system) {
  case CrystalSystem::triclinic:
    return "triclinic";
  case CrystalSystem::monoclinic:
    return "monoclinic";
  case CrystalSystem::orthorhombic:
    return "orthorhombic";
  case CrystalSystem::tetragonal:
    return "tetragonal";
  case CrystalSystem::trigonal:
    return "trigonal";
  case CrystalSystem::hexagonal:
    return "hexagonal";
  case CrystalSystem::cubic:
    return "cubic";
  }
  return "?";
}

TwoThetaResult refine_two_theta(const std::vector<TwoThetaObservation> &obs,
                                const UnitCell &start, CrystalSystem system) {
  TwoThetaResult out;
  out.system = system;
  const Parametrisation param = parametrise(system, start);
  std::vector<double> p = param.start(start);
  const std::size_t np = p.size(), n = obs.size();
  out.parameters = static_cast<int>(np);
  out.reflections = n;
  if (n <= np + 1) {
    out.why = "too few reflections";
    return out;
  }
  std::vector<double> r(n), J(n * np);
  double chi2 = 0.0;
  // Gauss-Newton: the derivatives by central differences, cheap at a handful
  // of parameters.
  for (int iteration = 0; iteration < 30; ++iteration) {
    const UnitCell c = param.cell(p);
    for (std::size_t i = 0; i < n; ++i)
      r[i] = obs[i].two_theta - two_theta_of(c, obs[i]);
    for (std::size_t k = 0; k < np; ++k) {
      const double h = 1e-6 * std::max(std::fabs(p[k]), 1.0);
      std::vector<double> up = p, down = p;
      up[k] += h;
      down[k] -= h;
      const UnitCell cu = param.cell(up), cd = param.cell(down);
      for (std::size_t i = 0; i < n; ++i)
        J[i * np + k] =
            (two_theta_of(cu, obs[i]) - two_theta_of(cd, obs[i])) / (2.0 * h);
    }
    std::vector<double> N(np * np, 0.0), g(np, 0.0);
    chi2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const double w = 1.0 / obs[i].variance;
      chi2 += w * r[i] * r[i];
      for (std::size_t a = 0; a < np; ++a) {
        g[a] += w * J[i * np + a] * r[i];
        for (std::size_t b = 0; b < np; ++b)
          N[a * np + b] += w * J[i * np + a] * J[i * np + b];
      }
    }
    std::vector<double> step = g, M = N;
    if (!solve_spd(M.data(), step.data(), np)) {
      out.why = "the normal matrix is singular";
      return out;
    }
    double largest = 0.0;
    for (std::size_t k = 0; k < np; ++k) {
      p[k] += step[k];
      largest = std::max(largest,
                         std::fabs(step[k]) / std::max(std::fabs(p[k]), 1.0));
    }
    if (largest < 1e-10)
      break;
  }
  // The covariance at the solution.
  const UnitCell c = param.cell(p);
  std::vector<double> N(np * np, 0.0);
  chi2 = 0.0;
  double sum2 = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    r[i] = obs[i].two_theta - two_theta_of(c, obs[i]);
    const double w = 1.0 / obs[i].variance;
    chi2 += w * r[i] * r[i];
    sum2 += r[i] * r[i];
    for (std::size_t a = 0; a < np; ++a)
      for (std::size_t b = 0; b < np; ++b)
        N[a * np + b] += w * J[i * np + a] * J[i * np + b];
  }
  const double scale = chi2 / static_cast<double>(n - np);
  std::vector<double> esd_p(np, 0.0);
  for (std::size_t k = 0; k < np; ++k) {
    std::vector<double> M = N, e(np, 0.0);
    e[k] = 1.0;
    if (!solve_spd(M.data(), e.data(), np)) {
      out.why = "the normal matrix is singular";
      return out;
    }
    esd_p[k] = std::sqrt(std::fmax(e[k] * scale, 0.0));
  }
  out.cell = c;
  for (std::size_t k = 0; k < 6; ++k)
    out.esd[k] = param.source[k] >= 0
                     ? esd_p[static_cast<std::size_t>(param.source[k])]
                     : 0.0;
  out.rmsd_deg = std::sqrt(sum2 / static_cast<double>(n)) / kDeg;
  out.valid = true;
  return out;
}

TwoThetaResult two_theta_refine(const ExperimentList &experiments,
                                const Table &reflections,
                                const SpaceGroup &group, double min_isigma) {
  TwoThetaResult none;
  const char *need[] = {"miller_index",
                        "id",
                        "flags",
                        "panel",
                        "xyzobs.mm.value",
                        "xyzobs.mm.variance",
                        "intensity.scale.value",
                        "intensity.scale.variance"};
  for (const char *c : need)
    if (!reflections.has(c)) {
      none.why = std::string("no ") + c + " column";
      return none;
    }
  const Column &hkl = reflections.at("miller_index"),
               &id = reflections.at("id"), &flags = reflections.at("flags"),
               &panel = reflections.at("panel"),
               &mm = reflections.at("xyzobs.mm.value"),
               &mmv = reflections.at("xyzobs.mm.variance"),
               &I = reflections.at("intensity.scale.value"),
               &V = reflections.at("intensity.scale.variance");
  const UnitCell *start = nullptr;
  UnitCell first;
  for (std::size_t e = 0; e < experiments.size(); ++e)
    if (experiments[e].crystal) {
      first = experiments[e].crystal->cell();
      start = &first;
      break;
    }
  if (!start) {
    none.why = "no crystal";
    return none;
  }
  std::vector<TwoThetaObservation> obs;
  for (std::size_t r = 0; r < reflections.nrows; ++r) {
    const std::int64_t f = flags.integer(r);
    if (((f >> 26) & 1) == 0 || ((f >> 23) & 1) != 0) // scaled, not an outlier
      continue;
    if (!(V.reals[r] > 0.0) || I.reals[r] / std::sqrt(V.reals[r]) < min_isigma)
      continue;
    const auto e = static_cast<std::size_t>(id.integer(r));
    if (e >= experiments.size() || !experiments[e].crystal)
      continue;
    const Experiment &x = experiments[e];
    const auto p = static_cast<std::size_t>(panel.integer(r));
    if (p >= x.detector.size())
      continue;
    const double fx = mm.reals[r * 3], fy = mm.reals[r * 3 + 1];
    // Each centroid's variance at least a pixel's own, size^2 / 12: counting
    // statistics alone put strong spots' at a tenth of a micron, and those,
    // mostly at low resolution where an angle says least about the cell,
    // carried 43 per cent of the weight -- insulin's a 77.931 A with chi^2 89
    // times its degrees of freedom; floored, 77.911, the median of each
    // reflection's own, and 2.6.
    const double fpx = x.detector[p].pixel_size[0],
                 fpy = x.detector[p].pixel_size[1];
    const double vx = std::fmax(mmv.reals[r * 3], fpx * fpx / 12.0),
                 vy = std::fmax(mmv.reals[r * 3 + 1], fpy * fpy / 12.0);
    if (!std::isfinite(fx) || !std::isfinite(fy) || !(mmv.reals[r * 3] > 0.0) ||
        !(mmv.reals[r * 3 + 1] > 0.0))
      continue;
    const Vec3 s0 = x.beam.s0();
    const auto angle = [&](double ux, double uy) {
      const Vec3 lab = x.detector[p].lab_coord_mm(ux, uy);
      const double c = lab.dot(s0) / (lab.norm() * s0.norm());
      return std::acos(std::fmax(-1.0, std::fmin(1.0, c)));
    };
    const double t = angle(fx, fy);
    const double h = 1e-3; // mm
    const double dx = (angle(fx + h, fy) - angle(fx - h, fy)) / (2.0 * h);
    const double dy = (angle(fx, fy + h) - angle(fx, fy - h)) / (2.0 * h);
    const double var = dx * dx * vx + dy * dy * vy;
    if (!(var > 0.0))
      continue;
    TwoThetaObservation o;
    for (int k = 0; k < 3; ++k)
      o.hkl[static_cast<std::size_t>(k)] =
          static_cast<int>(hkl.ints[r * 3 + static_cast<std::size_t>(k)]);
    o.two_theta = t;
    o.variance = var;
    o.wavelength = x.beam.wavelength;
    obs.push_back(o);
  }
  return refine_two_theta(obs, *start, crystal_system_of(group.number()));
}

} // namespace mxi
