// The cell from two-theta (src/two_theta.hh): planted cells recovered, and
// the e.s.d.s the scatter of the answer over many trials.

#include <cmath>
#include <random>
#include <vector>

#include "../src/two_theta.hh"
#include "check.hh"

namespace mxi {

namespace {

// Every reflection to 2 A of `cell`, its exact two-theta and Gaussian noise of
// `sigma` radians on it.
std::vector<TwoThetaObservation> planted(const UnitCell &cell, double sigma,
                                         std::mt19937 &rng) {
  std::normal_distribution<double> noise(0.0, sigma);
  std::vector<TwoThetaObservation> out;
  const double wavelength = 0.9537, kDeg = M_PI / 180.0;
  const double ca = std::cos(cell.alpha * kDeg),
               cb = std::cos(cell.beta * kDeg),
               cg = std::cos(cell.gamma * kDeg);
  const Mat3 G{cell.a * cell.a,      cell.a * cell.b * cg, cell.a * cell.c * cb,
               cell.a * cell.b * cg, cell.b * cell.b,      cell.b * cell.c * ca,
               cell.a * cell.c * cb, cell.b * cell.c * ca, cell.c * cell.c};
  const Mat3 S = G.inverse();
  for (int h = -12; h <= 12; ++h)
    for (int k = -12; k <= 12; ++k)
      for (int l = 0; l <= 12; ++l) {
        const double x[3] = {double(h), double(k), double(l)};
        double q = 0.0;
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            q += x[i] * S(i, j) * x[j];
        if (q <= 0.0 || q > 0.25) // to 2 A
          continue;
        TwoThetaObservation o;
        o.hkl = {h, k, l};
        o.wavelength = wavelength;
        o.two_theta =
            2.0 * std::asin(wavelength * std::sqrt(q) / 2.0) + noise(rng);
        o.variance = sigma * sigma;
        out.push_back(o);
      }
  return out;
}

} // namespace

TEST(a_planted_cubic_cell_is_recovered_and_its_esd_is_its_scatter) {
  std::mt19937 rng(7);
  const UnitCell truth{78.0, 78.0, 78.0, 90.0, 90.0, 90.0};
  const UnitCell start{78.3, 78.3, 78.3, 90.0, 90.0, 90.0};
  std::vector<double> a, esd;
  for (int trial = 0; trial < 200; ++trial) {
    const TwoThetaResult r = refine_two_theta(planted(truth, 2e-4, rng), start,
                                              CrystalSystem::cubic);
    check::is_true(r.valid, "refined");
    check::equal(r.parameters, 1, "cubic: one parameter");
    a.push_back(r.cell.a);
    esd.push_back(r.esd[0]);
    if (trial == 0) {
      check::is_true(std::fabs(r.cell.a - 78.0) < 4.0 * r.esd[0],
                     "a within four e.s.d.s of the truth");
      check::is_true(r.cell.b == r.cell.a && r.cell.c == r.cell.a &&
                         r.cell.alpha == 90.0 && r.esd[3] == 0.0,
                     "the cubic constraints held");
    }
  }
  double mean = 0.0, scatter = 0.0, mean_esd = 0.0;
  for (double v : a)
    mean += v / a.size();
  for (double v : a)
    scatter += (v - mean) * (v - mean) / (a.size() - 1);
  for (double v : esd)
    mean_esd += v / esd.size();
  scatter = std::sqrt(scatter);
  check::is_true(scatter > 0.8 * mean_esd && scatter < 1.2 * mean_esd,
                 "the e.s.d. the scatter: " + std::to_string(mean_esd) +
                     " against " + std::to_string(scatter));
}

TEST(a_planted_monoclinic_cell_is_recovered_with_its_unique_angle) {
  std::mt19937 rng(11);
  const UnitCell truth{40.0, 50.0, 60.0, 90.0, 105.0, 90.0};
  const UnitCell start{40.1, 49.9, 60.2, 90.0, 104.7, 90.0};
  const TwoThetaResult r = refine_two_theta(planted(truth, 2e-4, rng), start,
                                            CrystalSystem::monoclinic);
  check::is_true(r.valid, "refined");
  check::equal(r.parameters, 4, "monoclinic: a, b, c and beta");
  const double v[6] = {r.cell.a,     r.cell.b,    r.cell.c,
                       r.cell.alpha, r.cell.beta, r.cell.gamma};
  const double t[6] = {40.0, 50.0, 60.0, 90.0, 105.0, 90.0};
  for (int k : {0, 1, 2, 4})
    check::is_true(std::fabs(v[k] - t[k]) <
                       4.0 * r.esd[static_cast<std::size_t>(k)],
                   "parameter " + std::to_string(k) + " within four e.s.d.s");
  check::is_true(r.cell.alpha == 90.0 && r.cell.gamma == 90.0 &&
                     r.esd[3] == 0.0 && r.esd[5] == 0.0,
                 "alpha and gamma held at 90");
}

} // namespace mxi
