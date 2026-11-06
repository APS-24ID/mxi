// fit_together (src/reference.hh): reflections whose foregrounds overlap,
// fitted together on their pixels. On boxes made here: a little detector
// volume with a flat background and Gaussian spots of known intensity,
// Poisson-sampled once, each box copying its pixels from it.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "../src/reference.hh"
#include "check.hh"

namespace mxi {

namespace {

struct Spot {
  double x, y, z, intensity;
};
const double kSx = 1.2, kSy = 1.2, kSz = 0.9; // sigmas, pixels and frames

double gauss(const Spot &s, double x, double y, double z) {
  const double dx = (x - s.x) / kSx, dy = (y - s.y) / kSy, dz = (z - s.z) / kSz;
  return std::exp(-0.5 * (dx * dx + dy * dy + dz * dz));
}

// The volume's counts: background plus every spot, Poisson.
struct Volume {
  std::int32_t nx = 40, ny = 30, nz = 12;
  std::vector<float> counts;
  double background = 5.0;
};

Volume sample(const std::vector<Spot> &spots, std::mt19937 &random) {
  Volume v;
  v.counts.resize(static_cast<std::size_t>(v.nx) * v.ny * v.nz);
  std::vector<double> norm;
  for (const Spot &s : spots) {
    double n = 0.0;
    for (int z = 0; z < v.nz; ++z)
      for (int y = 0; y < v.ny; ++y)
        for (int x = 0; x < v.nx; ++x)
          n += gauss(s, x + 0.5, y + 0.5, z + 0.5);
    norm.push_back(n);
  }
  for (int z = 0; z < v.nz; ++z)
    for (int y = 0; y < v.ny; ++y)
      for (int x = 0; x < v.nx; ++x) {
        double mean = v.background;
        for (std::size_t k = 0; k < spots.size(); ++k)
          mean += spots[k].intensity *
                  gauss(spots[k], x + 0.5, y + 0.5, z + 0.5) / norm[k];
        std::poisson_distribution<int> poisson(mean);
        v.counts[(static_cast<std::size_t>(z) * v.ny + y) * v.nx + x] =
            static_cast<float>(poisson(random));
      }
  return v;
}

// A spot's box, +-5 sigma, from the volume; foreground within 3 sigma; its
// profile the spot's Gaussian on the box, normalised over the whole volume as
// the counts were; the background the true one.
OverlapMember box_of(const Spot &s, const Volume &v, Shoebox *box) {
  box->panel = 0;
  box->bbox[0] = std::max(0, static_cast<int>(std::floor(s.x - 5 * kSx)));
  box->bbox[1] = std::min(v.nx, static_cast<int>(std::ceil(s.x + 5 * kSx)));
  box->bbox[2] = std::max(0, static_cast<int>(std::floor(s.y - 5 * kSy)));
  box->bbox[3] = std::min(v.ny, static_cast<int>(std::ceil(s.y + 5 * kSy)));
  box->bbox[4] = std::max(0, static_cast<int>(std::floor(s.z - 5 * kSz)));
  box->bbox[5] = std::min(v.nz, static_cast<int>(std::ceil(s.z + 5 * kSz)));
  const std::size_t n = box->size();
  box->data.assign(n, 0.0f);
  box->background.assign(n, static_cast<float>(v.background));
  box->mask.assign(n, 0);
  double norm = 0.0;
  for (int z = 0; z < v.nz; ++z)
    for (int y = 0; y < v.ny; ++y)
      for (int x = 0; x < v.nx; ++x)
        norm += gauss(s, x + 0.5, y + 0.5, z + 0.5);
  OverlapMember m;
  m.box = box;
  m.profile.assign(n, 0.0);
  for (std::int32_t z = box->bbox[4]; z < box->bbox[5]; ++z)
    for (std::int32_t y = box->bbox[2]; y < box->bbox[3]; ++y)
      for (std::int32_t x = box->bbox[0]; x < box->bbox[1]; ++x) {
        const std::size_t i =
            box->at(x - box->bbox[0], y - box->bbox[2], z - box->bbox[4]);
        box->data[i] =
            v.counts[(static_cast<std::size_t>(z) * v.ny + y) * v.nx + x];
        const double g = gauss(s, x + 0.5, y + 0.5, z + 0.5);
        m.profile[i] = g / norm;
        const bool fg = g > std::exp(-0.5 * 9.0); // within 3 sigma
        box->mask[i] = static_cast<std::uint8_t>(
            shoebox_mask::kValid |
            (fg ? shoebox_mask::kForeground : shoebox_mask::kBackground));
      }
  return m;
}

} // namespace

TEST(one_reflection_fitted_together_is_fitted_as_alone) {
  std::mt19937 random(3);
  const Spot s{20, 15, 6, 3000};
  const Volume v = sample({s}, random);
  Shoebox box;
  const OverlapMember m = box_of(s, v, &box);
  const ProfileFit alone = fit_on_pixels(box, m.profile);
  const std::vector<ProfileFit> together = fit_together({m});
  check::is_true(alone.valid && together.size() == 1 && together[0].valid,
                 "both fitted");
  check::close(together[0].intensity, alone.intensity, 1e-9 * alone.intensity,
               "intensity");
  check::close(together[0].variance, alone.variance, 1e-9 * alone.variance,
               "variance");
  check::close(together[0].variance_background, alone.variance_background,
               1e-9 * alone.variance_background, "variance with no signal");
  check::close(together[0].correlation, alone.correlation, 1e-9, "correlation");
  check::close(together[0].measured, alone.measured, 1e-12, "measured");
}

TEST(overlapping_spots_fitted_together_are_unbiased_and_calibrated) {
  // A weak spot beside a strong one, three sigma apart in x: their 3-sigma
  // foregrounds overlap. Alone, the weak one counts its neighbour's photons;
  // together, neither does.
  const Spot weak{17, 15, 6, 1000}, strong{20.6, 15, 6, 8000};
  std::mt19937 random(11);
  const int trials = 400;
  double alone_sum = 0, together_sum = 0, z_sum = 0, z2_sum = 0, strong_sum = 0,
         isolated_sum = 0, isolated_z2 = 0, expected = weak.intensity,
         expected_strong = strong.intensity;
  for (int trial = 0; trial < trials; ++trial) {
    const Volume v = sample({weak, strong}, random);
    Shoebox a, b;
    const OverlapMember ma = box_of(weak, v, &a), mb = box_of(strong, v, &b);
    const ProfileFit alone = fit_on_pixels(a, ma.profile);
    const std::vector<ProfileFit> together = fit_together({ma, mb});
    check::is_true(alone.valid && together[0].valid && together[1].valid,
                   "fitted");
    alone_sum += alone.intensity;
    together_sum += together[0].intensity;
    strong_sum += together[1].intensity;
    // The same spot with no neighbour, fitted alone: the baseline.
    const Volume lone = sample({weak}, random);
    Shoebox c;
    const OverlapMember mc = box_of(weak, lone, &c);
    const ProfileFit isolated = fit_on_pixels(c, mc.profile);
    isolated_sum += isolated.intensity;
    // What a fit reports: the scale times the profile over the foreground,
    // here a 3-sigma sphere holding 0.971 of the Gaussian -- not the whole.
    double fraction = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
      if (a.mask[i] & shoebox_mask::kForeground)
        fraction += ma.profile[i];
    expected = weak.intensity * fraction;
    double strong_fraction = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i)
      if (b.mask[i] & shoebox_mask::kForeground)
        strong_fraction += mb.profile[i];
    expected_strong = strong.intensity * strong_fraction;
    isolated_z2 +=
        std::pow(isolated.intensity - expected, 2) / isolated.variance;
    const double z =
        (together[0].intensity - expected) / std::sqrt(together[0].variance);
    z_sum += z;
    z2_sum += z * z;
  }
  const double alone_mean = alone_sum / trials,
               together_mean = together_sum / trials;
  std::printf("        isolated %.1f (z^2 %.2f), alone %.1f, together %.1f, z "
              "mean %.3f var %.3f\n",
              isolated_sum / trials, isolated_z2 / trials, alone_sum / trials,
              together_sum / trials, z_sum / trials,
              z2_sum / trials - (z_sum / trials) * (z_sum / trials));
  const double z_mean = z_sum / trials,
               z_var = z2_sum / trials - z_mean * z_mean;
  // Alone, more than three times too strong: its neighbour's photons. Together,
  // as the same spot is with no neighbour at all.
  check::is_true(alone_mean > 1.5 * expected,
                 "alone, biased: " + std::to_string(alone_mean));
  // 400 fits, each sigma some 60: the mean's standard error 3.
  check::close(together_mean, expected, 0.015 * expected,
               "together, unbiased: " + std::to_string(together_mean));
  check::close(together_mean, isolated_sum / trials, 0.015 * expected,
               "together as with no neighbour");
  check::close(strong_sum / trials, expected_strong, 0.01 * expected_strong,
               "the strong one");
  check::close(z_mean, 0.0, 0.15,
               "deviations over sigma, mean: " + std::to_string(z_mean));
  check::close(z_var, 1.0, 0.2,
               "deviations over sigma, variance: " + std::to_string(z_var));
}

TEST(spots_that_cannot_be_told_apart_are_degenerate) {
  std::mt19937 random(5);
  const Spot s{20, 15, 6, 3000};
  const Volume v = sample({s}, random);
  Shoebox a, b;
  const OverlapMember ma = box_of(s, v, &a), mb = box_of(s, v, &b);
  const std::vector<ProfileFit> together = fit_together({ma, mb});
  check::is_true(!together[0].valid && !together[1].valid, "not fitted");
  check::is_true(together[0].why == FitFailure::degenerate, "degenerate");
}

} // namespace mxi
