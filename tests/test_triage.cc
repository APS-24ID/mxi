// mxi triage's statistics (src/triage.hh) on intensities drawn as Wilson's
// statistics have them, where each has a known value: untwinned, twinned at
// 0.3 and perfectly, and a symmetry missed (P 4 2 2 data taken as P 4).

#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <set>

#include "../src/triage.hh"
#include "check.hh"

namespace mxi {

namespace {

const UnitCell kCell{60, 60, 80, 90, 90, 90};

// Every unique reflection of `group` to 2 A, an intensity each: acentric
// exponential, centric eps Sigma x^2, Sigma falling off with B = 20; sigma
// 2 per cent of the shell's mean, noise added.
std::map<Miller, double> wilson(const SpaceGroup &group, std::uint32_t seed) {
  std::mt19937 random(seed);
  std::exponential_distribution<double> exponential(1.0);
  std::normal_distribution<double> normal(0.0, 1.0);
  std::map<Miller, double> out;
  for (std::int32_t h = -31; h <= 31; ++h)
    for (std::int32_t k = -31; k <= 31; ++k)
      for (std::int32_t l = -41; l <= 41; ++l) {
        const Miller m{h, k, l};
        if ((h == 0 && k == 0 && l == 0) || group.absent(m))
          continue;
        const double d = d_spacing(kCell, m);
        if (d < 2.0)
          continue;
        const Miller u = group.unique(m);
        if (out.count(u))
          continue;
        const double sigma =
            1000.0 * std::exp(-20.0 / (2.0 * d * d)) * epsilon(group, u);
        const double x = normal(random);
        out[u] = group.centric(u) ? sigma * x * x : sigma * exponential(random);
      }
  return out;
}

TriageData as_data(const std::map<Miller, double> &truth,
                   const SpaceGroup &group, std::uint32_t seed) {
  std::mt19937 random(seed);
  std::normal_distribution<double> normal(0.0, 1.0);
  TriageData data;
  for (const auto &[m, i] : truth) {
    const double d = d_spacing(kCell, m);
    const double s =
        0.02 * 1000.0 * std::exp(-20.0 / (2.0 * d * d)) * epsilon(group, m);
    data.hkl.push_back(m);
    data.intensity.push_back(i + s * normal(random));
    data.sigma.push_back(s);
  }
  return data;
}

// I(h) = (1 - a) I(h) + a I(T h): a twin of fraction a by the twin law T.
std::map<Miller, double> twinned(const std::map<Miller, double> &truth,
                                 const SpaceGroup &group, const Rotation &law,
                                 double a) {
  std::map<Miller, double> out;
  for (const auto &[m, i] : truth) {
    const auto mate = truth.find(group.unique(mxi::apply(law, m)));
    out[m] = mate == truth.end() ? i : (1 - a) * i + a * mate->second;
  }
  return out;
}

} // namespace

// Intensities with an anisotropic fall-off along the cell's axes,
// exp(-(B_a s_a^2 + B_b s_b^2 + B_c s_c^2) / 2) with s = (h/a, k/b, l/c), an
// orthogonal cell; Wilson's statistics; sigma a constant, so the data reach
// further where they fall off least. The fit zeroes negative intensities, so
// a noise floor flattens the fall-off where it is steepest: at sigma 3 the
// spread found is some 15 per cent short, at 30 half -- so the tensor is
// checked on precise data (0.3), the reach on noisier (3), where it ends
// inside the data. Every unique reflection to 1.5 A.
TriageData anisotropic(const UnitCell &cell, const SpaceGroup &group, double ba,
                       double bb, double bc, double noise, std::uint32_t seed) {
  std::mt19937 random(seed);
  std::exponential_distribution<double> exponential(1.0);
  std::normal_distribution<double> normal(0.0, 1.0);
  std::set<Miller> seen;
  TriageData data;
  const std::int32_t hm = static_cast<std::int32_t>(cell.a / 1.5),
                     km = static_cast<std::int32_t>(cell.b / 1.5),
                     lm = static_cast<std::int32_t>(cell.c / 1.5);
  for (std::int32_t h = -hm; h <= hm; ++h)
    for (std::int32_t k = -km; k <= km; ++k)
      for (std::int32_t l = -lm; l <= lm; ++l) {
        const Miller m{h, k, l};
        if ((h == 0 && k == 0 && l == 0) || group.absent(m) ||
            d_spacing(cell, m) < 1.5)
          continue;
        const Miller unique = group.unique(m);
        if (!seen.insert(unique).second)
          continue;
        const double sa = h / cell.a, sb = k / cell.b, sc = l / cell.c;
        const double mean =
            1000.0 * epsilon(group, unique) *
            std::exp(-0.5 * (ba * sa * sa + bb * sb * sb + bc * sc * sc));
        const double x = normal(random);
        const double i =
            group.centric(unique) ? mean * x * x : mean * exponential(random);
        data.hkl.push_back(unique);
        data.intensity.push_back(i + noise * normal(random));
        data.sigma.push_back(noise);
      }
  return data;
}

TEST(triage_anisotropy_along_three_axes_is_found_and_reached) {
  // Orthorhombic, B 10, 20 and 40 A^2 along a, b and c: three principal
  // directions, the weakest along c* and the strongest along a*; their
  // differences as made (the fit's own atoms' fall-off shifts each alike);
  // the data reaching further along a* than along c*.
  const UnitCell cell{60, 70, 80, 90, 90, 90};
  const SpaceGroup group = SpaceGroup::from_name("P 2 2 2");
  const Triage t =
      triage(anisotropic(cell, group, 10, 20, 40, 0.3, 7), cell, group, {});
  const Anisotropy &a = t.anisotropy;
  check::is_true(a.fitted, "fitted");
  check::equal(a.distinct, 3, "three directions");
  check::close(a.b[0] - a.b[2], 30.0, 3.0, "the spread");
  check::close(a.b[1] - a.b[2], 10.0, 2.0, "b's against a's");
  check::is_true(a.angle[0][2] < 5.0, "the weakest along c*");
  check::is_true(a.angle[2][0] < 5.0, "the strongest along a*");
  const Triage noisy =
      triage(anisotropic(cell, group, 10, 20, 40, 3.0, 7), cell, group, {});
  const Anisotropy &n = noisy.anisotropy;
  check::is_true(n.d_min[2] > 0.0 && n.d_min[0] > n.d_min[2] + 0.3,
                 "reaching further along a* than c*");
}

TEST(triage_anisotropy_of_a_tetragonal_crystal_is_an_axis_and_a_plane) {
  // B 15, 15, 40: two principal values one, so a unique axis, c*, and the
  // plane across it -- their directions in the plane arbitrary, not given.
  const UnitCell cell{60, 60, 80, 90, 90, 90};
  const SpaceGroup group = SpaceGroup::from_name("P 4");
  const Triage t =
      triage(anisotropic(cell, group, 15, 15, 40, 0.3, 11), cell, group, {});
  const Anisotropy &a = t.anisotropy;
  check::equal(a.distinct, 2, "an axis and a plane");
  check::is_true(a.angle[a.unique][2] < 5.0, "the axis c*");
  check::close(a.b[a.unique] - a.b_plane, 25.0, 3.0, "the axis's B the more");
  const Triage noisy =
      triage(anisotropic(cell, group, 15, 15, 40, 3.0, 11), cell, group, {});
  const Anisotropy &n = noisy.anisotropy;
  check::equal(n.distinct, 2, "an axis and a plane, noisier");
  check::is_true(n.d_min_plane > 0.0 && n.d_min[n.unique] > n.d_min_plane,
                 "reaching further across than along");
}

TEST(triage_anisotropy_of_a_cubic_crystal_is_none) {
  // Fitted to an asymmetric unit, the tensor came out 28, 21 and 19 A^2 for
  // cubic insulin until it was averaged over the point group.
  const UnitCell cell{60, 60, 60, 90, 90, 90};
  const SpaceGroup group = SpaceGroup::from_name("P 2 3");
  const Triage t =
      triage(anisotropic(cell, group, 20, 20, 20, 3.0, 13), cell, group, {});
  check::equal(t.anisotropy.distinct, 1, "isotropic");
  check::is_true(t.anisotropy.b[0] - t.anisotropy.b[2] < 0.01,
                 "the same B every way");
}

TEST(triage_twin_law_text_and_epsilon) {
  check::is_true(twin_law_text({0, 1, 0, 1, 0, 0, 0, 0, -1}) == "k,h,-l",
                 "k,h,-l");
  check::is_true(twin_law_text({-1, 0, 0, 0, 1, 0, 0, 0, -1}) == "-h,k,-l",
                 "-h,k,-l");
  const SpaceGroup p4 = SpaceGroup::from_name("P 4");
  check::equal(epsilon(p4, {0, 0, 4}), 4, "00l on the four-fold");
  check::equal(epsilon(p4, {3, 2, 0}), 1, "a general reflection");
  const std::vector<TwinLaw> laws = twin_laws(kCell, p4, 3.0);
  check::equal(static_cast<std::int64_t>(laws.size()), std::int64_t{1},
               "P 4 in a P 4 2 2 lattice: one law");
  check::is_true(laws[0].merohedral, "merohedral");
  check::equal(
      static_cast<std::int64_t>(
          twin_laws(kCell, SpaceGroup::from_name("P 4 2 2"), 3.0).size()),
      std::int64_t{0}, "the lattice's own group: none");
}

TEST(triage_untwinned_data_are_as_wilson_says) {
  const SpaceGroup p4 = SpaceGroup::from_name("P 4");
  const Triage t = triage(as_data(wilson(p4, 7), p4, 8), kCell, p4);
  check::close(t.d_min, 2.0, 0.05, "strong to the edge: the data's limit");
  check::close(t.acentric_moments.i2_over_i_squared, 2.0, 0.1, "<I^2>/<I>^2");
  check::close(t.acentric_moments.f_squared_over_f2, 0.785, 0.02,
               "<F>^2/<F^2>");
  check::close(t.acentric_moments.e2_minus_one, 0.736, 0.03, "<|E^2-1|>");
  check::close(t.centric_moments.i2_over_i_squared, 3.0, 0.4,
               "centric <I^2>/<I>^2");
  check::close(t.nz.max_deviation_acentric, 0.0, 0.02, "N(z) as Wilson's");
  check::close(t.l_test.mean_abs_l, 0.5, 0.02, "<|L|>");
  check::close(t.l_test.mean_l2, 1.0 / 3.0, 0.02, "<L^2>");
  check::is_true(!t.abnormal && !t.wide && t.l_test.z_score < 3.5,
                 "the L test normal, Z " + std::to_string(t.l_test.z_score));
  check::close(t.twin_laws[0].h_alpha, 0.0, 0.05, "no twin");
  check::close(t.twin_laws[0].r_obs, 0.5, 0.05, "unrelated pairs' R");
  check::is_true(t.verdict.find("no twinning is suspected") !=
                     std::string::npos,
                 t.verdict);
}

TEST(triage_with_no_acentric_reflections_says_it_cannot_judge) {
  // A centrosymmetric group: every reflection centric, the L test with no
  // acentric neighbours -- which had said "as expected: no twinning is
  // suspected" from statistics of nothing (a small molecule in R -3 m,
  // 13 November).
  const SpaceGroup p4m = SpaceGroup::from_name("P 4/m");
  const Triage t = triage(as_data(wilson(p4m, 7), p4m, 8), kCell, p4m);
  check::equal(static_cast<std::int64_t>(t.l_test.pairs), std::int64_t{0},
               "no acentric pairs");
  check::is_true(t.verdict.find("cannot judge twinning") != std::string::npos,
                 t.verdict);
  check::is_true(t.verdict.find("no twinning is suspected") ==
                     std::string::npos,
                 "and no verdict from nothing");
}

TEST(triage_twinned_data_are_said_so) {
  const SpaceGroup p4 = SpaceGroup::from_name("P 4");
  const Rotation law = twin_laws(kCell, p4, 3.0)[0].op;
  const auto truth = wilson(p4, 7);
  const Triage t =
      triage(as_data(twinned(truth, p4, law, 0.3), p4, 8), kCell, p4);
  check::close(t.twin_laws[0].h_alpha, 0.3, 0.04, "the H test's fraction");
  check::is_true(t.abnormal, "the L test abnormal, as twinning makes it");
  check::is_true(t.l_test.mean_abs_l < 0.47, "<|L|> down from 1/2");
  check::is_true(t.verdict.find("twinning could be why") != std::string::npos,
                 t.verdict);
  const Triage perfect =
      triage(as_data(twinned(truth, p4, law, 0.5), p4, 8), kCell, p4);
  // 1.5 with eps divided out of E^2, as it is defined; xtriage leaves it in,
  // and the 00l reflections (acentric in P 4, eps 4) had made it 1.58.
  check::close(perfect.acentric_moments.i2_over_i_squared, 1.5, 0.08,
               "perfect: 1.5");
  check::close(perfect.l_test.mean_abs_l, 0.375, 0.02, "perfect: <|L|> 3/8");
  check::close(perfect.twin_laws[0].h_alpha, 0.5, 0.03, "perfect: 1/2");
}

TEST(triage_a_symmetry_missed_is_not_a_twin) {
  // P 4 2 2 data taken as P 4: the pairs k,h,-l relates equal, the
  // intensities' statistics untwinned.
  const SpaceGroup p422 = SpaceGroup::from_name("P 4 2 2"),
                   p4 = SpaceGroup::from_name("P 4");
  const auto truth = wilson(p422, 7);
  std::map<Miller, double> low;
  for (const auto &[m, i] :
       wilson(p4, 9)) // the P 4 indices, each its P 4 2 2 value
    low[m] = truth.at(p422.unique(m));
  const Triage t = triage(as_data(low, p4, 8), kCell, p4);
  // Not twin-like: the zones P 4 2 2 makes centric, taken as acentric, widen
  // the spread -- <|L|> above 1/2, as xtriage's own example of a symmetry
  // missed (0.512).
  check::is_true(!t.abnormal, "not abnormal as twinning makes it, Z " +
                                  std::to_string(t.l_test.z_score) + " <|L|> " +
                                  std::to_string(t.l_test.mean_abs_l));
  check::is_true(t.twin_laws[0].h_alpha > 0.4 && t.twin_laws[0].r_obs < 0.1,
                 "pairs equal");
  check::is_true(t.verdict.find("symmetry element missed") != std::string::npos,
                 t.verdict);
}

} // namespace mxi
