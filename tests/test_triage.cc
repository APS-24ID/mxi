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
