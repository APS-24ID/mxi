// mxi merge's library: merging, French and Wilson, FreeR flags (26 November).
#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

#include "../src/french_wilson_tables.hh"
#include "../src/merge.hh"
#include "check.hh"

namespace mxi {

TEST(equivalents_merge_and_friedel_mates_stay_apart_when_asked) {
  const gemmi::SpaceGroup *p4 = gemmi::find_spacegroup_by_name("P 4");
  // (1 2 3) and its four-fold mate (-2 1 3) are one reflection; (-1 -2 -3) its
  // Friedel mate; (0 0 4) centric is not, P 4 being acentric... (h k 0) is
  // centric in P 4: (3 1 0).
  const std::vector<Miller> hkl = {
      {1, 2, 3}, {-2, 1, 3}, {-1, -2, -3}, {3, 1, 0}, {-3, -1, 0}};
  const std::vector<double> I = {100, 120, 80, 50, 70}, v = {4, 4, 4, 1, 1};
  const MergedEquivalents mean = merge_equivalents(hkl, I, v, *p4, false);
  check::equal(static_cast<std::int64_t>(mean.size()), 2, "two reflections");
  // (1 2 3) with its mates, all three: (100 + 120 + 80) / 3.
  std::int32_t found = 0;
  for (std::size_t i = 0; i < mean.size(); ++i) {
    if (mean.n[i] == 3) {
      check::close(mean.value[i], 100.0, 1e-9, "the mean of three");
      check::close(mean.sigma[i], 2.0 / std::sqrt(3.0), 1e-9, "1/sqrt(sum w)");
      ++found;
    }
    if (mean.n[i] == 2) {
      check::is_true(mean.centric[i], "(3 1 0) centric");
      ++found;
    }
  }
  check::equal(found, 2, "both found");
  const MergedEquivalents anom = merge_equivalents(hkl, I, v, *p4, true);
  check::equal(static_cast<std::int64_t>(anom.size()), 3,
               "I(+), I(-), and one centric");
  for (std::size_t i = 0; i < anom.size(); ++i) {
    if (anom.centric[i])
      check::equal(static_cast<std::int64_t>(anom.hemisphere[i]), 0,
                   "a centric has one entry");
    else if (anom.hemisphere[i] > 0)
      check::close(anom.value[i], 110.0, 1e-9, "I(+): (1 2 3), (-2 1 3)");
    else
      check::close(anom.value[i], 80.0, 1e-9, "I(-): (-1 -2 -3)");
  }
}

TEST(french_and_wilson_as_dials_computes_them) {
  namespace t = french_wilson_tables;
  // Many reflections of one shell, the expected intensity about 1000.
  std::vector<double> I, sigma, d3;
  std::vector<bool> centric;
  for (std::int32_t i = 0; i < 400; ++i) {
    I.push_back(1000.0 * std::log(400.0 / (i + 0.5)) * std::log(2.0) /
                std::log(2.0)); // exponential-ish
    sigma.push_back(20.0);
    d3.push_back(0.001 + 1e-6 * i);
    centric.push_back(false);
  }
  // A strong one, a negative one, one on a table point (h = 0.0 exactly is
  // unreachable with I free, so test the table's interpolation on its own).
  I[0] = 50000.0;
  I[1] = -40.0;
  const FrenchWilson fw = french_wilson(I, sigma, centric, d3);
  check::is_true(fw.valid[0] && fw.valid[1], "both kept");
  // Strong: F = sqrt(h sigma), h = I/sigma - sigma/Sigma.
  check::is_true(std::fabs(fw.F[0] - std::sqrt(50000.0)) < 0.05,
                 "strong: F near sqrt(I)");
  check::close(fw.sigF[0], 0.5 * 20.0 / fw.F[0], 1e-12, "sigF = sigma / 2F");
  check::is_true(fw.F[1] > 0.0 && fw.F[1] < std::sqrt(3 * 20.0),
                 "a negative intensity, a small positive F");
  // Refused: I/sigma below -3.7.
  std::vector<double> I2 = I;
  I2[2] = -100.0; // I/sigma -5
  check::is_true(!french_wilson(I2, sigma, centric, d3).valid[2],
                 "I/sigma -5 refused");
  check::close(t::ac_zf[40], 0.822, 1e-12, "the table as cctbx's, at h = 0");
}

TEST(free_flags_are_shared_by_lattice_mates_and_friedel_mates) {
  const UnitCell cubic{78.0, 78.0, 78.0, 90.0, 90.0, 90.0};
  const std::vector<Rotation> lattice = lattice_symmetry(cubic, 5.0);
  check::equal(static_cast<std::int64_t>(lattice.size()), 24,
               "the cube's 24 rotations");
  const Miller h{3, 5, 7};
  const std::int32_t f = free_flag(h, lattice, 20);
  for (const Rotation &r : lattice) {
    const Miller m = apply(r, h);
    check::equal(static_cast<std::int64_t>(free_flag(m, lattice, 20)),
                 static_cast<std::int64_t>(f), "a lattice mate, one flag");
    check::equal(static_cast<std::int64_t>(
                     free_flag(Miller{-m[0], -m[1], -m[2]}, lattice, 20)),
                 static_cast<std::int64_t>(f), "a Friedel mate, one flag");
  }
  // Over many distinct reflections, a twentieth in the test set.
  std::set<Miller> seen;
  std::int32_t test = 0, total = 0;
  for (std::int32_t a = 0; a <= 30; ++a)
    for (std::int32_t b = 0; b <= a; ++b)
      for (std::int32_t c = 0; c <= b; ++c) {
        if (a == 0)
          continue;
        ++total;
        test += free_flag(Miller{a, b, c}, lattice, 20) == 0;
      }
  const double fraction = static_cast<double>(test) / total;
  check::is_true(fraction > 0.04 && fraction < 0.06,
                 "a twentieth: " + std::to_string(fraction));
  check::equal(static_cast<std::int64_t>(free_flag(h, lattice, 20)),
               static_cast<std::int64_t>(f), "the same again");
}

} // namespace mxi
