#include <cmath>
#include <set>

#include "../src/symmetry.hh"
#include "check.hh"

namespace mxi {

TEST(every_equivalent_of_a_reflection_in_m3_is_one_reflection) {
  // I 2 3, Laue class m-3: 24 operations counting Friedel mates -- the three
  // cyclic permutations of (h, k, l), times sign changes in pairs, times
  // inversion. All 24 images of (1, 2, 3) are one reflection.
  const SpaceGroup g = SpaceGroup::from_name("I 2 3");
  check::is_true(g.laue() == "m-3", "the Laue class is m-3: " + g.laue());
  check::equal(static_cast<long long>(g.order()), 12,
               "12 operations without centring");
  const Miller h{1, 2, 3};
  const Miller u = g.unique(h);
  std::set<Miller> images;
  const int signs[4][3] = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
  for (int c = 0; c < 3; ++c) {
    const Miller p{h[c % 3], h[(c + 1) % 3], h[(c + 2) % 3]};
    for (const auto &s : signs) {
      for (int inv : {1, -1}) {
        const Miller e{inv * s[0] * p[0], inv * s[1] * p[1], inv * s[2] * p[2]};
        images.insert(e);
        check::is_true(g.unique(e) == u,
                       "every image of (1,2,3) is one reflection");
      }
    }
  }
  check::equal(static_cast<long long>(images.size()), 24,
               "and there are 24 of them");
  // A non-cyclic permutation is m-3m's symmetry, not m-3's: a different
  // reflection. This is what would catch the wrong Laue class.
  check::is_true(g.unique({1, 3, 2}) != u,
                 "(1,3,2) is another reflection in m-3");
  const SpaceGroup m3m = SpaceGroup::from_name("I 4 3 2");
  check::is_true(m3m.unique({1, 3, 2}) == m3m.unique(h),
                 "but the same one in m-3m");
}

TEST(centring_forbids_reflections_and_p1_pairs_only_friedel_mates) {
  const SpaceGroup g = SpaceGroup::from_name("I 2 3");
  check::is_true(g.absent({1, 0, 0}),
                 "h + k + l odd is absent under I centring");
  check::is_true(!g.absent({1, 1, 0}), "and even is not");
  const SpaceGroup p1 = SpaceGroup::from_hall(" P 1");
  check::equal(static_cast<long long>(p1.number()), 1,
               "an .expt's ' P 1' is P 1");
  check::is_true(p1.unique({1, 2, 3}) == p1.unique({-1, -2, -3}),
                 "Friedel mates merge");
  check::is_true(p1.unique({1, 2, 3}) != p1.unique({1, 2, -3}),
                 "and nothing else does");
  check::is_true(SpaceGroup::from_hall(g.hall()).name() == g.name(),
                 "a Hall symbol written out reads back as the same group");
}

TEST(a_change_of_basis_takes_the_primitive_cell_to_the_conventional) {
  // DIALS' b+c,a+c,a+b: the primitive cell of a body-centred cubic lattice to
  // the conventional one. hkl' = M^T hkl, so (1,2,3) -> (2+3, 1+3, 1+2).
  const ChangeOfBasis cb = ChangeOfBasis::parse("b+c,a+c,a+b");
  check::is_true(cb.apply({1, 2, 3}) == Miller{5, 4, 3}, "(1,2,3) -> (5,4,3)");
  const Miller n = cb.apply({1, 0, 0});
  check::equal(static_cast<long long>((n[0] + n[1] + n[2]) % 2), 0,
               "every reindexed index obeys I centring");
  // Its inverse, with fractions, undoes it; and a fractional index is refused.
  const ChangeOfBasis back =
      ChangeOfBasis::parse("-1/2a+1/2b+1/2c,1/2a-1/2b+1/2c,1/2a+1/2b-1/2c");
  check::is_true(back.apply({5, 4, 3}) == Miller{1, 2, 3},
                 "the inverse undoes it");
  bool refused = false;
  try {
    back.apply({1, 0, 0});
  } catch (const std::exception &) {
    refused = true;
  }
  check::is_true(refused,
                 "an index with no integral image is refused, not rounded");

  // A crystal: a cubic cell of 78 A, taken to its primitive cell and back.
  Crystal c;
  c.A = Mat3::identity() * (1.0 / 78.0);
  back.apply(c);
  const UnitCell primitive = c.cell();
  check::close(primitive.a, 78.0 * std::sqrt(3.0) / 2.0, 1e-9,
               "primitive a = a sqrt(3)/2");
  check::close(primitive.alpha, std::acos(-1.0 / 3.0) * 180.0 / std::acos(-1.0),
               1e-9, "and 109.47 degrees");
  cb.apply(c);
  const UnitCell cubic = c.cell();
  check::close(cubic.a, 78.0, 1e-9, "and back to 78 A");
  check::close(cubic.beta, 90.0, 1e-9, "cubic");
}

TEST(a_change_of_basis_is_refused_when_it_cannot_be_one) {
  for (const char *bad : {"a,b", "a+b,a+b,c", "a,b,d", "a b,b,c"}) {
    bool refused = false;
    try {
      ChangeOfBasis::parse(bad);
    } catch (const std::exception &) {
      refused = true;
    }
    check::is_true(refused, std::string("refused: ") + bad);
  }
}

} // namespace mxi

namespace mxi {

namespace {

UnitCell cell_of(const Mat3 &A) {
  Crystal c;
  c.A = A;
  return c.cell();
}

// U of A = U B, B upper triangular: what regularising must keep.
Mat3 orientation(const Mat3 &A) {
  const Mat3 g = A.transpose() * A;
  const auto at = [&](int i, int j) {
    return g.m[static_cast<std::size_t>(i * 3 + j)];
  };
  const double r00 = std::sqrt(at(0, 0)), r01 = at(0, 1) / r00,
               r02 = at(0, 2) / r00;
  const double r11 = std::sqrt(at(1, 1) - r01 * r01),
               r12 = (at(1, 2) - r01 * r02) / r11;
  const double r22 = std::sqrt(at(2, 2) - r02 * r02 - r12 * r12);
  const Mat3 B{r00, r01, r02, 0.0, r11, r12, 0.0, 0.0, r22};
  return A * B.inverse();
}

} // namespace

TEST(a_cubic_cell_is_made_cubic_and_its_orientation_kept) {
  // A cell near cubic, as refinement leaves one: after the group I 2 3, a = b
  // = c and every angle 90 -- at every scan point too -- and U as it was.
  const Mat3 real{77.90, 0.02, -0.01, 0.03, 77.95, 0.02, -0.02, 0.01, 77.86};
  const Mat3 U = rotation(Vec3{0.3, -0.5, 0.8}.normalized(), 0.7);
  Crystal crystal;
  crystal.A = U * real.inverse();
  crystal.A_points = {crystal.A, crystal.A};
  const Mat3 before = orientation(crystal.A);
  regularise_cell(crystal, SpaceGroup::from_name("I 2 3"));
  for (const Mat3 &A : {crystal.A, crystal.A_points[0], crystal.A_points[1]}) {
    const UnitCell c = cell_of(A);
    check::close(c.b, c.a, 1e-9, "b = a");
    check::close(c.c, c.a, 1e-9, "c = a");
    check::close(c.alpha, 90.0, 1e-9, "alpha 90");
    check::close(c.beta, 90.0, 1e-9, "beta 90");
    check::close(c.gamma, 90.0, 1e-9, "gamma 90");
  }
  const Mat3 after = orientation(crystal.A);
  for (std::size_t k = 0; k < 9; ++k)
    check::close(after.m[k], before.m[k], 1e-9, "the orientation kept");
}

TEST(a_monoclinic_cell_loses_only_what_its_two_fold_forbids) {
  // P 1 2 1, the two-fold along b: alpha and gamma made 90, a, b, c and beta
  // exactly as they were. And P 1 changes nothing.
  const Mat3 real{50.0, 0.3, 0.0, 0.0, 60.0, 0.4, -10.0, 0.2, 70.0};
  Crystal crystal;
  crystal.A = real.inverse();
  const UnitCell before = crystal.cell();
  Crystal p1 = crystal;
  regularise_cell(crystal, SpaceGroup::from_name("P 1 2 1"));
  const UnitCell after = crystal.cell();
  check::close(after.alpha, 90.0, 1e-9, "alpha 90");
  check::close(after.gamma, 90.0, 1e-9, "gamma 90");
  check::close(after.beta, before.beta, 1e-9, "beta kept");
  check::close(after.a, before.a, 1e-9, "a kept");
  check::close(after.b, before.b, 1e-9, "b kept");
  check::close(after.c, before.c, 1e-9, "c kept");
  regularise_cell(p1, SpaceGroup::from_name("P 1"));
  for (std::size_t k = 0; k < 9; ++k)
    check::close(p1.A.m[k], real.inverse().m[k], 0.0, "P 1 as it was");
}

} // namespace mxi
