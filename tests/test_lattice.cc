#include <cmath>

#include "../src/symmetry.hh"
#include "check.hh"

namespace mxi {

namespace {

//: A primitive cell from conventional axes and the vectors of its primitive
//: cell in terms of them.
UnitCell primitive_of(const UnitCell &conventional,
                      const char *primitive_axes) {
  return reindexed_cell(conventional, ChangeOfBasis::parse(primitive_axes));
}

std::optional<Setting> named(const UnitCell &cell, std::size_t order) {
  for (const auto &g : subgroups(lattice_symmetry(cell)))
    if (g.size() == order)
      return reference_setting(g, cell);
  return std::nullopt;
}

} // namespace

TEST(the_insulin_cells_lattice_is_body_centred_cubic) {
  // This pipeline's refined cell of the 300 image insulin sweep: 24 rotations,
  // and their 30 subgroups -- as many as S4, the rotation group of a cube, has.
  const UnitCell cell{67.4201, 67.4606, 67.4644, 109.443, 109.442, 109.459};
  const std::vector<Rotation> rots = lattice_symmetry(cell);
  check::equal(static_cast<long long>(rots.size()), 24, "24 rotations");
  check::equal(static_cast<long long>(subgroups(rots).size()), 30,
               "30 subgroups");
  const auto full = reference_setting(rots, cell);
  check::is_true(full && full->group.name() == "I m -3 m",
                 "the lattice is I m -3 m");
  const UnitCell cubic = reindexed_cell(cell, full->cb);
  check::close(cubic.a, 77.9, 0.1, "at 77.9 A");
  check::close(cubic.beta, 90.0, 0.1, "and 90 degrees");
  const auto m3 = named(cell, 12);
  check::is_true(m3 && m3->group.name() == "I m -3",
                 "its subgroup of 12 rotations is I m -3");
}

TEST(a_c_centred_monoclinic_lattice_is_found_from_its_primitive_cell) {
  const UnitCell conventional{90.0, 40.0, 60.0, 90.0, 105.0, 90.0};
  const UnitCell p = primitive_of(conventional, "1/2a-1/2b,1/2a+1/2b,c");
  const std::vector<Rotation> rots = lattice_symmetry(p);
  check::equal(static_cast<long long>(rots.size()), 2,
               "one two-fold and the identity");
  const auto s = reference_setting(rots, p);
  check::is_true(s && s->group.name() == "C 1 2/m 1", "C 1 2/m 1");
  const UnitCell back = reindexed_cell(p, s->cb);
  check::close(back.b, 40.0, 1e-6, "the unique axis is b, 40 A");
  check::close(
      back.a * back.c * std::sin(back.beta * std::acos(-1.0) / 180.0) * back.b,
      90.0 * 40.0 * 60.0 * std::sin(105.0 * std::acos(-1.0) / 180.0), 1e-3,
      "and the cell has the conventional volume, twice the primitive");
}

TEST(a_single_axis_finds_its_perpendicular_axes) {
  // A four-fold alone has no rotation fixing a or b; a three-fold does not
  // even reverse them. Tetragonal P, and rhombohedral from its primitive cell.
  const UnitCell tetragonal{50.0, 50.0, 80.0, 90.0, 90.0, 90.0};
  const auto t = named(tetragonal, 4);
  check::is_true(
      t && (t->group.name() == "P 4/m" || t->group.name() == "P m m m"),
      "an order-four subgroup of P 4/m m m is named");
  bool four = false;
  for (const auto &g : subgroups(lattice_symmetry(tetragonal)))
    if (g.size() == 4)
      if (auto s = reference_setting(g, tetragonal);
          s && s->group.name() == "P 4/m")
        four = true;
  check::is_true(four, "and P 4/m among them");
  const UnitCell rhombohedral{60.0, 60.0, 60.0, 80.0, 80.0, 80.0};
  const auto r = named(rhombohedral, 3);
  check::is_true(r && r->group.name() == "R -3:H",
                 "a three-fold alone is R -3 in hexagonal axes");
  const UnitCell hex = reindexed_cell(rhombohedral, r->cb);
  check::close(hex.gamma, 120.0, 1e-6, "gamma 120");
}

TEST(a_triclinic_lattice_is_itself) {
  const UnitCell cell{51.0, 62.0, 73.0, 81.0, 76.0, 69.0};
  const std::vector<Rotation> rots = lattice_symmetry(cell);
  check::equal(static_cast<long long>(rots.size()), 1, "the identity alone");
  const auto s = reference_setting(rots, cell);
  check::is_true(s && s->group.name() == "P -1" && s->cb_text == "a,b,c",
                 "P -1, as it is");
}

} // namespace mxi

namespace mxi {

TEST(a_primitive_monoclinic_cell_has_beta_obtuse_and_a_shorter_than_c) {
  // Graeme's 6S0R: dials.symmetry's P 1 2/m 1, 50.287 183.849 109.941 A,
  // beta 94.234 -- mxi had 110.075 183.899 50.214, a and c the wrong way
  // round, the convention's a <= c lost to a tie-break on minus signs. Here as
  // the reduced cell it came from: the two-fold along c, a^b 94.234.
  const UnitCell reduced{50.287, 109.941, 183.849, 90.0, 90.0, 94.234};
  const std::vector<Rotation> rots = lattice_symmetry(reduced);
  check::equal(static_cast<long long>(rots.size()), 2,
               "one two-fold and the identity");
  const auto s = reference_setting(rots, reduced);
  check::is_true(s && s->group.name() == "P 1 2/m 1", "P 1 2/m 1");
  const UnitCell c = reindexed_cell(reduced, s->cb);
  check::close(c.a, 50.287, 1e-3, "a the shorter, 50.287");
  check::close(c.b, 183.849, 1e-3, "b the unique axis");
  check::close(c.c, 109.941, 1e-3, "c the longer, 109.941");
  check::close(c.beta, 94.234, 1e-3, "beta obtuse");
}

} // namespace mxi
