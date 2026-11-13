#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "../src/index.hh"
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
  check::equal(static_cast<std::int64_t>(rots.size()), 24, "24 rotations");
  check::equal(static_cast<std::int64_t>(subgroups(rots).size()), 30,
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
  check::equal(static_cast<std::int64_t>(rots.size()), 2,
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

TEST(a_rhombohedral_lattice_in_another_primitive_cell_is_named) {
  // Graeme's small molecule, 13 November: its refined primitive cell, not the
  // rhombohedral one, in which the three-fold -- the hexagonal c -- is
  // [2,-1,-3]. The axes were sought with components -2..2 only, so R -3 m and
  // R -3, which scored best, had no c to lie along, could not be named, and
  // mxi symmetry stopped. The rhombohedral cell of the test above has the
  // three-fold at [1,1,1], which is why it passed.
  const UnitCell cell{9.511, 9.508, 10.248, 89.985, 62.363, 59.989};
  check::equal(static_cast<std::int64_t>(lattice_symmetry(cell).size()), 6,
               "six rotations: a three-fold and three two-folds");
  const auto full = named(cell, 6);
  check::is_true(full && full->group.name() == "R -3 m:H", "R -3 m:H");
  const auto three = named(cell, 3);
  check::is_true(three && three->group.name() == "R -3:H",
                 "the three-fold alone R -3:H");
  const UnitCell hex = reindexed_cell(cell, full->cb);
  check::close(hex.a, 9.51, 0.01, "a 9.51 A");
  check::close(hex.b, hex.a, 0.01, "b as a");
  check::close(hex.c, 25.96, 0.01, "c 25.96 A, along the three-fold");
  check::close(hex.gamma, 120.0, 0.05, "gamma 120");
  check::close(hex.alpha, 90.0, 0.05, "alpha 90");
}

namespace {

//: The real-space axes of a cell, as rows: a along x, b in the xy plane.
Mat3 axes_of(const UnitCell &c) {
  const double d = std::acos(-1.0) / 180.0;
  const double ca = std::cos(c.alpha * d), cb = std::cos(c.beta * d),
               cg = std::cos(c.gamma * d), sg = std::sin(c.gamma * d);
  const double cx = c.c * cb, cy = c.c * (ca - cb * cg) / sg;
  Mat3 m;
  const double v[9] = {c.a,
                       0.0,
                       0.0,
                       c.b * cg,
                       c.b * sg,
                       0.0,
                       cx,
                       cy,
                       std::sqrt(std::max(0.0, c.c * c.c - cx * cx - cy * cy))};
  for (std::size_t i = 0; i < 9; ++i)
    m.m[i] = v[i];
  return m;
}

//: The cell whose real-space axes are these rows.
UnitCell cell_of_rows(const Mat3 &r) {
  const auto row = [&](std::size_t i) {
    return Vec3{r.m[i * 3], r.m[i * 3 + 1], r.m[i * 3 + 2]};
  };
  const auto length = [](const Vec3 &x) {
    return std::sqrt(x.x * x.x + x.y * x.y + x.z * x.z);
  };
  const auto angle = [&](const Vec3 &x, const Vec3 &y) {
    return std::acos((x.x * y.x + x.y * y.y + x.z * y.z) /
                     (length(x) * length(y))) *
           180.0 / std::acos(-1.0);
  };
  return UnitCell{length(row(0)),        length(row(1)),
                  length(row(2)),        angle(row(1), row(2)),
                  angle(row(0), row(2)), angle(row(0), row(1))};
}

} // namespace

TEST(every_lattice_is_named_from_any_reduced_primitive_cell) {
  // Each Bravais lattice -- the rhombohedral in five shapes, flat to long --
  // its primitive cell in 25 random bases, each reduced as indexing reduces
  // a basis: the full lattice group named every time. Before the axes were
  // sought to components of 3, every rhombohedral lattice longer than about
  // c/a 1 failed in every basis, its group unnamed and mxi symmetry stopped
  // (13 November); the flat ones, and every other lattice, passed.
  const double h = 0.5, t = 1.0 / 3.0;
  using Fractions = std::vector<std::array<double, 3>>;
  const Fractions P = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
                  I = {{-h, h, h}, {h, -h, h}, {h, h, -h}},
                  F = {{0, h, h}, {h, 0, h}, {h, h, 0}},
                  C = {{h, h, 0}, {-h, h, 0}, {0, 0, 1}},
                  R = {{2 * t, t, t}, {-t, t, t}, {-t, -2 * t, t}};
  struct Lattice {
    const char *name;
    UnitCell conventional;
    Fractions primitive;
    const char *group;
  };
  const std::vector<Lattice> lattices = {
      {"cP", {41, 41, 41, 90, 90, 90}, P, "P m -3 m"},
      {"cI", {57, 57, 57, 90, 90, 90}, I, "I m -3 m"},
      {"cF", {83, 83, 83, 90, 90, 90}, F, "F m -3 m"},
      {"tP", {43, 43, 71, 90, 90, 90}, P, "P 4/m m m"},
      {"tI", {47, 47, 113, 90, 90, 90}, I, "I 4/m m m"},
      {"oP", {37, 53, 71, 90, 90, 90}, P, "P m m m"},
      {"oC", {37, 59, 71, 90, 90, 90}, C, "C m m m"},
      {"oI", {37, 53, 79, 90, 90, 90}, I, "I m m m"},
      {"oF", {41, 61, 89, 90, 90, 90}, F, "F m m m"},
      {"hP", {47, 47, 67, 90, 90, 120}, P, "P 6/m m m"},
      {"hR, flat", {83, 83, 37, 90, 90, 120}, R, "R -3 m:H"},
      {"hR, middling", {113, 113, 97, 90, 90, 120}, R, "R -3 m:H"},
      {"hR, Graeme's", {9.509, 9.509, 25.964, 90, 90, 120}, R, "R -3 m:H"},
      {"hR, long", {61, 61, 200, 90, 90, 120}, R, "R -3 m:H"},
      {"hR, longer", {57, 57, 480, 90, 90, 120}, R, "R -3 m:H"},
      {"mP", {41, 53, 67, 90, 103, 90}, P, "P 1 2/m 1"},
      {"mC", {91, 41, 61, 90, 107, 90}, C, "C 1 2/m 1"},
      {"aP", {41, 53, 67, 79, 97, 107}, P, "P -1"}};
  std::mt19937 rng(7);
  std::uniform_int_distribution<std::int32_t> entry(-2, 2);
  for (const Lattice &lattice : lattices) {
    const Mat3 axes = axes_of(lattice.conventional);
    Mat3 primitive;
    for (std::size_t i = 0; i < 3; ++i)
      for (std::size_t j = 0; j < 3; ++j)
        primitive.m[i * 3 + j] = lattice.primitive[i][0] * axes.m[j] +
                                 lattice.primitive[i][1] * axes.m[3 + j] +
                                 lattice.primitive[i][2] * axes.m[6 + j];
    std::int32_t named_right = 0;
    for (std::int32_t k = 0; k < 25; ++k) {
      std::int32_t u[9];
      std::int32_t det = 0;
      do {
        for (std::int32_t &x : u)
          x = entry(rng);
        det = u[0] * (u[4] * u[8] - u[5] * u[7]) -
              u[1] * (u[3] * u[8] - u[5] * u[6]) +
              u[2] * (u[3] * u[7] - u[4] * u[6]);
      } while (det != 1);
      Mat3 rows;
      for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
          rows.m[i * 3 + j] = u[i * 3] * primitive.m[j] +
                              u[i * 3 + 1] * primitive.m[3 + j] +
                              u[i * 3 + 2] * primitive.m[6 + j];
      const UnitCell cell = cell_of_rows(reduce_basis(rows));
      const auto s = reference_setting(lattice_symmetry(cell), cell);
      if (s && s->group.name() == lattice.group)
        ++named_right;
    }
    check::equal(static_cast<std::int64_t>(named_right), 25,
                 std::string(lattice.name) + ": " + lattice.group +
                     " in every basis");
  }
}

TEST(a_triclinic_lattice_is_itself) {
  const UnitCell cell{51.0, 62.0, 73.0, 81.0, 76.0, 69.0};
  const std::vector<Rotation> rots = lattice_symmetry(cell);
  check::equal(static_cast<std::int64_t>(rots.size()), 1, "the identity alone");
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
  check::equal(static_cast<std::int64_t>(rots.size()), 2,
               "one two-fold and the identity");
  const auto s = reference_setting(rots, reduced);
  check::is_true(s && s->group.name() == "P 1 2/m 1", "P 1 2/m 1");
  const UnitCell c = reindexed_cell(reduced, s->cb);
  check::close(c.a, 50.287, 1e-3, "a the shorter, 50.287");
  check::close(c.b, 183.849, 1e-3, "b the unique axis");
  check::close(c.c, 109.941, 1e-3, "c the longer, 109.941");
  check::close(c.beta, 94.234, 1e-3, "beta obtuse");
}

// An orthorhombic cell's reference setting has its axes in order of length,
// a <= b <= c, as cctbx's best cell and DIALS report: every permutation had
// tied, and a small molecule's 5.42 8.13 12.02 A, already in order, was
// reindexed c,a,b (20 October). In order, it is left as it is.
TEST(an_orthorhombic_setting_has_its_axes_in_order_of_length) {
  for (const UnitCell &cell : {UnitCell{12.02, 5.42, 8.13, 90.0, 90.0, 90.0},
                               UnitCell{5.42, 8.13, 12.02, 90.0, 90.0, 90.0},
                               UnitCell{8.13, 12.02, 5.42, 90.0, 90.0, 90.0}}) {
    const auto s = named(cell, 4);
    check::is_true(s && s->group.name() == "P m m m", "P m m m");
    const UnitCell r = reindexed_cell(cell, s->cb);
    check::close(r.a, 5.42, 1e-6, "a the shortest, from " + s->cb_text);
    check::close(r.b, 8.13, 1e-6, "then b");
    check::close(r.c, 12.02, 1e-6, "c the longest");
  }
  const auto same = named(UnitCell{5.42, 8.13, 12.02, 90.0, 90.0, 90.0}, 4);
  check::is_true(same && same->cb_text == "a,b,c",
                 "in order already: a,b,c, not " + (same ? same->cb_text : ""));
}

} // namespace mxi
