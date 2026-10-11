// cosym on planted data sets (docs/multi-lattice.md): a cubic crystal whose
// symmetry is P23, in its P432 lattice -- a twofold indexing ambiguity, as
// I23's in I432 is.

#include <cmath>
#include <map>
#include <random>
#include <set>
#include <vector>

#include "../src/cosym.hh"
#include "check.hh"

namespace mxi {

namespace {

// P23's twelve rotations: the axes' cyclic permutations, each with an even
// number of their signs flipped.
std::vector<Rotation> p23() {
  std::vector<Rotation> out;
  const int perms[3][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}};
  const int signs[4][3] = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
  for (const auto &p : perms)
    for (const auto &s : signs) {
      Rotation r{};
      for (int i = 0; i < 3; ++i)
        r[3 * i + p[i]] = s[i];
      out.push_back(r);
    }
  return out;
}

// A data set: intensities constant on P23's orbits and Friedel pairs, each
// reflection present with probability one half, a little noise; its indices
// those of the true ones under `twist` -- the setting it was indexed in.
P1Intensities planted(const Rotation &twist, unsigned seed) {
  std::mt19937 rng(seed);
  std::mt19937 truth(7); // the crystal's own intensities, the same every time
  std::exponential_distribution<double> wilson(1.0);
  std::normal_distribution<double> noise(0.0, 0.03);
  std::uniform_real_distribution<double> coin(0.0, 1.0);
  const std::vector<Rotation> group = p23();
  std::map<Miller, double> value;
  const auto orbit_value = [&](const Miller &h) {
    Miller rep = h;
    for (const Rotation &g : group)
      for (int s : {1, -1}) {
        Miller x = apply(g, h);
        x = {s * x[0], s * x[1], s * x[2]};
        rep = std::min(rep, x);
      }
    auto it = value.find(rep);
    if (it == value.end())
      it = value.emplace(rep, wilson(truth)).first;
    return it->second;
  };
  P1Intensities d;
  for (int h = -8; h <= 8; ++h)
    for (int k = -8; k <= 8; ++k)
      for (int l = -8; l <= 8; ++l) {
        if (h == 0 && k == 0 && l == 0)
          continue;
        const double I = orbit_value({h, k, l}); // drawn for every reflection
        if (coin(rng) < 0.5)
          continue;
        d.hkl.push_back(apply(twist, {h, k, l})); // as this setting indexes it
        d.i.push_back(I + noise(rng));
        d.sigma.push_back(0.03);
        d.d.push_back(50.0 / std::sqrt(double(h * h + k * k + l * l)));
      }
  return d;
}

} // namespace

TEST(cosym_brings_planted_data_sets_to_one_setting) {
  const std::vector<Rotation> lattice =
      lattice_symmetry(UnitCell{50, 50, 50, 90, 90, 90});
  check::equal(static_cast<long long>(lattice.size()), 24,
               "P432's 24 rotations");
  const Rotation identity{1, 0, 0, 0, 1, 0, 0, 0, 1};
  const Rotation four{0, -1, 0, 1, 0, 0, 0, 0, 1}; // about c: not in P23
  const std::vector<P1Intensities> sets = {
      planted(identity, 1), planted(four, 2), planted(identity, 3)};
  const CosymResult r = cosym(sets, lattice);
  check::equal(static_cast<long long>(r.chosen.size()), 3, "one choice a set");
  check::equal(static_cast<long long>(r.chosen[0]), 0,
               "the reference as it is");
  check::is_true(r.cc[1] > 0.95,
                 "the twisted set, reindexed, one with the first: " +
                     std::to_string(r.cc[1]));
  check::equal(static_cast<long long>(r.chosen[2]), 0,
               "a set in the reference's setting left as it is");
  check::is_true(r.cc[2] > 0.95, "and one with the first");
  // Not by luck: unreindexed, the twisted set agrees far less -- only where
  // the 4-fold takes a reflection's orbit to itself, on the axes and some
  // diagonals, about half.
  const CosymResult alone = cosym({sets[0], sets[1]}, {identity});
  check::is_true(alone.cc[1] < r.cc[1] - 0.3,
                 "the twisted set as it was: " + std::to_string(alone.cc[1]));
}

} // namespace mxi
