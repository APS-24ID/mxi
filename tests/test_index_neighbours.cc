// inverse_nearest_neighbours (src/index.hh), the max cell's nearest-neighbour
// analysis: against every pair, on point sets of every awkward kind -- a
// lattice, uniform scatter, a plane (cells small, the grid sparse), clusters,
// coincident pairs (a reflection seen twice: not neighbours), a few points --
// and the same on any number of threads.

#include <cmath>
#include <random>
#include <vector>

#include "../src/index.hh"
#include "../src/parallel.hh"
#include "check.hh"

namespace mxi {

namespace {

// Every pair: the reciprocal of each point's distance to its nearest other,
// coincident ones (within 1e-4) skipped, in members' order.
std::vector<double> every_pair(const std::vector<Vec3> &p,
                               const std::vector<std::size_t> &m) {
  std::vector<double> out;
  for (std::size_t a : m) {
    double best = 1e300;
    for (std::size_t b : m) {
      if (a == b)
        continue;
      const double d = (p[b] - p[a]).norm();
      if (d > 1e-4 && d < best)
        best = d;
    }
    if (best < 1e299)
      out.push_back(1.0 / best);
  }
  return out;
}

std::vector<std::size_t> all_of(std::size_t n) {
  std::vector<std::size_t> m(n);
  for (std::size_t i = 0; i < n; ++i)
    m[i] = i;
  return m;
}

} // namespace

TEST(nearest_neighbours_as_every_pair_finds_them) {
  std::mt19937 rng(9);
  std::uniform_real_distribution<double> u(-0.4, 0.4);
  std::normal_distribution<double> g(0, 0.0005);
  std::vector<std::vector<Vec3>> sets;
  { // a lattice, jittered
    std::vector<Vec3> s;
    for (int h = -9; h <= 9; ++h)
      for (int k = -9; k <= 9; ++k)
        for (int l = -9; l <= 9; ++l)
          s.push_back(
              Vec3{h * 0.02 + g(rng), k * 0.023 + g(rng), l * 0.019 + g(rng)});
    sets.push_back(s);
  }
  { // uniform scatter
    std::vector<Vec3> s;
    for (int i = 0; i < 3000; ++i)
      s.push_back(Vec3{u(rng), u(rng), u(rng)});
    sets.push_back(s);
  }
  { // a plane: z all one
    std::vector<Vec3> s;
    for (int i = 0; i < 2000; ++i)
      s.push_back(Vec3{u(rng), u(rng), 0.1});
    sets.push_back(s);
  }
  { // clusters, and coincident pairs
    std::vector<Vec3> s;
    for (int c = 0; c < 20; ++c) {
      const Vec3 centre{u(rng), u(rng), u(rng)};
      for (int i = 0; i < 60; ++i)
        s.push_back(
            Vec3{centre.x + g(rng), centre.y + g(rng), centre.z + g(rng)});
    }
    for (int i = 0; i < 200; ++i)
      s.push_back(s[static_cast<std::size_t>(i) * 3]);
    sets.push_back(s);
  }
  sets.push_back({Vec3{0, 0, 0}, Vec3{0.1, 0, 0}}); // two
  sets.push_back(
      {Vec3{0, 0, 0}, Vec3{0, 0, 0}, Vec3{0, 0, 0}}); // all coincident
  for (std::size_t k = 0; k < sets.size(); ++k) {
    const std::vector<Vec3> &s = sets[k];
    const std::vector<std::size_t> m = all_of(s.size());
    const std::vector<double> want = every_pair(s, m);
    for (std::size_t threads : {1, 3}) {
      set_parallel_threads(threads);
      check::is_true(inverse_nearest_neighbours(s, m) == want,
                     "set " + std::to_string(k) + ", " +
                         std::to_string(threads) + " thread(s)");
    }
    // A subset, in another order: only members count, in members' order.
    if (s.size() > 10) {
      std::vector<std::size_t> some;
      for (std::size_t i = s.size(); i-- > 0;)
        if (i % 3 != 1)
          some.push_back(i);
      check::is_true(inverse_nearest_neighbours(s, some) == every_pair(s, some),
                     "set " + std::to_string(k) + ", a subset");
    }
  }
  set_parallel_threads(0);
}

} // namespace mxi
