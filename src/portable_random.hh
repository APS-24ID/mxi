// Random choices the same on every platform. std::mt19937 gives the same
// numbers everywhere, but std::shuffle and std::uniform_int_distribution are
// each library's own algorithms: libc++ (macOS) and libstdc++ (Linux) dealt
// different observations from the same seed, so mxi scale fitted its model
// on different subsets on a Mac and on Linux -- and on Graeme's Mac drew a
// subset that broke the fit (24 November) where Linux's did not. These are
// written out here, so a seed means one subset everywhere.
#pragma once

#include <cstdint>
#include <iterator>
#include <random>
#include <utility>

namespace mxi {

//: An integer in [0, n), n at most 2^32, without bias: mt19937's 32-bit
//: draws below the largest multiple of n taken modulo n, the rest drawn
//: again.
inline std::uint64_t random_below(std::mt19937 &rng, std::uint64_t n) {
  if (n <= 1)
    return 0;
  const std::uint64_t span = std::uint64_t{1} << 32;
  const std::uint64_t limit = span - span % n;
  for (;;) {
    const std::uint64_t x = static_cast<std::uint64_t>(rng()) & 0xffffffffu;
    if (x < limit)
      return x % n;
  }
}

//: Fisher and Yates's shuffle, Durstenfeld's form, with random_below.
template <typename It>
void portable_shuffle(It first, It last, std::mt19937 &rng) {
  const auto n = static_cast<std::uint64_t>(std::distance(first, last));
  for (std::uint64_t i = n; i > 1; --i) {
    const std::uint64_t j = random_below(rng, i);
    using std::swap;
    swap(*(first + static_cast<std::ptrdiff_t>(i - 1)),
         *(first + static_cast<std::ptrdiff_t>(j)));
  }
}

} // namespace mxi
