// Merging scaled observations, and amplitudes from them: what dials.merge
// does (DIALS' dials/algorithms/merging/merge.py and french_wilson.py, read
// 26 November), for mxi merge.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <gemmi/symmetry.hpp>

#include "symmetry.hh"

namespace mxi {

//: Observations merged into the asymmetric unit (gemmi's, the CCP4 one), as
//: cctbx's merge_equivalents with use_internal_variance False, which
//: dials.merge uses: the inverse-variance mean, sigma 1/sqrt(sum w). With
//: Friedel mates apart (anomalous), an acentric reflection's two
//: hemispheres are merged apart; a centric one has one.
struct MergedSet {
  std::vector<Miller> hkl; //: in the asymmetric unit
  std::vector<double> value, sigma;
  std::vector<std::int32_t> n;
  std::vector<bool> centric;
  //: +1 the positive hemisphere, I(+); -1 the negative, I(-); 0 merged.
  std::vector<std::int8_t> hemisphere;
  std::size_t size() const { return hkl.size(); }
};

MergedSet merge_equivalents(const std::vector<Miller> &hkl,
                            const std::vector<double> &intensity,
                            const std::vector<double> &variance,
                            const gemmi::SpaceGroup &group, bool anomalous);

//: Posterior amplitudes, French and Wilson (1978), as dials.merge's own
//: implementation computes them: the expected intensity the standardised
//: median (median over ln 2) in bins of d*^3 of equal count (at most 60, at
//: least 20 a bin), interpolated linearly at each reflection; the posterior
//: moments interpolated in cctbx's tables for h below 3 (acentric) or 4
//: (centric), and their asymptotic forms above; a reflection refused where
//: I/sigma is below -3.7 or h below -4, or its sigma not positive.
struct FrenchWilson {
  std::vector<double> F, sigF, J, sigJ;
  std::vector<bool> valid;
};

FrenchWilson french_wilson(const std::vector<double> &intensity,
                           const std::vector<double> &sigma,
                           const std::vector<bool> &centric,
                           const std::vector<double> &d_star_cubed,
                           std::size_t max_bins = 60);

//: A FreeR flag in CCP4's convention, 0 to groups - 1, 0 the test set: a
//: hash of the reflection's index taken in a canonical form under the
//: lattice's rotations and inversion -- so that reflections a twin law would
//: relate share a flag (as cctbx's use_lattice_symmetry), and the same index
//: has the same flag in every data set of one lattice. Not cctbx's random
//: draw, and no cap of max_free.
std::int32_t free_flag(const Miller &hkl, const std::vector<Rotation> &lattice,
                       std::int32_t groups);

} // namespace mxi
