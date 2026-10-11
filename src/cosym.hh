// cosym: several data sets brought to one setting under their lattice's
// symmetry -- Gildea & Winter (2018), Acta Cryst. D74, 405-410, after Brehm &
// Diederichs (2014); docs/multi-lattice.md.

#ifndef MXI_COSYM_HH
#define MXI_COSYM_HH

#include <cstddef>
#include <vector>

#include "laue.hh"
#include "symmetry.hh"

namespace mxi {

struct CosymResult {
  //: The lattice group's rotations, identity first.
  std::vector<Rotation> operations;
  //: Each data set's reindexing, an index into `operations`: its Miller
  //: indices h to apply(operations[k], h). The first data set's is the
  //: identity, the reference.
  std::vector<std::size_t> chosen;
  //: The share of the variance on each principal component, largest first
  //: (the first few): how many clusters there are to see.
  std::vector<double> explained;
  std::size_t dimensions = 0; //: components kept
  //: Each data set's chosen point's cosine with the reference point.
  std::vector<double> cosine;
  //: Each data set's correlation with the first, reindexed, over the
  //: reflections both have; and how many.
  std::vector<double> cc;
  std::vector<std::size_t> common;
};

//: The data sets, each merged in P1 and normalised, and the rotations of their
//: lattice in the first's basis. The correlation of data set i under rotation
//: k with data set j under rotation l, over the reflections both have -- each
//: data set's Friedel mates averaged -- for every i, j, k, l; each a point in
//: as many dimensions as the matrix's leading eigenvectors hold, the matrix's
//: best fit of low rank, which is what the paper's L-BFGS fits where entries
//: are missing; each data set reindexed by the rotation whose point lies
//: nearest in direction to the first's identity point -- in its cluster, so
//: its coset of the true Patterson group -- keeping the identity where it is
//: in that cluster, so a data set in the reference's setting is not reindexed
//: into an equivalent one for nothing.
CosymResult cosym(const std::vector<P1Intensities> &datasets,
                  const std::vector<Rotation> &operations);

} // namespace mxi

#endif // MXI_COSYM_HH
