// The two views of an experiment list (docs/multi-lattice.md): its lattices,
// the experiments with a crystal, which every program that predicts takes; and
// its sweeps, each a crystal-less experiment carried through from import, the
// experiment its unindexed spots belong to (dials/dials#1029). Made here, in
// one place, so that no program asks of an experiment whether it has a
// crystal: it takes the lattices, and gives them back.

#ifndef MXI_VIEWS_HH
#define MXI_VIEWS_HH

#include <cstddef>
#include <vector>

#include "expt.hh"
#include "refl.hh"

namespace mxi {

struct Lattices {
  //: The experiments with a crystal, in order, as a list of their own -- their
  //: sharing, profile models and source kept -- for a program to work on.
  ExperimentList list;
  //: Each one's place in the whole list.
  std::vector<std::size_t> position;
  //: The whole list's crystal-less experiments, its sweeps' own.
  std::vector<std::size_t> crystal_less;
  //: Whether the whole list has any: a list of lattices alone, as before,
  //: passes through every function here unchanged.
  bool any_crystal_less() const { return !crystal_less.empty(); }
};

//: The whole list's lattices.
Lattices lattices_of(const ExperimentList &whole);

//: A table's rows of the lattices, `id` and the identifier map in the
//: lattices' own numbering; and its other rows -- the crystal-less
//: experiments', the unindexed spots -- as they are.
Table lattice_rows(const Table &table, const Lattices &lattices);
Table crystal_less_rows(const Table &table, const Lattices &lattices);

//: The whole list again, the lattices as a program left them put back in
//: their places, and each crystal-less experiment given any model it shares
//: with a lattice of its sweep -- a detector refined once for the sweep is
//: the sweep's detector -- as the lattice now has it.
ExperimentList whole_again(const ExperimentList &whole,
                           const Lattices &lattices,
                           const ExperimentList &processed);

//: The lattices' rows, as a program left them, with the crystal-less
//: experiments' after them: `id` and the identifier map in the whole list's
//: numbering again, and any column the lattices' rows have and the others do
//: not given them as zeros, as an unindexed spot has no prediction.
Table rows_again(const Table &lattice_rows, const Table &crystal_less_rows,
                 const Lattices &lattices, const ExperimentList &whole);

} // namespace mxi

#endif // MXI_VIEWS_HH
