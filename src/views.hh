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

//: The whole list of a program that finds lattices -- indexing -- as
//: dials/dials#1029 has it: every experiment given carried through,
//: crystal-less, first and in its place, with its identifier but the models
//: its sweep's lattices refined; each lattice after them, a fresh identifier
//: if it had one given, referring to its sweep's shared models; so the count
//: out is the count in plus the lattices. `sweep_of_lattice` is each
//: lattice's experiment given. `table`'s ids made the whole list's: an
//: indexed spot's its lattice's, an unindexed one's -- Miller index 0 0 0 --
//: its sweep's.
ExperimentList add_sweeps(const ExperimentList &given,
                          const ExperimentList &lattices,
                          const std::vector<std::size_t> &sweep_of_lattice,
                          Table *table);

//: An experiment list, as a document, in the form DIALS reads: its
//: crystal-less experiments left out -- every model, image set and profile
//: still there, each kept experiment's references to them unchanged -- and
//: `table`, if given, numbered to match: a lattice's spots its new index, a
//: crystal-less experiment's -1, as dials.index leaves an unindexed spot; the
//: identifier map renumbered. A list without crystal-less experiments, as it
//: is; and one without lattices -- imported, before indexing -- likewise.
//: mxi_convert's, writing the pair (docs/multi-lattice.md, stage 2).
json::Value dials_form(const json::Value &document, Table *table);

} // namespace mxi

#endif // MXI_VIEWS_HH
