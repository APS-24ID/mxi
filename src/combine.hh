// Experiment lists and their reflections split apart and joined together:
// mxi_split and mxi_combine, and mxi_correlation's clusters.
//
// An experiment list refers to its models -- beam, crystal, detector,
// goniometer, image set, profile, scaling model, scan -- by index into arrays
// of them; the reflections refer to their experiments by `id`, and to their
// identifiers through the table's map. Splitting keeps the models and
// renumbers the experiments; joining appends each list's models and offsets
// every reference by what came before, as dials.combine_experiments does.

#ifndef MXI_COMBINE_HH
#define MXI_COMBINE_HH

#include <string>
#include <vector>

#include "json.hh"
#include "refl.hh"

namespace mxi {

//: An experiment list and its reflections.
struct ExperimentsAndReflections {
  json::Value document;
  Table reflections;
};

//: The experiments `keep` marks (by index), each with the crystal-less sweep
//: experiment it shares an image set and scan with, if there is one -- without
//: it the data set could not be integrated or converted again -- and their
//: reflections, the ids and identifiers renumbered to match, the shoeboxes
//: with them -- integration's profile model wants a sweep's own spots'.
ExperimentsAndReflections select_experiments(const json::Value &document,
                                             const Table &reflections,
                                             const std::vector<bool> &keep);

//: Several experiment lists and their reflections joined: each list's models
//: appended, its experiments' references to them offset, its reflections'
//: ids and identifiers offset by the experiments before. The same experiment
//: twice -- the same identifier in two inputs -- is refused. The columns kept
//: are those every table has, of the same type, and the shoeboxes where every
//: input has them; any other is left out and named in `dropped`.
struct Combined : ExperimentsAndReflections {
  std::vector<std::string> dropped;
};
Combined combine_experiments(const std::vector<json::Value> &documents,
                             const std::vector<Table> &tables);

//: The experiments in groups, each to be refined on its own: those on one
//: sweep (one image set and scan) together -- several lattices on one sweep
//: share its detector -- and, `by_models`, also any sharing a detector, a beam
//: or a crystal model. mxi_refine groups by sweep alone: each sweep refined on
//: its own, a shared beam or detector a copy each, gave better results even
//: for one crystal's sweeps (Graeme, 18 October); --joint for one fit.
std::vector<std::vector<std::size_t>>
independent_groups(const json::Value &document, bool by_models = false);

} // namespace mxi

#endif // MXI_COMBINE_HH
