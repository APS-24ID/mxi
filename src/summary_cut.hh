// A processing summary cut down to some of its experiments: a cluster's, as
// mxi correlation --write-clusters writes each. Carried whole, a cluster's
// summary was every crystal's -- its report the whole run's (Graeme,
// 5 November).
//
// Of each step's record: the parts a sweep or an experiment at a time --
// tables with a "sweep" column, lists of entries with "sweep" or
// "experiment", refinement's groups by their "experiments", spot finding's
// spots_per_image and first_image, import's table of sweeps -- kept for the
// experiments kept, renumbered; spot finding's spots and signal pixels the
// sums of the sweeps kept; the rest, figures of the whole run that cannot be
// cut, left out and named in "whole_run_left_out". Symmetry's record whole:
// its group and reindexing are the cluster's too. Scaling's left out: its
// statistics are of every crystal, and a cluster is scaled afresh.

#ifndef MXI_SUMMARY_CUT_HH
#define MXI_SUMMARY_CUT_HH

#include <vector>

#include "json.hh"

namespace mxi {

//: The summary {"steps": [...]} cut to the experiments `keep` marks of the
//: experiment list `document`, renumbered as the kept experiments are, in
//: order; their sweeps numbered as their image sets first appear. `own`, if
//: not null, the record of the step that cut it, appended.
json::Value cut_summary(const json::Value &summary, const json::Value &document,
                        const std::vector<bool> &keep,
                        const json::Value *own = nullptr);

} // namespace mxi

#endif // MXI_SUMMARY_CUT_HH
