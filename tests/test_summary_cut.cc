// cut_summary (src/summary_cut.hh): a cluster's summary, its own crystals'.

#include "../src/json.hh"
#include "../src/summary_cut.hh"
#include "check.hh"

namespace mxi {

namespace {

json::Value parse(const char *text) { return json::parse(text); }

// Three sweeps, each its crystal-less experiment then its lattice: 0..5.
const char *kDocument = R"({"experiment": [
  {"imageset": 0}, {"imageset": 0, "crystal": 0},
  {"imageset": 1}, {"imageset": 1, "crystal": 1},
  {"imageset": 2}, {"imageset": 2, "crystal": 2}]})";

const char *kSummary = R"({"steps": [
  {"program": "mxi_import", "time": "t0",
   "sweeps": {"columns": ["file", "first_image"], "file": ["a", "b", "c"],
              "first_image": [1, 1, 1]}},
  {"program": "mxi_find", "time": "t1", "spots": 60, "signal_pixels": 600,
   "first_image": [1, 1, 1], "spots_per_image": [[1, 2], [3, 4], [5, 6]],
   "each_sweep": {"columns": ["sweep", "spots", "signal_pixels"],
                  "sweep": [0, 1, 2], "spots": [10, 20, 30],
                  "signal_pixels": [100, 200, 300]},
   "spots_by_resolution": {"columns": ["d_min", "spots"], "d_min": [2.0], "spots": [60]}},
  {"program": "mxi_refine", "time": "t2", "apart": 3,
   "groups": [{"program": "mxi_refine", "experiments": [0, 1]},
              {"program": "mxi_refine", "experiments": [2, 3]},
              {"program": "mxi_refine", "experiments": [4, 5]}],
   "cell_over_scan": [{"experiment": 1, "a": [1]}, {"experiment": 3, "a": [2]},
                      {"experiment": 5, "a": [3]}],
   "cycles": {"columns": ["cycle"], "cycle": [1]}},
  {"program": "mxi_integrate", "time": "t3",
   "each_sweep": [{"sweep": 0, "n": 1}, {"sweep": 1, "n": 2}, {"sweep": 2, "n": 3}],
   "by_resolution": {"columns": ["d_min"], "d_min": [2.0]}},
  {"program": "mxi_symmetry", "time": "t4", "space_group": "I 2 3"},
  {"program": "mxi_scale", "time": "t5", "unique": 9}]})";

} // namespace

TEST(a_clusters_summary_keeps_its_crystals_and_leaves_the_rest) {
  const json::Value document = parse(kDocument);
  // Sweeps 0 and 2: their crystal-less experiments and lattices.
  const std::vector<bool> keep = {true, true, false, false, true, true};
  const json::Value own =
      parse(R"({"program": "mxi_correlation", "cluster": 1})");
  const json::Value cut = cut_summary(parse(kSummary), document, keep, &own);
  const json::Array &steps = cut["steps"].as_array();
  check::equal(static_cast<long long>(steps.size()), 6LL,
               "scaling left out, ours added");
  check::is_true(steps.back()["program"].as_string() == "mxi_correlation",
                 "its own last");

  const json::Value &import = steps[0];
  check::is_true(import["sweeps"]["file"][0].as_string() == "a" &&
                     import["sweeps"]["file"][1].as_string() == "c" &&
                     import["sweeps"]["file"].size() == 2,
                 "import's sweeps 0 and 2");

  const json::Value &find = steps[1];
  check::equal(static_cast<long long>(find["spots"].as_number()), 40LL,
               "10 + 30");
  check::equal(static_cast<long long>(find["signal_pixels"].as_number()), 400LL,
               "100 + 300");
  check::equal(static_cast<long long>(find["spots_per_image"].size()), 2LL,
               "two sweeps");
  check::equal(find["spots_per_image"][1][0].as_number(), 5.0,
               "the third's, now second");
  check::equal(find["each_sweep"]["sweep"][1].as_number(), 1.0,
               "renumbered: 2 to 1");
  check::is_true(!find.contains("spots_by_resolution"),
                 "the whole run's: left out");
  check::is_true(find["whole_run_left_out"][0].as_string() ==
                     "spots_by_resolution",
                 "and said so");

  const json::Value &refine = steps[2];
  check::equal(static_cast<long long>(refine["groups"].size()), 2LL,
               "two groups");
  check::equal(refine["groups"][1]["experiments"][1].as_number(), 3.0,
               "5 to 3");
  check::equal(static_cast<long long>(refine["cell_over_scan"].size()), 2LL,
               "two crystals");
  check::equal(refine["cell_over_scan"][1]["experiment"].as_number(), 3.0,
               "5 to 3");
  check::equal(refine["cell_over_scan"][1]["a"][0].as_number(), 3.0,
               "its own cell");
  check::is_true(!refine.contains("cycles") && !refine.contains("apart"),
                 "whole-run");

  const json::Value &integrate = steps[3];
  check::equal(static_cast<long long>(integrate["each_sweep"].size()), 2LL,
               "two sweeps");
  check::equal(integrate["each_sweep"][1]["n"].as_number(), 3.0, "the third's");
  check::equal(integrate["each_sweep"][1]["sweep"].as_number(), 1.0,
               "renumbered");
  check::is_true(!integrate.contains("by_resolution"), "whole-run");

  check::is_true(steps[4]["space_group"].as_string() == "I 2 3",
                 "symmetry whole");
}

} // namespace mxi
