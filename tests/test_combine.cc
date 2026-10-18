// mxi_combine's and mxi_split's library half (src/combine.hh): lists joined
// with every model reference offset, reflections renumbered, the same
// experiment refused twice, uncommon columns left out and named; and a data
// set taken out with its sweep, and put back.

#include <cstdint>
#include <string>
#include <vector>

#include "../src/combine.hh"
#include "../src/expt.hh"
#include "../src/json.hh"
#include "../src/postrefine.hh"
#include "../src/refl.hh"
#include "check.hh"

using namespace mxi;

namespace {

// A sweep (crystal-less) and a data set of it: two experiments sharing an
// image set and a scan, with one of each model and the given identifiers.
json::Value list(const std::string &sweep, const std::string &data) {
  return json::parse(R"({"__id__": "ExperimentList",
    "experiment": [
      {"beam": 0, "detector": 0, "goniometer": 0, "scan": 0, "imageset": 0,
       "identifier": ")" +
                     sweep + R"("},
      {"beam": 0, "detector": 0, "goniometer": 0, "scan": 0, "imageset": 0,
       "crystal": 0, "identifier": ")" +
                     data + R"("}],
    "beam": [{"b": 1}], "detector": [{"d": 1}], "goniometer": [{"g": 1}],
    "scan": [{"s": 1}], "imageset": [{"i": 1}], "crystal": [{"c": 1}],
    "profile": [], "scaling_model": []})");
}

// Rows: one unindexed spot of the sweep, two reflections of the data set.
Table rows(const std::string &sweep, const std::string &data, bool extra) {
  Table t;
  t.nrows = 3;
  Column &id = t.int_column("id", "int", 1);
  id.ints = {0, 1, 1};
  Column &x = t.real_column("intensity.sum.value", "double", 1);
  x.reals = {1.0, 2.0, 3.0};
  if (extra)
    t.real_column("only.here", "double", 1).reals = {0.0, 0.0, 0.0};
  t.identifiers = {{0, sweep}, {1, data}};
  return t;
}

} // namespace

TEST(combine_offsets_every_reference_and_renumbers_the_reflections) {
  const Combined c = combine_experiments(
      {list("s-a", "d-a"), list("s-b", "d-b")},
      {rows("s-a", "d-a", false), rows("s-b", "d-b", true)});
  const json::Array &e = c.document["experiment"].as_array();
  check::equal(static_cast<long long>(e.size()), 4, "four experiments");
  check::equal(static_cast<long long>(e[3]["crystal"].as_number()), 1,
               "the second data set's crystal, offset by the first's");
  check::equal(static_cast<long long>(e[2]["imageset"].as_number()), 1,
               "the second sweep's image set, offset");
  check::equal(static_cast<long long>(c.document["crystal"].as_array().size()),
               2, "both crystals");
  check::is_true(c.reflections.at("id").ints ==
                     std::vector<std::int64_t>{0, 1, 1, 2, 3, 3},
                 "the second input's reflections numbered after the first's");
  check::is_true(c.reflections.identifiers.at(3) == "d-b",
                 "and its identifiers with them");
  check::is_true(c.dropped == std::vector<std::string>{"only.here"},
                 "a column not in every input, left out and named");
}

TEST(combine_refuses_the_same_experiment_twice) {
  bool refused = false;
  try {
    combine_experiments({list("s-a", "d-a"), list("s-b", "d-a")},
                        {rows("s-a", "d-a", false), rows("s-b", "d-a", false)});
  } catch (const std::exception &) {
    refused = true;
  }
  check::is_true(refused, "one identifier in two inputs");
}

TEST(a_data_set_taken_out_with_its_sweep_and_put_back) {
  const Combined both = combine_experiments(
      {list("s-a", "d-a"), list("s-b", "d-b")},
      {rows("s-a", "d-a", false), rows("s-b", "d-b", false)});
  // The second data set alone: experiment 3, which brings its sweep, 2.
  const ExperimentsAndReflections b = select_experiments(
      both.document, both.reflections, {false, false, false, true});
  const json::Array &e = b.document["experiment"].as_array();
  check::equal(static_cast<long long>(e.size()), 2,
               "the data set and its sweep");
  check::is_true(e[0]["identifier"].as_string() == "s-b" &&
                     e[1]["identifier"].as_string() == "d-b",
                 "the second input's two");
  check::is_true(b.reflections.at("id").ints ==
                     std::vector<std::int64_t>{0, 1, 1},
                 "its rows, numbered from 0 again");
  check::is_true(b.reflections.at("intensity.sum.value").reals ==
                     std::vector<double>{1.0, 2.0, 3.0},
                 "the values unchanged");
}

// Where a list's images are is its experiments' image set, not the document's
// first: a split or combined list keeps every image set, and taking the first
// integrated one data set from another's file (17 October).
TEST(a_lists_images_are_its_experiments_own_image_set) {
  const json::Value document = json::parse(R"({"__id__": "ExperimentList",
    "experiment": [{"beam": 0, "detector": 0, "goniometer": 0, "scan": 0,
                    "imageset": 1}],
    "beam": [], "detector": [], "goniometer": [], "scan": [],
    "imageset": [{"__id__": "ImageSequence", "template": "/a/first.nxs"},
                 {"__id__": "ImageSequence", "template": "/b/second.nxs"}],
    "crystal": [], "profile": [], "scaling_model": []})");
  const std::string got = experiments_from_json(document).image_template;
  check::is_true(got == "/b/second.nxs",
                 "the second file, its experiment's: " + got);
}

// Which experiments refinement must fit together: one sweep's, or any sharing
// a detector, beam or crystal. Two data sets sharing one beam had been refined
// as one, each pulled from where it lies alone (17 October).
TEST(experiments_sharing_nothing_refinement_fits_are_independent) {
  // Two sweeps, each a crystal-less experiment and a lattice, nothing shared.
  json::Value apart = json::parse(R"({"experiment": [
    {"beam": 0, "detector": 0, "imageset": 0, "scan": 0},
    {"beam": 0, "detector": 0, "imageset": 0, "scan": 0, "crystal": 0},
    {"beam": 1, "detector": 1, "imageset": 1, "scan": 1},
    {"beam": 1, "detector": 1, "imageset": 1, "scan": 1, "crystal": 1}]})");
  const auto g = independent_groups(apart);
  check::equal(static_cast<long long>(g.size()), 2, "two groups");
  check::is_true(g[0] == std::vector<std::size_t>{0, 1} &&
                     g[1] == std::vector<std::size_t>{2, 3},
                 "each sweep's two together");
  // The same, the beam shared: one group.
  json::Value shared = json::parse(R"({"experiment": [
    {"beam": 0, "detector": 0, "imageset": 0, "scan": 0, "crystal": 0},
    {"beam": 0, "detector": 1, "imageset": 1, "scan": 1, "crystal": 1}]})");
  check::equal(static_cast<long long>(independent_groups(shared).size()), 2,
               "by sweep, the default: a shared beam does not join them");
  check::equal(static_cast<long long>(independent_groups(shared, true).size()),
               1, "by models: a shared beam joins them");
  // One crystal across two sweeps: one group.
  json::Value crystal = json::parse(R"({"experiment": [
    {"beam": 0, "detector": 0, "imageset": 0, "scan": 0, "crystal": 0},
    {"beam": 1, "detector": 1, "imageset": 1, "scan": 1, "crystal": 0}]})");
  check::equal(static_cast<long long>(independent_groups(crystal).size()), 2,
               "by sweep: a shared crystal does not join them");
  check::equal(static_cast<long long>(independent_groups(crystal, true).size()),
               1, "by models: a shared crystal joins them");
  // Two lattices on one sweep: together either way, its detector theirs.
  json::Value lattices = json::parse(R"({"experiment": [
    {"beam": 0, "detector": 0, "imageset": 0, "scan": 0, "crystal": 0},
    {"beam": 0, "detector": 0, "imageset": 0, "scan": 0, "crystal": 1}]})");
  check::equal(static_cast<long long>(independent_groups(lattices).size()), 1,
               "one sweep's two lattices together");
}

// Control points one per `interval` degrees of a scan, at least `minimum`:
// the default one per 10, at least five; mxi_refine --interval 30 on a 170
// degree sweep, six -- for sparse data, a small molecule's (18 October).
TEST(control_points_follow_the_interval) {
  Scan scan;
  scan.osc_width = 0.1;
  scan.first_image = 1;
  scan.last_image = 1700; // 170 degrees
  check::equal(static_cast<long long>(scan_varying_points(scan)), 17,
               "one per 10 degrees");
  check::equal(static_cast<long long>(scan_varying_points(scan, 30.0, 4)), 6,
               "one per 30");
  scan.last_image = 300; // 30 degrees
  check::equal(static_cast<long long>(scan_varying_points(scan)), 5,
               "at least five by default");
  check::equal(static_cast<long long>(scan_varying_points(scan, 30.0, 4)), 4,
               "at least four with an interval");
}
