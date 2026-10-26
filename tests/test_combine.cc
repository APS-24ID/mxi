// mxi_combine's and mxi_split's library half (src/combine.hh): lists joined
// with every model reference offset, reflections renumbered, the same
// experiment refused twice, uncommon columns left out and named; and a data
// set taken out with its sweep, and put back.

#include <cstdint>
#include <string>
#include <vector>

#include "../src/args.hh"
#include "../src/combine.hh"
#include "../src/expt.hh"
#include "../src/index.hh"
#include "../src/json.hh"
#include "../src/postrefine.hh"
#include "../src/refl.hh"
#include "../src/untrusted.hh"
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

// A known cell for indexing (mxi_index --known-cell, item 74): insulin's 78 A
// I-centred cube is, primitive and reduced, a rhombohedron of 67.55 A and
// 109.47 deg; it matches a refined reduced cell in any order of its edges and
// angles, and a wrong one does not.
TEST(a_known_centred_cell_reduces_to_its_primitive_and_matches) {
  KnownCell known;
  known.cell = {78.0, 78.0, 78.0, 90.0, 90.0, 90.0};
  known.centring = 'I';
  const UnitCell r = reduced_primitive(known);
  check::close(r.a, 67.55, 0.01, "the edge");
  check::close(r.alpha, 109.47, 0.01, "the angle");
  const UnitCell refined = {67.428, 67.383, 67.374, 109.466, 109.533, 109.520};
  check::is_true(cells_match(refined, r, known), "a refined cell, reordered");
  const UnitCell wrong = {50.0, 50.0, 50.0, 90.0, 90.0, 90.0};
  check::is_true(!cells_match(wrong, r, known), "a wrong cell");
  KnownCell primitive;
  primitive.cell = {67.4, 67.4, 67.4, 109.47, 109.47, 109.47};
  check::is_true(cells_match(reduced_primitive(primitive), r, known),
                 "the same lattice given primitive");
}

// The real-space grid search (Gildea et al. 2014; mxi_index --real-space-grid):
// a planted lattice in an arbitrary orientation, its reciprocal lattice points
// to 2 A; the search finds each real-space axis, length and direction, among
// its candidates, and the known cell then chooses that basis.
TEST(the_real_space_grid_search_finds_a_planted_lattice) {
  // An orthorhombic cell, turned by an arbitrary rotation.
  const double t = 0.7, u = 1.1, v = 0.4;
  const Mat3 Rx{1, 0,           0,          0, std::cos(t), -std::sin(t),
                0, std::sin(t), std::cos(t)};
  const Mat3 Ry{std::cos(u),  0, std::sin(u), 0, 1, 0,
                -std::sin(u), 0, std::cos(u)};
  const Mat3 Rz{
      std::cos(v), -std::sin(v), 0, std::sin(v), std::cos(v), 0, 0, 0, 1};
  const Mat3 R = Rz * Ry * Rx;
  const Vec3 a = R * Vec3{11.0, 0, 0}, b = R * Vec3{0, 17.0, 0},
             c = R * Vec3{0, 0, 23.0};
  const Mat3 real = Mat3::from_rows(a, b, c);
  const Mat3 recip = real.inverse(); // columns a*, b*, c*
  std::vector<Vec3> points;
  for (int h = -6; h <= 6; ++h)
    for (int k = -9; k <= 9; ++k)
      for (int l = -12; l <= 12; ++l) {
        const Vec3 r = recip * Vec3{double(h), double(k), double(l)};
        if ((h || k || l) && r.norm() < 0.5)
          points.push_back(r);
      }
  KnownCell known;
  known.cell = {11.0, 17.0, 23.0, 90.0, 90.0, 90.0};
  const std::vector<Vec3> found =
      real_space_grid_candidates(points, known, 1.15, 30, 5.0);
  for (const Vec3 &axis : {a, b, c}) {
    double best = 180.0;
    for (const Vec3 &f : found)
      if (std::fabs(f.norm() - axis.norm()) < 0.01 * axis.norm()) {
        const double cosine =
            std::fabs(f.x * axis.x + f.y * axis.y + f.z * axis.z) /
            (f.norm() * axis.norm());
        best = std::fmin(best, std::acos(std::fmin(1.0, cosine)) * 180 / M_PI);
      }
    check::is_true(best < 1.5,
                   "an axis found within 1.5 deg: " + std::to_string(best));
  }
  Crystal crystal;
  std::size_t indexed = 0;
  check::is_true(
      choose_basis(found, points, 0.3, &crystal, &indexed, nullptr, &known),
      "a basis chosen");
  check::is_true(cells_match(crystal.cell(), reduced_primitive(known), known),
                 "the planted cell");
  check::is_true(indexed > points.size() * 9 / 10, "most points indexed");
}

// Untrusted regions, as mxi_mask draws them: a pixel inside when its centre
// is, DIALS's pixel convention; clipped at the frame's edge; a concave polygon
// by the even-odd rule; a shape on another panel refused.
TEST(untrusted_regions_are_the_pixels_whose_centres_they_hold) {
  const auto shapes = [](const char *text) { return json::parse(text); };
  // A 2 x 2 rectangle, corners either way round: four pixels.
  auto r = untrusted_pixels(shapes(R"([{"rectangle": [3, 2, 5, 4]}])"), 10, 10);
  check::equal(static_cast<long long>(r.size()), 4, "a 2 x 2 rectangle");
  check::equal(static_cast<long long>(r.front()), 2 * 10 + 3, "its first");
  check::equal(
      static_cast<long long>(
          untrusted_pixels(shapes(R"([{"rectangle": [5, 4, 3, 2]}])"), 10, 10)
              .size()),
      4, "corners the other way round");
  // A circle of radius 1.5 about a pixel's centre: that pixel and its eight
  // neighbours' centres within 1.5 -- the four edge ones, the corners at
  // 1.41 too.
  auto c = untrusted_pixels(shapes(R"([{"circle": [5.5, 5.5, 1.5]}])"), 10, 10);
  check::equal(static_cast<long long>(c.size()), 9, "nine pixels");
  // Clipped at the edge: a circle about a corner keeps only what is inside.
  auto e = untrusted_pixels(shapes(R"([{"circle": [0, 0, 2]}])"), 10, 10);
  for (std::uint64_t i : e)
    check::is_true(i % 10 < 2 && i / 10 < 2, "inside the frame");
  // A concave polygon, an L: its notch is not inside.
  auto l = untrusted_pixels(
      shapes(R"([{"polygon": [[0,0],[4,0],[4,2],[2,2],[2,4],[0,4]]}])"), 10,
      10);
  check::equal(static_cast<long long>(l.size()), 12, "the L's twelve pixels");
  check::is_true(!std::binary_search(l.begin(), l.end(), 3 * 10 + 3),
                 "the notch left out");
  // Overlaps counted once.
  auto o = untrusted_pixels(
      shapes(R"([{"rectangle": [0,0,2,2]}, {"rectangle": [1,1,3,3]}])"), 10,
      10);
  check::equal(static_cast<long long>(o.size()), 7,
               "four and four, one shared");
  bool refused = false;
  try {
    untrusted_pixels(shapes(R"([{"panel": 1, "circle": [1, 1, 1]}])"), 10, 10);
  } catch (const std::runtime_error &) {
    refused = true;
  }
  check::is_true(refused, "another panel refused");
  check::is_true(describe_untrusted(shapes(R"([{"circle": [1,1,1]},
      {"polygon": [[0,0],[1,0],[0,1]]}, {"circle": [2,2,1]}])")) ==
                     "2 circles, 1 polygon",
                 "described");
}

// An option given more than once keeps every value, in order, for mxi_scale's
// --exclude; value() keeps the last, as it always has.
TEST(a_repeated_option_keeps_every_value) {
  const char *argv[] = {"prog",    "--exclude", "6:600:800",
                        "in.rflx", "--exclude", "0:550-700"};
  const Arguments a = parse_arguments(6, const_cast<char **>(argv),
                                      {"--exclude"}, {"--exclude"});
  check::is_true(a.ok, "parsed");
  check::equal(static_cast<long long>(a.values("--exclude").size()), 2, "both");
  check::is_true(a.values("--exclude")[0] == "6:600:800" &&
                     a.values("--exclude")[1] == "0:550-700",
                 "in order");
  check::is_true(a.value("--exclude") == "0:550-700", "value() the last");
  check::is_true(a.positional.size() == 1, "the file");
}
