// The two views of an experiment list (docs/multi-lattice.md, stage 2): a
// sweep's crystal-less experiment and its lattices, split for a program and
// put back.

#include <vector>

#include "../src/views.hh"
#include "check.hh"

namespace mxi {

namespace {

Experiment lattice(double cell) {
  Experiment e;
  e.scan.osc_width = 0.1;
  e.crystal = Crystal::from_real_space({cell, 0.0, 0.0}, {0.0, cell, 0.0},
                                       {0.0, 0.0, cell});
  Panel p;
  p.origin = {-100.0, 100.0, -200.0};
  e.detector.panels.push_back(p);
  return e;
}

// A sweep's own experiment, crystal-less, then two lattices of it; all three
// sharing one detector and one beam, as one sweep's do.
ExperimentList one_sweep() {
  ExperimentList whole;
  Experiment sweep = lattice(50.0);
  sweep.crystal.reset();
  sweep.identifier = "sweep";
  Experiment a = lattice(50.0), b = lattice(60.0);
  a.identifier = "a";
  b.identifier = "b";
  whole.experiments = {sweep, a, b};
  whole.beam_of = {0, 0, 0};
  whole.detector_of = {0, 0, 0};
  whole.sweep_of = {0, 0, 0};
  return whole;
}

// An unindexed spot of the sweep, then one of each lattice.
Table spots() {
  Table t;
  t.nrows = 3;
  Column &id = t.int_column("id", "int", 1);
  id.ints = {0, 1, 2};
  Column &x = t.real_column("xyzobs.px.value", "vec3<double>", 3);
  x.reals = {1, 1, 1, 2, 2, 2, 3, 3, 3};
  t.identifiers = {{0, "sweep"}, {1, "a"}, {2, "b"}};
  return t;
}

} // namespace

TEST(a_sweeps_lattices_are_taken_out_and_put_back) {
  const ExperimentList whole = one_sweep();
  const Table table = spots();
  const Lattices l = lattices_of(whole);
  check::equal(static_cast<long long>(l.list.size()), 2, "two lattices");
  check::is_true(l.position == std::vector<std::size_t>{1, 2} &&
                     l.crystal_less == std::vector<std::size_t>{0},
                 "in their places, the sweep's own apart");
  const Table rows = lattice_rows(table, l);
  check::equal(static_cast<long long>(rows.nrows), 2, "the lattices' rows");
  check::is_true(rows.at("id").ints == std::vector<std::int64_t>{0, 1} &&
                     rows.identifiers.at(0) == "a" &&
                     rows.identifiers.at(1) == "b",
                 "numbered as the lattices are");
  const Table rest = crystal_less_rows(table, l);
  check::equal(static_cast<long long>(rest.nrows), 1, "the unindexed spot");

  // A program: the shared detector moved, a column added.
  ExperimentList done = l.list;
  for (Experiment &e : done.experiments)
    e.detector.panels[0].origin.z = -201.0;
  Table out = rows;
  Column &cal = out.real_column("xyzcal.px", "vec3<double>", 3);
  cal.reals = {9, 9, 9, 8, 8, 8};

  const ExperimentList back = whole_again(whole, l, done);
  check::equal(static_cast<long long>(back.size()), 3, "all three again");
  check::is_true(!back[0].crystal && back[1].crystal && back[2].crystal,
                 "the sweep's own still crystal-less, first");
  check::close(back[0].detector.panels[0].origin.z, -201.0, 1e-12,
               "with the detector the lattices refined");
  const Table joined = rows_again(out, rest, l, back);
  check::equal(static_cast<long long>(joined.nrows), 3, "every row again");
  check::is_true(joined.at("id").ints == std::vector<std::int64_t>{1, 2, 0},
                 "numbered as the whole list is");
  check::close(joined.at("xyzcal.px").real(2, 0), 0.0, 0.0,
               "the unindexed spot no prediction");
  check::is_true(joined.identifiers.size() == 3 &&
                     joined.identifiers.at(0) == "sweep",
                 "every experiment's identifier");
}

TEST(a_list_of_lattices_alone_passes_through) {
  ExperimentList whole = one_sweep();
  whole.experiments.erase(whole.experiments.begin());
  whole.beam_of = whole.detector_of = whole.sweep_of = {0, 0};
  Table table = spots();
  table.nrows = 2;
  table.int_column("id", "int", 1).ints = {0, 1};
  table.real_column("xyzobs.px.value", "vec3<double>", 3).reals = {2, 2, 2,
                                                                   3, 3, 3};
  table.identifiers = {{0, "a"}, {1, "b"}};
  const Lattices l = lattices_of(whole);
  check::is_true(!l.any_crystal_less(), "nothing crystal-less");
  const Table rows = lattice_rows(table, l);
  check::is_true(rows.at("id").ints == table.at("id").ints &&
                     rows.nrows == table.nrows,
                 "the table as it was");
  check::equal(static_cast<long long>(whole_again(whole, l, l.list).size()), 2,
               "the list as it was");
}

TEST(the_experiments_of_one_image_set_and_scan_are_one_sweep) {
  ExperimentList whole = one_sweep();
  whole.sweep_of.clear();
  const json::Value written = experiments_to_json(whole);
  const auto on = [](const json::Value &e, long set, const char *name) {
    json::Object o = e.as_object();
    o["imageset"] = json::Value(set);
    if (name != nullptr)
      o["identifier"] = json::Value(name);
    return json::Value(std::move(o));
  };
  json::Array x;
  for (const json::Value &e : written["experiment"].as_array())
    x.push_back(on(e, 0, nullptr));
  // A fourth: another sweep, its own image set.
  x.push_back(on(written["experiment"].as_array()[1], 1, "other"));
  json::Object o = written.as_object();
  o["experiment"] = json::Value(std::move(x));
  const json::Value doc(std::move(o));
  const ExperimentList back = experiments_from_json(doc);
  check::is_true(back.sweep_of == std::vector<std::size_t>{0, 0, 0, 1},
                 "three of one sweep, one of another");
}

} // namespace mxi

namespace mxi {

TEST(indexing_carries_the_sweep_given_through_and_adds_its_lattices) {
  // dials/dials#1029, as mxi_index --keep-sweeps writes it: the sweep given,
  // crystal-less, first, with the detector its lattices refined; the
  // lattices after it, the first's identifier -- the sweep's, as indexing
  // gives it -- made fresh; an indexed spot its lattice's, an unindexed one
  // its sweep's.
  ExperimentList given = one_sweep();
  given.experiments.resize(1); // the sweep alone, as imported
  given.beam_of = given.detector_of = given.sweep_of = {0};
  ExperimentList found;
  Experiment a = lattice(50.0), b = lattice(60.0);
  a.identifier = "sweep"; // the first lattice keeps the sweep's, in indexing
  b.identifier = "b";
  a.detector.panels[0].origin.z = b.detector.panels[0].origin.z = -201.0;
  found.experiments = {a, b};
  found.beam_of = found.detector_of = {0, 0};
  Table t;
  t.nrows = 3;
  t.int_column("id", "int", 1).ints = {0, 1, 0};
  t.int_column("miller_index", "cctbx::miller::index<>", 3).ints = {
      1, 2, 3, 4, 5, 6, 0, 0, 0};
  const ExperimentList whole = add_sweeps(given, found, {0, 0}, &t);
  check::equal(static_cast<long long>(whole.size()), 3, "one in, two found");
  check::is_true(!whole[0].crystal && whole[1].crystal && whole[2].crystal,
                 "the sweep first, crystal-less");
  check::is_true(whole[0].identifier == "sweep" &&
                     whole[1].identifier != "sweep" &&
                     whole[2].identifier == "b",
                 "the sweep's identifier its own, the lattice's made fresh");
  check::close(whole[0].detector.panels[0].origin.z, -201.0, 1e-12,
               "the detector its lattices refined");
  check::is_true(whole.detector_of[0] == whole.detector_of[1] &&
                     whole.detector_of[1] == whole.detector_of[2] &&
                     whole.sweep_of == std::vector<std::size_t>{0, 0, 0},
                 "one detector, one sweep");
  check::is_true(t.at("id").ints == std::vector<std::int64_t>{1, 2, 0},
                 "indexed spots their lattices', the unindexed the sweep's");
  check::is_true(t.identifiers.size() == 3 &&
                     t.identifiers.at(1) == whole[1].identifier,
                 "every experiment's identifier");
}

} // namespace mxi

namespace mxi {

TEST(dials_form_leaves_the_crystal_less_out_and_unindexed_spots_minus_one) {
  // mxi_convert's pair for DIALS: the sweep's crystal-less experiment left
  // out, the lattices renumbered from 0 with their references unchanged, an
  // unindexed spot -1 as dials.index leaves it, the identifiers to match.
  const ExperimentList whole = one_sweep();
  const json::Value doc = experiments_to_json(whole);
  Table t = spots();
  const json::Value dials = dials_form(doc, &t);
  const json::Array &x = dials["experiment"].as_array();
  check::equal(static_cast<long long>(x.size()), 2, "the two lattices");
  check::is_true(
      x[0]["crystal"].is_number() && x[1]["crystal"].is_number() &&
          x[0]["detector"].as_number() ==
              doc["experiment"].as_array()[1]["detector"].as_number(),
      "their references unchanged");
  check::is_true(t.at("id").ints == std::vector<std::int64_t>{-1, 0, 1},
                 "the unindexed spot -1, the lattices' 0 and 1");
  check::is_true(t.identifiers.size() == 2 && t.identifiers.at(0) == "a" &&
                     t.identifiers.at(1) == "b",
                 "the identifiers renumbered");
  // A list of lattices alone: as it is.
  ExperimentList lattices = whole;
  lattices.experiments.erase(lattices.experiments.begin());
  lattices.beam_of = lattices.detector_of = lattices.sweep_of = {0, 0};
  const json::Value plain = experiments_to_json(lattices);
  Table u = spots();
  const json::Value same = dials_form(plain, &u);
  check::equal(static_cast<long long>(same["experiment"].as_array().size()), 2,
               "nothing left out");
  check::is_true(u.at("id").ints == std::vector<std::int64_t>{0, 1, 2},
                 "the table as it was");
  // An imported list -- crystal-less alone, before indexing -- as it is: the
  // first form mxi_convert wrote it with no experiments at all.
  ExperimentList imported = whole;
  imported.experiments.resize(1);
  imported.beam_of = imported.detector_of = imported.sweep_of = {0};
  const json::Value before = experiments_to_json(imported);
  check::equal(static_cast<long long>(
                   dials_form(before, nullptr)["experiment"].as_array().size()),
               1, "an imported list kept whole");
}

} // namespace mxi
