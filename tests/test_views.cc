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
