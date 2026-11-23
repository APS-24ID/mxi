// The beam centre found from the spots (src/beam_search.hh), on an exact
// lattice: a cubic 78 A crystal, 170 mm, 0.9537 A, 40 degrees -- as insulin's
// geometry -- its reflections predicted, their positions the strong spots.
// The detector moved off by a known amount across the beam, the search must
// move it back. Within half the alias spacing: a detector moved by D lambda
// |g| shifts the reciprocal lattice by its own vector g, about 3 mm here,
// which only the Ewald sphere's curvature tells from the truth.
#include <cmath>
#include <vector>

#include "../src/beam_search.hh"
#include "../src/predict.hh"
#include "check.hh"

namespace mxi {
namespace {

Experiment lattice_experiment() {
  Experiment e;
  e.beam.direction = {0.0, 0.0, 1.0};
  e.beam.wavelength = 0.9537;
  Panel p;
  p.fast = {1.0, 0.0, 0.0};
  p.slow = {0.0, -1.0, 0.0};
  p.pixel_size[0] = p.pixel_size[1] = 0.075;
  p.image_size[0] = 2068;
  p.image_size[1] = 2162;
  p.origin = {-0.075 * 1034, 0.075 * 1081, -170.0};
  e.detector.panels.push_back(p);
  e.goniometer.axis = {1.0, 0.0, 0.0};
  e.scan.first_image = 1;
  e.scan.last_image = 400;
  e.scan.osc_start = 0.0;
  e.scan.osc_width = 0.1;
  const Mat3 r = rotation({0.3, -0.5, 0.81}, 0.7);
  e.crystal = Crystal::from_real_space(r * Vec3{78, 0, 0}, r * Vec3{0, 78, 0},
                                       r * Vec3{0, 0, 78});
  return e;
}

// The reflections to 2 A where they are predicted: the strong spots.
Table spots_of(const Experiment &e) {
  PredictOptions po;
  po.d_min = 2.0;
  const std::vector<Prediction> predictions = predict(e, po);
  Table t;
  t.nrows = predictions.size();
  Column &xyz = t.real_column("xyzobs.px.value", "vec3<double>", 3);
  Column &panel = t.int_column("panel", "std::size_t", 1);
  Column &id = t.int_column("id", "int", 1);
  for (std::size_t i = 0; i < predictions.size(); ++i) {
    xyz.reals[i * 3 + 0] = predictions[i].px_fast;
    xyz.reals[i * 3 + 1] = predictions[i].px_slow;
    xyz.reals[i * 3 + 2] = predictions[i].z;
    panel.ints[i] = predictions[i].panel;
    id.ints[i] = 0;
  }
  return t;
}

// The experiment without its crystal, as imported: what the search is given.
ExperimentList imported(const Experiment &e, const Vec3 &moved) {
  ExperimentList list;
  Experiment bare = e;
  bare.crystal.reset();
  for (Panel &p : bare.detector.panels)
    p.origin = p.origin + moved;
  list.experiments.push_back(bare);
  return list;
}

} // namespace

TEST(beam_search_from_the_true_centre_stays_there) {
  const Experiment truth = lattice_experiment();
  const Table spots = spots_of(truth);
  ExperimentList list = imported(truth, {0.0, 0.0, 0.0});
  const BeamSearchResult r = search_beam(list, spots, {});
  check::is_true(r.found, "found: " + r.why);
  check::is_true(r.candidates >= 3, "periodicities enough");
  check::close(r.offset.norm(), 0.0, 0.05, "moved by no more than 0.05 mm");
}

TEST(beam_search_moves_a_displaced_detector_back) {
  // 1.0 mm across and 0.7 mm down: 13 and 9 pixels.
  const Experiment truth = lattice_experiment();
  const Table spots = spots_of(truth);
  const Vec3 moved{1.0, -0.7, 0.0};
  ExperimentList list = imported(truth, moved);
  const BeamSearchResult r = search_beam(list, spots, {});
  check::is_true(r.found, "found: " + r.why);
  check::close(r.offset.x, -moved.x, 0.05, "back across, mm");
  check::close(r.offset.y, -moved.y, 0.05, "back down, mm");
  const Vec3 placed = list.experiments[0].detector.panels[0].origin;
  check::close(placed.x, truth.detector.panels[0].origin.x, 0.05,
               "the panel where it was");
}

} // namespace mxi
