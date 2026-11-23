// The beam centre found from the strong spots, as dials.search_beam_position's
// default method -- LABELIT's (Sauter, Grosse-Kunstleve & Adams 2004, J. Appl.
// Cryst. 37, 399-409; DIALS's labelit_method.py).
//
// The idea: with the beam centre right, the reciprocal lattice's planes pass
// through the origin, so along each of the lattice's strong periodicities t
// every spot's reciprocal lattice point r projects to an integer, r.t, and
// the phase of sum exp(2 pi i r.t) is zero. A wrong centre displaces the
// lattice from the origin and the phases move off zero. So: the strongest
// periodicities found once (the three-dimensional FFT's candidate vectors,
// as indexing finds them -- DIALS uses LABELIT's one-dimensional search for
// the same thing); each trial centre, the detector moved across the beam,
// scored by the sum over them of cos(phase), at most one each; the best
// found by a grid over the search's scope, then refined.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#include "expt.hh"
#include "geometry.hh"
#include "linalg.hh"
#include "refl.hh"

namespace mxi {

struct BeamSearchOptions {
  double scope_mm = 4.0;   //: the grid's reach either way across the beam
  double binning = 2.0;    //: its step, in pixels
  double d_min = 0.0;      //: spots finer than this left out (0: none)
  double max_cell = 0.0;   //: the largest cell, A (0: from the spots)
  std::size_t max_reflections = 10000; //: a random subset beyond, as DIALS
  std::uint32_t seed = 42;
  std::size_t candidates = 20; //: the strongest periodicities scored, as DIALS
};

struct BeamSearchResult {
  bool found = false;
  std::string why; //: when not found
  Vec3 offset;     //: the detector's shift, lab mm, across the beam
  double score_start = 0.0, score_best = 0.0;
  std::size_t candidates = 0, spots = 0, trials = 0;
  std::size_t peaks = 0; //: the grid's local maxima refined
  double max_cell = 0.0, d_min = 0.0;
  //: Where the beam meets the detector, before and after: panel, mm, px.
  std::size_t panel_before = 0, panel_after = 0;
  double mm_before[2] = {0, 0}, mm_after[2] = {0, 0};
  double px_before[2] = {0, 0}, px_after[2] = {0, 0};
};

//: Where the first experiment's beam meets its detector: the panel, and the
//: point in mm and in pixels (no parallax). False if it misses every panel.
bool beam_centre(const Experiment &e, std::size_t *panel, double mm[2],
                 double px[2]);

//: Every experiment's detector moved by `offset` (lab mm).
void move_detectors(ExperimentList &experiments, const Vec3 &offset);

//: The search. On success the experiments' detectors are moved to the centre
//: found; otherwise left as they were, `why` saying what failed.
BeamSearchResult search_beam(ExperimentList &experiments, const Table &strong,
                             const BeamSearchOptions &options);

} // namespace mxi
