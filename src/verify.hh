// How much of what a lattice predicts was seen: every reflection the lattice
// predicts over its scan, to the resolution its indexed spots reach, against
// the spots indexed to it. Indexing counts the spots a lattice explains; this
// counts the reflections it should explain and does not -- a wrong cell, a
// missed twin or a second crystal can leave the first high and this low
// (Graeme, 30 October). mxi verify reports it, after indexing or any later
// step.

#ifndef MXI_VERIFY_HH
#define MXI_VERIFY_HH

#include <cstddef>
#include <cstdint>
#include <vector>

#include "expt.hh"
#include "geometry.hh"
#include "refl.hh"

namespace mxi {

struct VerifyOptions {
  //: To this resolution, in Å; 0 for the indexed spots' own limit: the d that
  //: 99.5 per cent of them reach, so that a few misindexed spots far out do not
  //: carry it.
  double d_min = 0.0;
  //: Shells, of equal volume in reciprocal space.
  std::size_t shells = 10;
  //: A prediction is seen where a strong spot of its sweep -- indexed or not,
  //: to any lattice -- has fractional indices, by this lattice's matrix,
  //: within this of its own, as indexing judges a spot (mxi index --tolerance,
  //: 0.3); each spot taken by the nearest in rotation of its reflection's
  //: crossings. Judged by indexing's assignment instead, as at first, a spot
  //: in place but outside its tolerance was missing, and the outer shells
  //: fell for it (Graeme, 30 October).
  double tolerance = 0.3;
};

struct VerifyShell {
  double d_max = 0.0, d_min = 0.0; //: Å
  std::size_t predicted = 0;       //: every reflection predicted
  std::size_t unobservable = 0;    //: its centre on a masked pixel
  std::size_t observed = 0;        //: of the rest, seen
  std::size_t indexed = 0;         //: of those, by a spot indexed to it
  double fraction() const {
    const std::size_t could = predicted - unobservable;
    return could > 0
               ? static_cast<double>(observed) / static_cast<double>(could)
               : 0.0;
  }
};

struct VerifyResult {
  double d_min = 0.0;           //: the limit used, Å
  bool d_min_from_spots = true; //: the indexed spots' own, not given
  std::size_t indexed = 0;      //: the lattice's indexed spots
  std::size_t spots = 0;        //: the sweep's strong spots looked at
  bool masked = false;          //: a mask was given
  std::vector<VerifyShell> shells;
  VerifyShell all;
};

//: Lattice `id` of `experiments` (its crystal, scan and detector) against the
//: strong spots of its sweep in `table`: rows flagged strong or indexed with a
//: measured position, whatever their experiment. `mask`, where not null, is
//: the detector's: a byte a pixel, row by row, non-zero where masked;
//: predictions on it are unobservable.
VerifyResult verify_lattice(const ExperimentList &experiments, std::size_t id,
                            const Table &table, const VerifyOptions &options,
                            const std::vector<std::uint8_t> *mask = nullptr,
                            std::size_t mask_width = 0);

} // namespace mxi

#endif // MXI_VERIFY_HH
