// The spot finder's spots as a reflection table: what mxi_find writes, into a
// .rflx or a DIALS .refl, through the library's writers. The columns are those
// dials.find_spots writes -- bbox, flags (strong), id, intensity.sum.value and
// variance, n_signal, panel, the shoeboxes (optionally), xyzobs.px.value and
// variance -- and the library writes them in DIALS's own order, so that a
// .refl from here is byte for byte what DIALS writes for the same spots.

#ifndef MXI_SPOT_TABLE_HH
#define MXI_SPOT_TABLE_HH

#include <cstddef>
#include <string>
#include <vector>

#include "refl.hh"
#include "spots.hh"

namespace mxi {

//: One experiment's spots, from the spot finder: the spots, their pixels, the
//: image's width, and the experiment's id and identifier in the table.
struct SpotPart {
  const std::vector<mxi::spots::Spot> *spots = nullptr;
  const std::vector<mxi::spots::Pixel> *pixels = nullptr;
  std::size_t width = 0;
  int id = 0;
  std::string identifier;
};

//: Every part's spots as one table, a part's rows together, in order; the
//: shoeboxes as DIALS's Shoebox<> column (version 2, a uint8 mask) or left
//: out. One panel: a segmented detector is not handled. Throws on a part of
//: width zero, or on a big-endian machine, where the shoeboxes' raw bytes
//: would not read back.
Table spot_table(const std::vector<SpotPart> &parts, std::size_t panel,
                 bool shoeboxes);

//: One experiment's spots written to `path` (.refl), for the tests.
struct SpotTableOptions {
  std::size_t panel = 0;
  int id = 0;
  std::string identifier;
  bool shoeboxes = true;
};
void write_spots(const std::string &path,
                 const std::vector<mxi::spots::Spot> &spots,
                 const std::vector<mxi::spots::Pixel> &pixels,
                 std::size_t width, const SpotTableOptions &options);

} // namespace mxi

#endif // MXI_SPOT_TABLE_HH
