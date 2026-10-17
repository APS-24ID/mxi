// The spot finder's spots as a reflection table: see spot_table.hh.

#include "spot_table.hh"

#include <bit>
#include <cstdint>
#include <cstring>
#include <stdexcept>

namespace mxi {

namespace {

// DIALS's flags and shoebox mask bits, as dials.find_spots sets them.
constexpr std::int64_t kStrong = 1 << 5;
constexpr std::uint8_t kValid = 1 << 0;
constexpr std::uint8_t kForeground = 1 << 2;

template <typename T> void append(std::string &out, T value) {
  char bytes[sizeof(T)];
  std::memcpy(bytes, &value, sizeof(T));
  out.append(bytes, sizeof(T));
}

} // namespace

Table spot_table(const std::vector<SpotPart> &parts, std::size_t panel,
                 bool shoeboxes) {
  if constexpr (std::endian::native != std::endian::little)
    throw std::runtime_error(
        "spot_table: the shoeboxes are raw little-endian bytes, and this "
        "machine is not little-endian");
  std::size_t rows = 0;
  for (const SpotPart &part : parts) {
    if (part.width == 0)
      throw std::runtime_error("refl: a frame width of zero");
    rows += part.spots->size();
  }
  Table table;
  table.nrows = rows;
  table.version = 2;
  for (const SpotPart &part : parts)
    if (!part.identifier.empty())
      table.identifiers[static_cast<std::size_t>(part.id)] = part.identifier;

  Column &bbox = table.int_column("bbox", "int6", 6);
  Column &flags = table.int_column("flags", "std::size_t", 1);
  Column &id = table.int_column("id", "int", 1);
  Column &value = table.real_column("intensity.sum.value", "double", 1);
  Column &variance = table.real_column("intensity.sum.variance", "double", 1);
  Column &n_signal = table.int_column("n_signal", "int", 1);
  Column &panels = table.int_column("panel", "std::size_t", 1);
  Column &xyz = table.real_column("xyzobs.px.value", "vec3<double>", 3);
  Column &xyz_variance =
      table.real_column("xyzobs.px.variance", "vec3<double>", 3);
  // The columns come sized to the table's rows (Table::int_column): filled
  // row by row.
  std::size_t r = 0;
  for (const SpotPart &part : parts)
    for (const dials_spots::Spot &spot : *part.spots) {
      for (std::size_t i = 0; i < 6; ++i)
        bbox.ints[r * 6 + i] = spot.bbox[i];
      flags.ints[r] = kStrong;
      id.ints[r] = part.id;
      value.reals[r] = spot.intensity;
      variance.reals[r] = spot.intensity_variance;
      n_signal.ints[r] = static_cast<std::int32_t>(spot.n_signal);
      panels.ints[r] = static_cast<std::int64_t>(panel);
      for (std::size_t axis = 0; axis < 3; ++axis) {
        xyz.reals[r * 3 + axis] = spot.position[axis];
        xyz_variance.reals[r * 3 + axis] = spot.variance[axis];
      }
      ++r;
    }

  if (shoeboxes) {
    // Each spot's shoebox as DIALS's Shoebox<> serialises it: the panel, the
    // bbox, the version (2: a uint8 mask), then the box's data, mask and
    // background, c_grid<3>(z, y, x) -- x fastest. Signal pixels only, each
    // valid and foreground; the background each one's local background under
    // --subtract-background, zero otherwise, as DIALS's spot finder writes it.
    Table::Opaque boxes;
    boxes.type = "Shoebox<>";
    boxes.rows = rows;
    std::vector<float> data, background;
    std::vector<std::uint8_t> mask;
    for (const SpotPart &part : parts)
      for (const dials_spots::Spot &spot : *part.spots) {
        const auto xsize =
            static_cast<std::size_t>(spot.bbox[1] - spot.bbox[0]);
        const auto ysize =
            static_cast<std::size_t>(spot.bbox[3] - spot.bbox[2]);
        const auto zsize =
            static_cast<std::size_t>(spot.bbox[5] - spot.bbox[4]);
        const std::size_t n = xsize * ysize * zsize;
        data.assign(n, 0.0f);
        background.assign(n, 0.0f);
        mask.assign(n, 0);
        for (std::uint32_t k = 0; k < spot.n_signal; ++k) {
          const dials_spots::Pixel &pixel = (*part.pixels)[spot.first + k];
          const std::size_t x =
              static_cast<std::size_t>(pixel.index % part.width) -
              static_cast<std::size_t>(spot.bbox[0]);
          const std::size_t y =
              static_cast<std::size_t>(pixel.index / part.width) -
              static_cast<std::size_t>(spot.bbox[2]);
          const std::size_t z = static_cast<std::size_t>(
              static_cast<std::int64_t>(pixel.frame) - spot.bbox[4]);
          const std::size_t at = (z * ysize + y) * xsize + x;
          data[at] = static_cast<float>(pixel.value);
          background[at] = pixel.background;
          mask[at] = kValid | kForeground;
        }
        append<std::uint32_t>(boxes.bytes, static_cast<std::uint32_t>(panel));
        for (int i = 0; i < 6; ++i)
          append<std::int32_t>(boxes.bytes, spot.bbox[i]);
        boxes.bytes.push_back(static_cast<char>(2));
        boxes.bytes.append(reinterpret_cast<const char *>(data.data()),
                           n * sizeof(float));
        boxes.bytes.append(reinterpret_cast<const char *>(mask.data()), n);
        boxes.bytes.append(reinterpret_cast<const char *>(background.data()),
                           n * sizeof(float));
      }
    table.set_opaque("shoebox", std::move(boxes));
  }
  return table;
}

void write_spots(const std::string &path,
                 const std::vector<dials_spots::Spot> &spots,
                 const std::vector<dials_spots::Pixel> &pixels,
                 std::size_t width, const SpotTableOptions &options) {
  SpotPart part;
  part.spots = &spots;
  part.pixels = &pixels;
  part.width = width;
  part.id = options.id;
  part.identifier = options.identifier;
  write_reflections(path, spot_table({part}, options.panel, options.shoeboxes));
}

} // namespace mxi
