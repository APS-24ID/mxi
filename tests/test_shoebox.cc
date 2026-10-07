// The shoebox column, which this package skipped for a long time while saying
// it had no pixels. It has 13.6 MB of them for insulin.

#include <cstdio>
#include <string>
#include <vector>

#include "../src/shoebox.hh"
#include "check.hh"

namespace mxi {

namespace {

Shoebox made(std::int32_t x0, std::int32_t nx, std::int32_t ny,
             std::int32_t nz) {
  Shoebox box;
  box.panel = 0;
  box.bbox[0] = x0;
  box.bbox[1] = x0 + nx;
  box.bbox[2] = 100;
  box.bbox[3] = 100 + ny;
  box.bbox[4] = 0;
  box.bbox[5] = nz;
  box.flag = 2;
  const std::size_t n = box.size();
  for (std::size_t i = 0; i < n; ++i) {
    box.data.push_back(static_cast<float>(i % 13));
    box.mask.push_back(i % 3 == 0
                           ? shoebox_mask::kValid | shoebox_mask::kForeground
                           : shoebox_mask::kValid);
    box.background.push_back(0.0f);
  }
  return box;
}

Table table_with(const std::vector<Shoebox> &boxes) {
  Table t;
  t.nrows = boxes.size();
  Table::Opaque column;
  column.type = "Shoebox<>";
  column.bytes = encode_shoeboxes(boxes);
  column.rows = boxes.size();
  t.set_opaque("shoebox", std::move(column));
  return t;
}

} // namespace

TEST(shoeboxes_survive_an_encode_and_a_decode) {
  const std::vector<Shoebox> boxes = {made(2010, 4, 4, 2), made(1913, 2, 2, 1),
                                      made(1775, 4, 5, 2)};
  const std::vector<Shoebox> back = decode_shoeboxes(table_with(boxes));

  check::equal(static_cast<long long>(back.size()), 3, "all three");
  for (std::size_t i = 0; i < boxes.size(); ++i) {
    check::equal(static_cast<long long>(back[i].panel), boxes[i].panel,
                 "panel");
    for (int k = 0; k < 6; ++k) {
      check::equal(static_cast<long long>(back[i].bbox[k]),
                   static_cast<long long>(boxes[i].bbox[k]), "bbox");
    }
    check::equal(static_cast<long long>(back[i].flag), boxes[i].flag, "flag");
    check::is_true(back[i].data == boxes[i].data, "data");
    check::is_true(back[i].mask == boxes[i].mask, "mask");
    check::is_true(back[i].background == boxes[i].background, "background");
  }
}

TEST(the_voxel_order_is_z_then_y_then_x) {
  // Slowest to fastest, which is the order the bytes are in and the order the
  // bounding box names its axes. Getting this backwards would transpose every
  // profile without changing any total, so no test of a summed intensity would
  // ever see it.
  Shoebox box = made(0, 4, 3, 2);
  for (std::size_t i = 0; i < box.size(); ++i)
    box.data[i] = static_cast<float>(i);
  check::close(box.data[box.at(0, 0, 0)], 0.0, 0.0, "first voxel");
  check::close(box.data[box.at(1, 0, 0)], 1.0, 0.0, "x is fastest");
  check::close(box.data[box.at(0, 1, 0)], 4.0, 0.0, "then y, by nx");
  check::close(box.data[box.at(0, 0, 1)], 12.0, 0.0, "then z, by nx * ny");
  check::equal(static_cast<long long>(box.size()), 24,
               "and the size is nx ny nz");
}

TEST(a_blob_that_does_not_parse_exactly_is_refused) {
  // Bytes left over mean the layout is wrong somewhere earlier, and the boxes
  // already read are wrong with it. Reporting them would be worse than
  // refusing: a profile model fitted to misparsed pixels would look plausible.
  const std::vector<Shoebox> boxes = {made(0, 3, 3, 1)};
  Table t = table_with(boxes);
  Table::Opaque longer = t.opaque().at("shoebox");
  longer.bytes += std::string(8, '\0');
  t.set_opaque("shoebox", longer);

  bool threw = false;
  try {
    decode_shoeboxes(t);
  } catch (const ReflError &e) {
    threw = true;
    const std::string what = e.what();
    check::is_true(what.find("left over") != std::string::npos,
                   "says what is wrong");
  }
  check::is_true(threw, "must refuse trailing bytes");

  // And a record cut short.
  Table shorter = table_with(boxes);
  Table::Opaque cut = shorter.opaque().at("shoebox");
  cut.bytes.resize(cut.bytes.size() - 4);
  shorter.set_opaque("shoebox", cut);
  threw = false;
  try {
    decode_shoeboxes(shorter);
  } catch (const ReflError &) {
    threw = true;
  }
  check::is_true(threw, "must refuse a truncated record");
}

TEST(a_table_without_shoeboxes_decodes_to_nothing) {
  Table t;
  t.nrows = 5;
  check::is_true(decode_shoeboxes(t).empty(), "not an error, just empty");
}

} // namespace mxi

namespace mxi {

TEST(selecting_rows_keeps_each_rows_own_shoebox) {
  // Several sweeps' strong spots divided among them: rows 2 and 0, in that
  // order, with their own shoeboxes -- told apart by width and position.
  Table t =
      table_with({made(10, 3, 2, 1), made(20, 4, 2, 2), made(30, 5, 3, 1)});
  Column &x = t.real_column("intensity.sum.value", "double", 1);
  x.reals = {1.0, 2.0, 3.0};
  const Table s = select_rows_with_shoeboxes(t, {2, 0});
  check::equal(static_cast<long long>(s.nrows), 2LL, "two rows");
  check::close(s.at("intensity.sum.value").real(0), 3.0, 0.0, "row 2 first");
  check::close(s.at("intensity.sum.value").real(1), 1.0, 0.0, "row 0 second");
  const std::vector<Shoebox> boxes = decode_shoeboxes(s);
  check::equal(static_cast<long long>(boxes.size()), 2LL, "two shoeboxes");
  check::equal(static_cast<long long>(boxes[0].bbox[0]), 30LL, "row 2's box");
  check::equal(static_cast<long long>(boxes[0].nx()), 5LL, "row 2's width");
  check::equal(static_cast<long long>(boxes[1].bbox[0]), 10LL, "row 0's box");
  check::equal(static_cast<long long>(boxes[1].nx()), 3LL, "row 0's width");
}

} // namespace mxi

namespace mxi {

TEST(tables_join_with_their_shoeboxes_each_rows_own) {
  // Several sweeps' integrated tables under --save-shoeboxes: rows and
  // shoeboxes in order, each row its own; a table without joined to one with,
  // refused.
  Table a = table_with({made(10, 3, 2, 1), made(20, 4, 2, 2)});
  a.real_column("intensity.sum.value", "double", 1).reals = {1.0, 2.0};
  Table b = table_with({made(30, 5, 3, 1)});
  b.real_column("intensity.sum.value", "double", 1).reals = {3.0};
  const Table joined = concat_rows_with_shoeboxes({a, b});
  check::equal(static_cast<long long>(joined.nrows), 3LL, "three rows");
  const std::vector<Shoebox> boxes = decode_shoeboxes(joined);
  check::equal(static_cast<long long>(boxes.size()), 3LL, "three shoeboxes");
  const long long first[3] = {10, 20, 30}, width[3] = {3, 4, 5};
  for (std::size_t i = 0; i < 3; ++i) {
    check::equal(static_cast<long long>(boxes[i].bbox[0]), first[i],
                 "each its own box");
    check::equal(static_cast<long long>(boxes[i].nx()), width[i],
                 "each its own width");
    check::close(joined.at("intensity.sum.value").real(i),
                 static_cast<double>(i + 1), 0.0, "and its own row");
  }
  Table bare = b;
  bare.remove_opaque("shoebox");
  bool refused = false;
  try {
    concat_rows_with_shoeboxes({a, bare});
  } catch (const ReflError &) {
    refused = true;
  }
  check::is_true(refused, "with and without, refused");
}

} // namespace mxi
