// mxi_untrusted: each image set's untrusted regions, listed, set from a mask
// file, or cleared -- the shapes mxi_view draws (src/untrusted.hh), written
// into the experiment list where every later program finds them. mxi_view
// calls it to save; on its own it applies a mask kept from before to another
// run's sweeps, with no window.

#include <cstdio>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "args.hh"
#include "json.hh"
#include "log_mirror.hh"
#include "processing.hh"
#include "refl.hh"
#include "rflx.hh"
#include "untrusted.hh"

namespace {

const char *const kUsage =
    "usage: mxi_untrusted FILE.rflx [--set MASK.json | --clear] [-o OUT.rflx]\n"
    "\n"
    "Each image set's untrusted regions -- circles, rectangles and polygons\n"
    "on the detector that spot finding and integration leave out, as they do\n"
    "the detector's own bad pixels. mxi_view draws them; this lists, sets and\n"
    "clears them. With neither --set nor --clear, lists them.\n"
    "\n"
    "  --set MASK.json  the regions from a mask file: {\"imagesets\": {\"0\":\n"
    "                   [shapes], \"1\": [shapes]}} for each image set named, "
    "or\n"
    "                   a list of shapes alone for every image set. A shape "
    "is\n"
    "                   {\"circle\": [x, y, r]}, {\"rectangle\": [x0, y0, x1,\n"
    "                   y1]} or {\"polygon\": [[x, y], ...]}, in pixels, panel "
    "0\n"
    "  --clear          every image set's regions removed\n"
    "  -o PATH          where to write (masked.rflx): the input is left as it\n"
    "                   is, to compare with; name it here to write over it\n"
    "  --help, --version\n";

std::string base_name(const std::string &path) {
  return std::filesystem::path(path).filename().string();
}

} // namespace

int main(int argc, char **argv) {
  using namespace mxi;
  if (!only_asks_for_help(argc, argv))
    mirror_to_log("mxi_untrusted.log");
  const Arguments args = parse_arguments(
      argc, argv, {"-o", "--set", "--clear", "--help"}, {"-o", "--set"});
  if (args.help) {
    std::printf("%s", kUsage);
    return 0;
  }
  if (!args.ok || args.positional.size() != 1 ||
      (args.has("--set") && args.has("--clear"))) {
    std::fprintf(stderr, "mxi_untrusted: %s\n%s",
                 !args.ok                      ? args.error.c_str()
                 : args.positional.size() != 1 ? "one experiment list"
                                               : "--set or --clear, not both",
                 kUsage);
    return 2;
  }
  try {
    const std::string &in = args.positional[0];
    std::optional<json::Value> document = rflx::read_experiments(in);
    if (!document)
      throw std::runtime_error(in + ": a .rflx holding an experiment list");
    json::Value &imagesets = document->as_object()["imageset"];
    if (!imagesets.is_array() || imagesets.as_array().empty())
      throw std::runtime_error(in + " has no image sets");
    json::Array &sets = imagesets.as_array();

    const bool changing = args.has("--set") || args.has("--clear");
    if (args.has("--set")) {
      const json::Value mask = json::parse_file(args.value("--set", ""));
      const auto apply = [&](std::size_t k, const json::Value &shapes) {
        if (k >= sets.size())
          throw std::runtime_error("the mask names image set " +
                                   std::to_string(k) + "; the file has " +
                                   std::to_string(sets.size()));
        (void)untrusted_pixels(shapes, 1, 1); // each shape's form, checked
        sets[k].as_object()["untrusted"] = shapes;
      };
      if (mask.is_array()) {
        for (std::size_t k = 0; k < sets.size(); ++k)
          apply(k, mask);
      } else if (mask["imagesets"].is_object()) {
        for (const auto &[key, shapes] : mask["imagesets"].as_object())
          apply(static_cast<std::size_t>(std::stoul(key)), shapes);
      } else {
        throw std::runtime_error(
            "a mask file is a list of shapes, or {\"imagesets\": {...}}");
      }
    } else if (args.has("--clear")) {
      for (json::Value &s : sets)
        if (s.is_object())
          s.as_object().erase("untrusted");
    }

    for (std::size_t k = 0; k < sets.size(); ++k) {
      const json::Value &s = sets[k];
      std::printf("image set %zu, %s: %s\n", k,
                  s["template"].is_string()
                      ? base_name(s["template"].as_string()).c_str()
                      : "?",
                  describe_untrusted(s["untrusted"]).c_str());
    }
    if (!changing)
      return 0;

    // A new file, the input left as it was, for comparison (Graeme, 22
    // October); the input itself only when named.
    const std::string out = args.value("-o", "masked.rflx");
    std::vector<std::string> skipped;
    const std::optional<Table> table = rflx::read_reflections(in, &skipped);
    Summary step("mxi_untrusted", in);
    step.set("imagesets", json::Value(static_cast<long>(sets.size())));
    const json::Value steps = step.document();
    const std::string temporary = out + ".mxi_untrusted.tmp";
    rflx::write(temporary, &*document, table ? &*table : nullptr,
                "mxi_untrusted", false, &steps);
    std::filesystem::rename(temporary, out);
    std::printf("Wrote %s\n", out.c_str());
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi_untrusted: %s\n", e.what());
    return 1;
  }
}
