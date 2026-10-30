// mxi verify: how much of what each lattice predicts was seen. Every
// reflection the lattice predicts over its scan, to the resolution its indexed
// spots reach, against the spots indexed to it, in shells of resolution --
// where indexing counts the spots a lattice explains, this counts what it
// should and does not (src/verify.hh). A prediction is observed where a strong
// spot of its sweep, indexed or not, has fractional indices within 0.3 of its
// own, as indexing judges one; the column "indexed" counts those also indexed
// to the lattice. After indexing, or any step after it:
// an integrated file's own rows are predictions, and its flags say which were
// seen. Predictions on masked pixels -- the detector's mask, the regions marked
// untrusted -- cannot be seen, and are counted apart where the images can be
// read (Graeme, 30 October).

#include <cmath>
#include <cstdio>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "args.hh"
#include "experiment_images.hh"
#include "expt.hh"
#include "image_series.hh"
#include "log_mirror.hh"
#include "refl.hh"
#include "rflx.hh"
#include "untrusted.hh"
#include "verify.hh"

namespace {

const char *const kUsage =
    "usage: mxi verify FILE [options]\n"
    "  -j, --threads N   threads for everything this does in parallel; 0, "
    "every\n"
    "                    core this job may use (0)\n"
    "How much of what each lattice predicts was seen: every reflection it\n"
    "predicts over its scan to the resolution its indexed spots reach, "
    "against\n"
    "the spots indexed to it, in shells. FILE is a .rflx -- indexed, refined,\n"
    "integrated or scaled -- or an .expt and its .refl.\n"
    "  --d-min D         to D A, not the indexed spots' own limit (the d 99.5\n"
    "                    per cent of them reach)\n"
    "  --shells N        shells of equal reciprocal volume (10)\n"
    "  --tolerance T     a prediction is observed where a strong spot's\n"
    "                    fractional indices are within T of its own, as\n"
    "                    indexing judges one (0.3, mxi index's)\n"
    "  --no-mask         count predictions on masked pixels as observable, "
    "the\n"
    "                    images not read\n";

// The detector's mask for experiment `index`: its own pixel mask, read with a
// frame, and the regions marked untrusted; empty where the images cannot be
// read, and `why` says so.
std::vector<std::uint8_t> mask_of(const std::string &path, std::size_t index,
                                  std::size_t *width, std::string *why) {
  try {
    const mxi::ExperimentImages info = mxi::experiment_images(path, index);
    if (info.image_file.empty()) {
      *why = "no image file named";
      return {};
    }
    std::unique_ptr<mxi::series::Series> series =
        mxi::series::nxmx(info.image_file);
    if (info.untrusted.is_array() && !info.untrusted.as_array().empty())
      series->add_masked(mxi::untrusted_pixels(info.untrusted, info.image_fast,
                                               info.image_slow));
    mxi::series::Info opened;
    if (!series->try_open(&opened)) {
      *why = info.image_file + " cannot be read";
      return {};
    }
    const std::vector<std::string> keys = series->ready();
    if (keys.empty()) {
      *why = info.image_file + " has no frames";
      return {};
    }
    mxi::series::Frame frame;
    if (!series->reader()->read(keys.front(), &frame)) {
      *why = info.image_file + ": the first frame cannot be read";
      return {};
    }
    std::vector<std::uint8_t> mask(opened.width * opened.height, 0);
    if (frame.masked)
      for (std::uint64_t at : *frame.masked)
        if (at < mask.size())
          mask[at] = 1;
    *width = opened.width;
    return mask;
  } catch (const std::exception &e) {
    *why = e.what();
    return {};
  }
}

} // namespace

// The program's entry, run by mxi (verify) -- apps/mxi.cc.
int mxi_verify_main(int argc, char **argv) {
  using namespace mxi;
  // Its table kept beside the data, as every analysis's is.
  if (!only_asks_for_help(argc, argv))
    mirror_to_log("mxi_verify.log");
  const std::set<std::string> known = {"--d-min",   "--shells", "--tolerance",
                                       "--no-mask", "--help",   "-h"};
  const std::set<std::string> takes_value = {"--d-min", "--shells",
                                             "--tolerance"};
  Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (!args.ok) {
    std::fprintf(stderr, "mxi verify: %s\n", args.error.c_str());
    std::fputs(kUsage, stderr);
    return 2;
  }
  if (args.has("--help") || args.has("-h")) {
    std::fputs(kUsage, stdout);
    return 0;
  }
  args.positional = rflx::as_pair(args.positional);
  if (args.positional.size() != 2) {
    std::fprintf(stderr,
                 "mxi verify: expected a .rflx, or an .expt and a .refl\n");
    std::fputs(kUsage, stderr);
    return 2;
  }
  VerifyOptions options;
  options.d_min = args.number("--d-min", 0.0);
  options.shells = static_cast<std::size_t>(args.number("--shells", 10.0));
  options.tolerance = args.number("--tolerance", 0.3);
  try {
    const ExperimentList experiments = read_experiments(args.positional[0]);
    const Table table = read_reflections(args.positional[1]);
    std::size_t lattices = 0;
    for (std::size_t i = 0; i < experiments.size(); ++i) {
      const Experiment &e = experiments[i];
      if (!e.crystal)
        continue;
      ++lattices;
      std::vector<std::uint8_t> mask;
      std::size_t width = 0;
      std::string why;
      if (!args.has("--no-mask"))
        mask = mask_of(args.positional[0], i, &width, &why);
      const VerifyResult r =
          verify_lattice(experiments, i, table, options,
                         mask.empty() ? nullptr : &mask, width);
      std::printf("\nExperiment %zu: %zu indexed spots", i, r.indexed);
      if (r.d_min <= 0.0) {
        std::printf(", so no limit to predict to: give --d-min\n");
        continue;
      }
      std::printf("; predicted to %.2f A (%s)\n", r.d_min,
                  r.d_min_from_spots ? "the indexed spots' limit" : "--d-min");
      if (!r.masked)
        std::printf("  masked pixels not taken into account: %s\n",
                    args.has("--no-mask") ? "--no-mask" : why.c_str());
      std::printf("  %-15s %10s %13s %9s %9s %9s\n", "d (A)", "predicted",
                  "unobservable", "observed", "fraction", "indexed");
      // A shell with nothing it could see has no fraction, not 0 per cent.
      const auto line = [](const char *label, const VerifyShell &s) {
        std::printf("  %-15s %10zu %13zu %9zu", label, s.predicted,
                    s.unobservable, s.observed);
        if (s.predicted > s.unobservable)
          std::printf(" %8.1f%%", 100.0 * s.fraction());
        else
          std::printf(" %9s", "-");
        std::printf(" %9zu\n", s.indexed);
      };
      for (const VerifyShell &s : r.shells) {
        char label[32];
        std::snprintf(label, sizeof label, "%6.2f - %-6.2f", s.d_max, s.d_min);
        line(label, s);
      }
      line("all", r.all);
    }
    if (lattices == 0) {
      std::fprintf(stderr, "mxi verify: %s has no lattice: index it first\n",
                   args.positional[0].c_str());
      return 1;
    }
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi verify: %s\n", e.what());
    return 1;
  }
  return 0;
}
