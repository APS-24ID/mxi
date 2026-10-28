// mxi_split: a file's data sets apart, each into a file of its own -- or the
// ones named into one -- to be processed or judged alone.
//
//   mxi_split scaled.rflx                       # split_0.rflx, split_1.rflx,
//   ... mxi_split scaled.rflx --experiments 0,3 -o two.rflx
//
// Each data set goes with its crystal-less sweep experiment, if it has one,
// and its reflections renumbered (src/combine.hh), as mxi_correlation's
// clusters are written.

#include <cstdio>
#include <cstdlib>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "args.hh"
#include "combine.hh"
#include "log_mirror.hh"
#include "processing.hh"
#include "refl.hh"
#include "rflx.hh"

namespace {

const char *const kUsage =
    "usage: mxi_split FILE.rflx [--experiments 0,3,...] [-o PATH | --prefix "
    "P]\n"
    "\n"
    "A file's data sets apart: each data set -- each experiment with a "
    "crystal,\n"
    "with its crystal-less sweep experiment if it has one -- into a file of "
    "its\n"
    "own, P_0.rflx, P_1.rflx, ... numbered by experiment; or, with\n"
    "--experiments, the ones named together into one. Before indexing, with "
    "no\n"
    "crystals, each sweep. mxi_combine does the reverse.\n"
    "\n"
    "  --experiments I,J  those experiments, by index, into one file\n"
    "  -o PATH            where that one file goes (split.rflx)\n"
    "  --prefix P         the files' names, each apart (split)\n"
    "  --help, --version\n";

} // namespace

// The program's entry, run by mxi (split) -- apps/mxi.cc.
int mxi_split_main(int argc, char **argv) {
  using namespace mxi;
  if (!only_asks_for_help(argc, argv))
    mirror_to_log("mxi_split.log");
  const Arguments args =
      parse_arguments(argc, argv, {"-o", "--prefix", "--experiments", "--help"},
                      {"-o", "--prefix", "--experiments"});
  if (args.help) {
    std::printf("%s", kUsage);
    return 0;
  }
  if (!args.ok || args.positional.size() != 1) {
    std::fprintf(stderr, "mxi_split: %s\n%s",
                 args.ok ? "one file to split" : args.error.c_str(), kUsage);
    return 2;
  }
  try {
    const std::string &in = args.positional[0];
    const std::optional<json::Value> document = rflx::read_experiments(in);
    std::vector<std::string> skipped;
    const std::optional<Table> table = rflx::read_reflections(in, &skipped);
    if (!document || !table)
      throw std::runtime_error(in + ": a .rflx with experiments and "
                                    "reflections, to split");
    const std::optional<json::Value> summary = rflx::read_summary(in);
    const json::Array &all = (*document)["experiment"].as_array();
    bool any_crystal = false;
    for (const json::Value &e : all)
      any_crystal = any_crystal || e["crystal"].is_number();

    auto write = [&](const std::vector<bool> &keep, const std::string &path) {
      const ExperimentsAndReflections part =
          select_experiments(*document, *table, keep);
      Summary step("mxi_split", in);
      step.set("experiments",
               json::Value(static_cast<long>(
                   part.document["experiment"].as_array().size())));
      const json::Value steps = step.document();
      rflx::write(path, &part.document, &part.reflections, "mxi_split", false,
                  &steps);
      std::printf("Wrote %s: %zu experiments, %zu reflections\n", path.c_str(),
                  part.document["experiment"].as_array().size(),
                  part.reflections.nrows);
    };

    if (args.has("--experiments")) {
      std::vector<bool> keep(all.size(), false);
      std::stringstream list(args.value("--experiments", ""));
      std::string item;
      while (std::getline(list, item, ',')) {
        char *end = nullptr;
        const long i = std::strtol(item.c_str(), &end, 10);
        if (item.empty() || *end != '\0' || i < 0 ||
            static_cast<std::size_t>(i) >= all.size())
          throw std::runtime_error("--experiments: '" + item +
                                   "' is not an experiment of " + in +
                                   ", which has " + std::to_string(all.size()));
        keep[static_cast<std::size_t>(i)] = true;
      }
      write(keep, args.value("-o", "split.rflx"));
      return 0;
    }
    const std::string prefix = args.value("--prefix", "split");
    for (std::size_t e = 0; e < all.size(); ++e) {
      if (any_crystal && !all[e]["crystal"].is_number())
        continue; // a sweep goes with its data sets
      std::vector<bool> keep(all.size(), false);
      keep[e] = true;
      write(keep, prefix + "_" + std::to_string(e) + ".rflx");
    }
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi_split: %s\n", e.what());
    return 1;
  }
}
