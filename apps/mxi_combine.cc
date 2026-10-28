// mxi_combine: several runs' files joined into one, to be taken on together --
// data sets processed separately, then symmetry and scaling over them all.
//
//   mxi_combine a/integrated.rflx b/integrated.rflx -o combined.rflx
//   mxi_symmetry combined.rflx && mxi_scale symmetrized.rflx
//
// Each input's models are appended and every reference to them offset, its
// reflections renumbered to their experiments in the whole (src/combine.hh):
// as dials.combine_experiments does.

#include <cstdio>
#include <optional>
#include <set>
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
    "usage: mxi combine A.rflx B.rflx [...] [-o combined.rflx]\n"
    "\n"
    "Several runs' files joined into one: every experiment and its models, "
    "and\n"
    "every reflection, renumbered to its experiment in the whole. Typically\n"
    "integrated data sets, each processed on its own, to go through\n"
    "mxi symmetry and mxi scale together; mxi split does the reverse.\n"
    "\n"
    "The same experiment twice (one identifier in two inputs) is refused. The\n"
    "columns kept are those every input has, shoeboxes included where every\n"
    "input has them; any other is left out, and named. A DIALS pair is made a\n"
    ".rflx first by mxi convert.\n"
    "\n"
    "  -o PATH           where to write (combined.rflx)\n"
    "  --help, --version\n";

} // namespace

// The program's entry, run by mxi (combine) -- apps/mxi.cc.
int mxi_combine_main(int argc, char **argv) {
  using namespace mxi;
  if (!only_asks_for_help(argc, argv))
    mirror_to_log("mxi_combine.log");
  const Arguments args = parse_arguments(argc, argv, {"-o", "--help"}, {"-o"});
  if (args.help) {
    std::printf("%s", kUsage);
    return 0;
  }
  if (!args.ok || args.positional.size() < 2) {
    std::fprintf(stderr, "mxi combine: %s\n%s",
                 args.ok ? "two files or more to join" : args.error.c_str(),
                 kUsage);
    return 2;
  }
  try {
    std::vector<json::Value> documents;
    std::vector<Table> tables;
    json::Array inputs;
    for (const std::string &path : args.positional) {
      if (!rflx::has_experiments(path))
        throw std::runtime_error(path +
                                 " is not a .rflx with an experiment "
                                 "list (a DIALS pair: mxi convert first)");
      const std::optional<json::Value> experiments =
          rflx::read_experiments(path);
      std::vector<std::string> skipped;
      std::optional<Table> table = rflx::read_reflections(path, &skipped);
      if (!experiments || !table)
        throw std::runtime_error(path + " has no reflections to join");
      const std::size_t n = (*experiments)["experiment"].as_array().size();
      std::printf("%s: %zu experiments, %zu reflections\n", path.c_str(), n,
                  table->nrows);
      json::Object input;
      input["path"] = json::Value(path);
      input["experiments"] = json::Value(static_cast<long>(n));
      input["reflections"] = json::Value(static_cast<long>(table->nrows));
      if (const auto s = rflx::read_summary(path)) {
        const std::vector<json::Value> steps = Summary::steps(*s);
        input["steps"] = json::Value(json::Array(steps.begin(), steps.end()));
      }
      inputs.push_back(json::Value(std::move(input)));
      documents.push_back(*experiments);
      tables.push_back(std::move(*table));
    }
    const Combined combined = combine_experiments(documents, tables);
    for (const std::string &name : combined.dropped)
      std::printf("Left out the column %s: not in every input alike\n",
                  name.c_str());
    // The summary starts again: one step, naming each input with its own steps.
    Summary summary("mxi_combine", "");
    summary.set("inputs", json::Value(std::move(inputs)));
    const json::Value document = summary.document();
    const std::string out = args.value("-o", "combined.rflx");
    rflx::write(out, &combined.document, &combined.reflections, "mxi_combine",
                false, &document);
    std::printf("Wrote %s: %zu experiments, %zu reflections\n", out.c_str(),
                combined.document["experiment"].as_array().size(),
                combined.reflections.nrows);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi combine: %s\n", e.what());
    return 1;
  }
}
