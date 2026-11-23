// mxi search-beam: the beam centre found from the strong spots, as
// dials.search_beam_position's default method, LABELIT's (src/beam_search.hh,
// docs/search_beam.md). Where the header's
// centre is too far out for indexing to find the lattice, which then refines
// it no further.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <stdexcept>
#include <string>

#include "args.hh"
#include "beam_search.hh"
#include "expt.hh"
#include "log_mirror.hh"
#include "output.hh"
#include "refl.hh"
#include "rflx.hh"

namespace mxi {
namespace {

const char *const kUsage =
    "usage: mxi search-beam STRONG.rflx [options]\n"
    "       mxi search-beam IMPORTED.expt STRONG.refl [options]\n"
    "The beam centre found from the strong spots, as dials.search_beam_position\n"
    "finds it by default (LABELIT's method): with the centre right, the\n"
    "reciprocal lattice's planes pass through the origin, so along each of the\n"
    "lattice's strong periodicities every spot projects to an integer. The\n"
    "detector moved across the beam over a grid, each trial scored by how\n"
    "nearly that holds, the best refined. For when the header's centre is too\n"
    "far out for indexing to find the lattice.\n"
    "  --scope MM        the search's reach either way across the beam (4)\n"
    "  --binning N       its grid's step, in pixels (2)\n"
    "  --d-min D         spots finer than D A left out\n"
    "  --max-cell A      the largest cell (from the spots)\n"
    "  --max-reflections N   a random subset of at most N spots (10000)\n"
    "  -o, --output F    optimised.rflx; .expt and .refl for the pair\n";

int run_search_beam(int argc, char **argv) {
  const std::set<std::string> known = {
      "--help",     "-h",          "--scope",       "--binning",
      "--d-min",    "--max-cell",  "--max-reflections", "-o",
      "--output",   "--output-expt", "--output-refl"};
  const std::set<std::string> takes_value = {
      "--scope", "--binning",  "--d-min",       "--max-cell",
      "--max-reflections", "-o", "--output", "--output-expt",
      "--output-refl"};
  Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (!args.ok) {
    mxi::eprint("mxi search-beam: %s\n", args.error.c_str());
    mxi::put(mxi::Stream::err, kUsage);
    return 2;
  }
  if (args.has("--help") || args.has("-h")) {
    mxi::put(mxi::Stream::out, kUsage);
    return 0;
  }
  args.positional = rflx::as_pair(args.positional);
  if (args.positional.size() != 2) {
    mxi::put(mxi::Stream::err, kUsage);
    return 2;
  }
  if (!only_asks_for_help(argc, argv))
    mirror_to_log("mxi_search_beam.log");
  try {
    ExperimentList experiments = read_experiments(args.positional[0]);
    const Table strong = read_reflections(args.positional[1]);
    BeamSearchOptions options;
    options.scope_mm = args.number("--scope", options.scope_mm);
    options.binning = args.number("--binning", options.binning);
    options.d_min = args.number("--d-min", 0.0);
    options.max_cell = args.number("--max-cell", 0.0);
    options.max_reflections = static_cast<std::size_t>(std::max(
        50.0, args.number("--max-reflections", 10000.0)));
    if (!(options.scope_mm > 0.0) || !(options.binning > 0.0)) {
      mxi::eprint("mxi search-beam: --scope and --binning are above 0\n");
      return 2;
    }
    const BeamSearchResult r = search_beam(experiments, strong, options);
    mxi::print("Beam centre as given: %.2f, %.2f px (%.2f, %.2f mm)\n",
               r.px_before[0], r.px_before[1], r.mm_before[0], r.mm_before[1]);
    if (r.spots > 0)
      mxi::print("  %zu spots; %zu periodicities, the strongest (the largest "
                 "cell %.1f A, to %.2f A)\n",
                 r.spots, r.candidates, r.max_cell, r.d_min);
    if (!r.found) {
      mxi::eprint("mxi search-beam: no centre found: %s\n", r.why.c_str());
      return 1;
    }
    mxi::print("  searched %.1f mm either way, %zu trials, %zu peak%s "
               "refined: the score %.2f as given, %.2f at the best (of %zu "
               "at most)\n",
               options.scope_mm, r.trials, r.peaks, r.peaks == 1 ? "" : "s",
               r.score_start, r.score_best, r.candidates);
    mxi::print("Beam centre found:    %.2f, %.2f px (%.2f, %.2f mm)\n",
               r.px_after[0], r.px_after[1], r.mm_after[0], r.mm_after[1]);
    mxi::print("  shift %+.2f, %+.2f px (%+.3f, %+.3f mm)\n",
               r.px_after[0] - r.px_before[0], r.px_after[1] - r.px_before[1],
               r.mm_after[0] - r.mm_before[0], r.mm_after[1] - r.mm_before[1]);
    if (r.panel_after != r.panel_before)
      mxi::print("  on panel %zu, where it was %zu\n", r.panel_after,
                 r.panel_before);
    const rflx::Outputs out = rflx::outputs(args, "optimised");
    const json::Value document = experiments_to_json(experiments);
    rflx::write_outputs(out, &document, &strong, "mxi_search_beam");
    mxi::print("Wrote %s\n", rflx::describe(out).c_str());
    return 0;
  } catch (const std::exception &e) {
    mxi::eprint("mxi search-beam: %s\n", e.what());
    return 1;
  }
}

} // namespace
} // namespace mxi

int mxi_search_beam_main(int argc, char **argv) {
  return mxi::run_search_beam(argc, argv);
}
