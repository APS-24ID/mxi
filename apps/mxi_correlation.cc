// mxi_correlation: intensity-based clustering of isomorphous data sets, as
// dials.correlation_matrix (Thompson et al. 2025). Stage 1: the sigma-weighted
// correlation of every pair of data sets and its effective sample size, from
// a scaled file whose lattices are each a data set (docs/correlation.md).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <set>
#include <string>
#include <vector>

#include "args.hh"
#include "correlation.hh"
#include "expt.hh"
#include "json.hh"
#include "log_mirror.hh"
#include "refl.hh"
#include "rflx.hh"
#include "symmetry.hh"

namespace {

const char *kUsage =
    "usage: mxi_correlation SCALED.rflx [options]\n"
    "       mxi_correlation SCALED.expt SCALED.refl [options]\n"
    "\n"
    "The correlation of every pair of data sets -- each lattice of a scaled\n"
    "file -- weighted by the intensities' uncertainties, as\n"
    "dials.correlation_matrix: the first stage of clustering them.\n"
    "\n"
    "  --min-partiality P  reflections at least this partial (0.4)\n"
    "  --min-pairs N       a pair's correlation needs N reflections in common "
    "(3)\n"
    "  --json FILE         the matrix, its weights and the data sets, as JSON\n"
    "  -h, --help          this\n";

} // namespace

int main(int argc, char **argv) {
  using namespace mxi;
  if (!only_asks_for_help(argc, argv))
    mirror_to_log("mxi_correlation.log");
  const std::set<std::string> known = {"--min-partiality", "--min-pairs",
                                       "--json", "-h", "--help"};
  std::set<std::string> takes_value = known;
  takes_value.erase("-h");
  takes_value.erase("--help");
  Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (args.has("--help") || args.has("-h") || args.positional.empty()) {
    std::fputs(kUsage, args.has("--help") || args.has("-h") ? stdout : stderr);
    return args.positional.empty() && !args.has("--help") && !args.has("-h")
               ? 2
               : 0;
  }
  args.positional = rflx::as_pair(args.positional);
  if (args.positional.size() != 2) {
    std::fputs(kUsage, stderr);
    return 2;
  }
  try {
    const ExperimentList experiments = read_experiments(args.positional[0]);
    const Table reflections = read_reflections(args.positional[1]);
    const Experiment *first = nullptr;
    for (const Experiment &e : experiments)
      if (e.crystal) {
        first = &e;
        break;
      }
    if (first == nullptr) {
      std::fprintf(stderr, "mxi_correlation: no experiment has a crystal\n");
      return 1;
    }
    const SpaceGroup group =
        SpaceGroup::from_hall(first->crystal->space_group_hall);
    std::vector<std::size_t> which;
    const std::vector<MergedSet> sets =
        merged_by_experiment(experiments, reflections, group,
                             args.number("--min-partiality", 0.4), &which);
    const auto min_pairs =
        static_cast<std::size_t>(args.number("--min-pairs", 3.0));
    const CorrelationMatrix m = correlation_matrix(sets, min_pairs);
    const std::size_t n = m.n;
    std::printf("%zu data sets, in %s\n", n, group.name().c_str());
    std::vector<std::size_t> sizes;
    for (const MergedSet &s : sets)
      sizes.push_back(s.hkl.size());
    if (!sizes.empty()) {
      std::vector<std::size_t> sorted = sizes;
      std::sort(sorted.begin(), sorted.end());
      std::printf("merged reflections a data set: %zu to %zu, median %zu\n",
                  sorted.front(), sorted.back(), sorted[sorted.size() / 2]);
    }
    std::vector<double> rs, ns;
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = i + 1; j < n; ++j)
        if (m.w[i * n + j] > 0.0)
          rs.push_back(m.r[i * n + j]);
    const std::size_t pairs = n * (n - 1) / 2;
    std::printf("%zu of %zu pairs have a correlation\n", rs.size(), pairs);
    if (!rs.empty()) {
      std::sort(rs.begin(), rs.end());
      std::printf("correlations: %.3f to %.3f, median %.3f\n", rs.front(),
                  rs.back(), rs[rs.size() / 2]);
    }
    if (n <= 16 && n > 0) {
      std::printf("\nthe correlation matrix, by experiment:\n%6s", "");
      for (std::size_t j = 0; j < n; ++j)
        std::printf(" %6zu", which[j]);
      std::printf("\n");
      for (std::size_t i = 0; i < n; ++i) {
        std::printf("%6zu", which[i]);
        for (std::size_t j = 0; j < n; ++j) {
          if (i == j)
            std::printf(" %6s", "-");
          else if (m.w[i * n + j] > 0.0)
            std::printf(" %6.3f", m.r[i * n + j]);
          else
            std::printf(" %6s", "none");
        }
        std::printf("\n");
      }
    }
    if (args.has("--json")) {
      json::Array exps, r_rows, w_rows, merged;
      for (std::size_t i = 0; i < n; ++i) {
        exps.push_back(json::Value(static_cast<long>(which[i])));
        merged.push_back(json::Value(static_cast<long>(sizes[i])));
        json::Array rr, ww;
        for (std::size_t j = 0; j < n; ++j) {
          rr.push_back(json::Value(m.r[i * n + j]));
          ww.push_back(json::Value(m.w[i * n + j]));
        }
        r_rows.push_back(json::Value(std::move(rr)));
        w_rows.push_back(json::Value(std::move(ww)));
      }
      json::Object doc;
      doc["experiments"] = json::Value(std::move(exps));
      doc["merged_reflections"] = json::Value(std::move(merged));
      doc["r"] = json::Value(std::move(r_rows));
      doc["w"] = json::Value(std::move(w_rows));
      json::dump_file(args.value("--json", "correlation.json"),
                      json::Value(std::move(doc)));
      std::printf("\nWrote the matrix to %s\n",
                  args.value("--json", "correlation.json").c_str());
    }
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi_correlation: %s\n", e.what());
    return 1;
  }
  return 0;
}
