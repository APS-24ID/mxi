// mxi_correlation: intensity-based clustering of isomorphous data sets, as
// dials.correlation_matrix (Thompson et al. 2025). Stage 1: the sigma-weighted
// correlation of every pair of data sets and its effective sample size, from
// a scaled file whose lattices are each a data set (docs/correlation.md).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <optional>
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
#include "views.hh"

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
    "  --dimensions N      the cosym coordinates in N dimensions (by default\n"
    "                      the elbow of the functional over 1 to\n"
    "                      min(n, 50), as dials.correlation_matrix)\n"
    "  --by-projection     the functional for each dimension from one\n"
    "                      minimisation, projected (DIALS #3248): faster\n"
    "  --xi X              OPTICS's steepness for a cluster's edge (0.05)\n"
    "  --max-distance D    OPTICS's max_eps (0.5, as DIALS's)\n"
    "  --min-samples N     OPTICS's min_samples, not searched for\n"
    "  --write-clusters    each cluster as its own .rflx, PREFIX_1.rflx and\n"
    "                      on: its data sets and their reflections, for\n"
    "                      mxi_scale\n"
    "  --prefix PREFIX     those files' names (cluster)\n"
    "  --json FILE         the matrix, its weights, the data sets, the\n"
    "                      coordinates and their principal components\n"
    "  -h, --help          this\n";

} // namespace

std::size_t members_count(const std::vector<int> &labels, int label) {
  std::size_t n = 0;
  for (int l : labels)
    n += l == label;
  return n;
}

int main(int argc, char **argv) {
  using namespace mxi;
  if (!only_asks_for_help(argc, argv))
    mirror_to_log("mxi_correlation.log");
  const std::set<std::string> known = {
      "--write-clusters", "--prefix", "--min-partiality",
      "--min-pairs",      "--json",   "--dimensions",
      "--by-projection",  "--xi",     "--max-distance",
      "--min-samples",    "-h",       "--help"};
  std::set<std::string> takes_value = known;
  takes_value.erase("--by-projection");
  takes_value.erase("--write-clusters");
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
    // Stage 2: the coordinates, their dimension, their principal components.
    CosymOptions options;
    options.dimensions =
        static_cast<std::size_t>(args.number("--dimensions", 0.0));
    options.by_projection = args.has("--by-projection");
    const CosymCoordinates c = cosym_coordinates(m, options);
    if (c.dim > 0) {
      if (!c.functional.empty()) {
        std::printf("\nthe cosym functional by dimension%s:\n",
                    options.by_projection ? " (by projection)" : "");
        for (std::size_t k = 0; k < c.functional.size(); ++k)
          std::printf("  %3.0f  %.6g\n", c.dimensions[k], c.functional[k]);
      }
      std::printf("\n%zu dimensions%s; the functional there %.6g\n", c.dim,
                  options.dimensions > 0 ? " (as asked)" : ", at the elbow",
                  c.final_functional);
      std::printf("the variance each principal component explains:");
      for (std::size_t k = 0; k < c.explained_variance_ratio.size(); ++k)
        std::printf(" %.3f", c.explained_variance_ratio[k]);
      std::printf("\n");
    }
    // Stage 3: OPTICS on the coordinates, the clusters and the outliers.
    SignificantClusters clusters;
    if (c.dim > 0) {
      ClusterOptions co;
      co.xi = args.number("--xi", 0.05);
      co.max_eps = args.number("--max-distance", 0.5);
      if (args.has("--min-samples")) {
        co.optimise = false;
        co.min_points_buffer = 0.0; // the initial guess replaced below
      }
      CosymCoordinates for_clusters = c;
      clusters = significant_clusters(for_clusters, co);
      if (args.has("--min-samples")) {
        std::vector<double> points(c.n * c.dim);
        for (std::size_t i = 0; i < c.n; ++i)
          for (std::size_t k = 0; k < c.dim; ++k)
            points[i * c.dim + k] = c.x[k * c.n + i];
        const auto ms =
            static_cast<std::size_t>(args.number("--min-samples", 5));
        clusters.model = optics(points, c.n, c.dim, ms, co.max_eps, co.xi);
        clusters.labels = clusters.model.labels;
        clusters.min_samples = ms;
      }
      std::printf(
          "\nOPTICS (Ankerst et al. 1999), xi %.3g, max distance %.3g: ", co.xi,
          co.max_eps);
      if (args.has("--min-samples"))
        std::printf("min_samples %zu, as asked\n", clusters.min_samples);
      else if (c.n < 5)
        std::printf("fewer than five data sets, min_samples %zu, as DIALS's\n",
                    clusters.min_samples);
      else if (clusters.single_cluster)
        std::printf("no min_samples gave two clusters; one cluster, as DIALS's "
                    "single-cluster heuristic finds it, min_samples %zu\n",
                    clusters.min_samples);
      else
        std::printf("min_samples %zu (the initial guess %zu), silhouette %.3f, "
                    "score %.3f\n",
                    clusters.min_samples, clusters.initial_min_samples,
                    clusters.silhouette, clusters.score);
      std::map<int, std::vector<std::size_t>> members;
      for (std::size_t i = 0; i < clusters.labels.size(); ++i)
        members[clusters.labels[i]].push_back(which[i]);
      for (const auto &[label, list] : members) {
        if (label == -1)
          continue;
        std::printf("  cluster %d, %zu data sets:", label + 1, list.size());
        for (std::size_t e : list)
          std::printf(" %zu", e);
        std::printf("\n");
      }
      if (members.count(-1)) {
        std::printf("  outliers, %zu:", members[-1].size());
        for (std::size_t e : members[-1])
          std::printf(" %zu", e);
        std::printf("\n");
      }
    }
    // Stage 4: each cluster written as its own file -- its data sets, any
    // crystal-less sweep experiment of theirs (the same image set and scan),
    // and their reflections, the ids renumbered (keep_experiments).
    if (args.has("--write-clusters") && !clusters.labels.empty()) {
      const json::Value document = read_experiment_document(args.positional[0]);
      const json::Array &all = document["experiment"].as_array();
      std::optional<json::Value> summary;
      if (rflx::has_experiments(args.positional[0]))
        summary = rflx::read_summary(args.positional[0]);
      const std::string prefix = args.value("--prefix", "cluster");
      std::set<int> labels(clusters.labels.begin(), clusters.labels.end());
      labels.erase(-1);
      std::printf("\n");
      for (int label : labels) {
        std::vector<bool> keep(all.size(), false);
        std::set<std::pair<long, long>> sweeps; // image set, scan
        for (std::size_t i = 0; i < clusters.labels.size(); ++i)
          if (clusters.labels[i] == label) {
            keep[which[i]] = true;
            const json::Value &x = all[which[i]];
            if (x["imageset"].is_number() && x["scan"].is_number())
              sweeps.insert({static_cast<long>(x["imageset"].as_number()),
                             static_cast<long>(x["scan"].as_number())});
          }
        for (std::size_t e = 0; e < all.size(); ++e)
          if (!all[e]["crystal"].is_number() &&
              all[e]["imageset"].is_number() && all[e]["scan"].is_number() &&
              sweeps.count({static_cast<long>(all[e]["imageset"].as_number()),
                            static_cast<long>(all[e]["scan"].as_number())}))
            keep[e] = true;
        std::vector<std::size_t> rows;
        if (reflections.has("id")) {
          const Column &id = reflections.at("id");
          for (std::size_t r = 0; r < reflections.nrows; ++r) {
            const std::int64_t e = id.integer(r);
            if (e >= 0 && static_cast<std::size_t>(e) < keep.size() &&
                keep[static_cast<std::size_t>(e)])
              rows.push_back(r);
          }
        }
        Table table = select_rows(reflections, rows);
        const json::Value sub = keep_experiments(document, keep, &table);
        const std::string path =
            prefix + "_" + std::to_string(label + 1) + ".rflx";
        rflx::write(path, &sub, &table, "mxi_correlation", false,
                    summary ? &*summary : nullptr);
        std::printf("Wrote %s: cluster %d, %zu data sets, %zu reflections\n",
                    path.c_str(), label + 1,
                    members_count(clusters.labels, label), table.nrows);
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
      json::Array dims, func, comps, ratio;
      for (std::size_t k = 0; k < c.dimensions.size(); ++k) {
        dims.push_back(json::Value(c.dimensions[k]));
        func.push_back(json::Value(c.functional[k]));
      }
      for (std::size_t i = 0; i < c.n && c.dim > 0; ++i) {
        json::Array row;
        for (std::size_t k = 0; k < c.dim; ++k)
          row.push_back(json::Value(c.components[i * c.dim + k]));
        comps.push_back(json::Value(std::move(row)));
      }
      for (double v : c.explained_variance_ratio)
        ratio.push_back(json::Value(v));
      json::Array coords;
      for (double v : c.x)
        coords.push_back(json::Value(v));
      doc["dimension"] = json::Value(static_cast<long>(c.dim));
      doc["scan_dimensions"] = json::Value(std::move(dims));
      doc["scan_functional"] = json::Value(std::move(func));
      doc["coordinates"] = json::Value(std::move(coords));
      doc["components"] = json::Value(std::move(comps));
      doc["explained_variance_ratio"] = json::Value(std::move(ratio));
      doc["functional"] = json::Value(c.final_functional);
      json::Array labels, order, reach, search;
      for (int l : clusters.labels)
        labels.push_back(json::Value(static_cast<long>(l)));
      for (std::size_t k : clusters.model.ordering) {
        order.push_back(json::Value(static_cast<long>(k)));
        reach.push_back(json::Value(clusters.model.reachability[k]));
      }
      for (const auto &[ms, score] : clusters.scores) {
        json::Array pair;
        pair.push_back(json::Value(static_cast<long>(ms)));
        pair.push_back(json::Value(score));
        search.push_back(json::Value(std::move(pair)));
      }
      doc["labels"] = json::Value(std::move(labels));
      doc["ordering"] = json::Value(std::move(order));
      doc["reachability_in_order"] = json::Value(std::move(reach));
      doc["min_samples"] = json::Value(static_cast<long>(clusters.min_samples));
      doc["min_samples_scores"] = json::Value(std::move(search));
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
