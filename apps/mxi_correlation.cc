// mxi_correlation: intensity-based clustering of isomorphous data sets, as
// dials.correlation_matrix (Thompson et al. 2025). Stage 1: the sigma-weighted
// correlation of every pair of data sets and its effective sample size, from
// a scaled file whose lattices are each a data set (docs/correlation.md).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
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
    "  --html FILE         the report (mxi_correlation.html): the matrices "
    "with\n"
    "                      their dendrograms, the reachability plot, the\n"
    "                      principal components, the dimension functional\n"
    "  --write-clusters    each cluster as its own .rflx, PREFIX_1.rflx and\n"
    "                      on: its data sets and their reflections, for\n"
    "                      mxi_scale\n"
    "  --prefix PREFIX     those files' names (cluster)\n"
    "  --json FILE         the matrix, its weights, the data sets, the\n"
    "                      coordinates and their principal components\n"
    "  -h, --help          this\n";

} // namespace

// A Plotly figure, as mxi_report draws them.
std::string figure(int number, const mxi::json::Array &traces,
                   const mxi::json::Object &layout, int height = 420) {
  using namespace mxi;
  json::Object l = layout;
  l["height"] = json::Value(static_cast<long>(height));
  l["margin"] = json::Value(json::Object{{"t", json::Value(40L)}});
  std::ostringstream o;
  o << "<div class=\"plot\" id=\"plot" << number
    << "\"></div><script>Plotly.newPlot('plot" << number << "', "
    << json::dump(json::Value(traces), 0) << ", "
    << json::dump(json::Value(std::move(l)), 0)
    << ", {responsive: true});</script>\n";
  return o.str();
}

mxi::json::Value text(const std::string &s) { return mxi::json::Value(s); }

// A matrix as a heat map in its Ward dendrogram's order, the dendrogram
// above it: as dials.correlation_matrix draws the correlations and the cosine
// angles. `values` n x n; the distances 1 - value.
std::string matrix_figure(int number, const std::string &title,
                          const std::vector<double> &values, std::size_t n,
                          const std::vector<std::size_t> &which) {
  using namespace mxi;
  std::vector<double> distance(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      distance[i * n + j] =
          i == j ? 0.0 : std::max(0.0, 1.0 - values[i * n + j]);
  const Dendrogram tree = dendrogram(ward_linkage(distance, n), n);
  json::Array traces, places, labels;
  for (std::size_t k = 0; k < tree.leaves.size(); ++k) {
    places.push_back(json::Value(10.0 * static_cast<double>(k) + 5.0));
    labels.push_back(text(std::to_string(which[tree.leaves[k]])));
  }
  json::Array z;
  for (std::size_t a : tree.leaves) {
    json::Array row;
    for (std::size_t b : tree.leaves)
      row.push_back(json::Value(a == b ? 1.0 : values[a * n + b]));
    z.push_back(json::Value(std::move(row)));
  }
  traces.push_back(json::Value(json::Object{
      {"type", text("heatmap")},
      {"x", json::Value(places)},
      {"y", json::Value(places)},
      {"z", json::Value(std::move(z))},
      {"colorscale", text("Viridis")},
      {"colorbar", json::Value(json::Object{{"len", json::Value(0.7)},
                                            {"y", json::Value(0.35)}})},
      {"xaxis", text("x")},
      {"yaxis", text("y")}}));
  for (std::size_t k = 0; k < tree.x.size(); ++k) {
    json::Array xs, ys;
    for (int q = 0; q < 4; ++q) {
      xs.push_back(json::Value(tree.x[k][q]));
      ys.push_back(json::Value(tree.y[k][q]));
    }
    traces.push_back(json::Value(json::Object{
        {"type", text("scatter")},
        {"mode", text("lines")},
        {"x", json::Value(std::move(xs))},
        {"y", json::Value(std::move(ys))},
        {"line", json::Value(json::Object{{"color", text("#444")},
                                          {"width", json::Value(1L)}})},
        {"hoverinfo", text("skip")},
        {"showlegend", json::Value(false)},
        {"xaxis", text("x")},
        {"yaxis", text("y2")}}));
  }
  const json::Object ticks{
      {"tickvals", json::Value(places)},
      {"ticktext", json::Value(labels)},
      {"tickfont", json::Value(json::Object{{"size", json::Value(9L)}})}};
  json::Object x = ticks, y = ticks;
  x["title"] = text("data set (experiment), in the dendrogram's order");
  y["domain"] = json::Value(json::Array{json::Value(0.0), json::Value(0.72)});
  y["autorange"] = text("reversed");
  json::Object y2{
      {"domain", json::Value(json::Array{json::Value(0.75), json::Value(1.0)})},
      {"title", text("Ward distance")},
      {"zeroline", json::Value(false)}};
  json::Object layout{{"title", text(title)},
                      {"xaxis", json::Value(std::move(x))},
                      {"yaxis", json::Value(std::move(y))},
                      {"yaxis2", json::Value(std::move(y2))},
                      {"showlegend", json::Value(false)}};
  return figure(number, traces, layout, 760);
}

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
      "--html",      "--write-clusters", "--prefix",      "--min-partiality",
      "--min-pairs", "--json",           "--dimensions",  "--by-projection",
      "--xi",        "--max-distance",   "--min-samples", "-h",
      "--help"};
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
    // Stage 4: the report.
    if (c.dim > 0) {
      const std::size_t n = c.n;
      std::ostringstream html;
      html
          << "<!DOCTYPE html><html><head><meta charset=\"utf-8\"><title>"
             "mxi_correlation: "
          << args.positional[0]
          << "</title><script src=\"https://cdn.plot.ly/plotly-2.27.0.min.js\">"
             "</script><style>body{font-family:sans-serif;margin:2em;"
             "max-width:1100px}table{border-collapse:collapse}td,th{padding:"
             "2px 10px;text-align:right}th{border-bottom:1px solid #888}"
             ".plot{margin:1em 0}</style></head><body>\n"
             "<h1>mxi_correlation</h1>\n<p>"
          << n << " data sets of " << args.positional[0] << ", in "
          << group.name()
          << ". Clustered as dials.correlation_matrix (Thompson et al. "
             "2025): the sigma-weighted correlation of every pair, cosym's "
             "coordinates in "
          << c.dim << " dimensions, OPTICS (min_samples "
          << clusters.min_samples << ").</p>\n";
      std::map<int, std::vector<std::size_t>> members;
      for (std::size_t i = 0; i < clusters.labels.size(); ++i)
        members[clusters.labels[i]].push_back(which[i]);
      html << "<h2>Clusters</h2><table><tr><th>cluster</th><th>data sets"
              "</th><th style=\"text-align:left\">experiments</th></tr>";
      for (const auto &[label, list] : members) {
        html << "<tr><td>"
             << (label == -1 ? std::string("outliers")
                             : std::to_string(label + 1))
             << "</td><td>" << list.size()
             << "</td><td style=\"text-align:left\">";
        for (std::size_t k = 0; k < list.size(); ++k)
          html << (k ? " " : "") << list[k];
        html << "</td></tr>";
      }
      html << "</table>\n";
      if (clusters.single_cluster)
        html << "<p>No min_samples gave two clusters: one cluster, as "
                "DIALS's single-cluster heuristic finds it.</p>\n";
      else if (clusters.score > -1.0)
        html << "<p>Silhouette " << clusters.silhouette << ", score "
             << clusters.score << ".</p>\n";
      int number = 0;
      std::vector<double> r = m.r;
      for (std::size_t i = 0; i < n; ++i)
        r[i * n + i] = 1.0;
      html << "<h2>The correlation matrix</h2>\n"
           << matrix_figure(++number,
                            "Correlation between data sets, sigma-weighted", r,
                            n, which);
      std::vector<double> cosine(n * n, 1.0);
      for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
          double ab = 0, aa = 0, bb = 0;
          for (std::size_t k = 0; k < c.dim; ++k) {
            ab += c.x[k * n + i] * c.x[k * n + j];
            aa += c.x[k * n + i] * c.x[k * n + i];
            bb += c.x[k * n + j] * c.x[k * n + j];
          }
          cosine[i * n + j] = aa > 0 && bb > 0 ? ab / std::sqrt(aa * bb) : 0.0;
        }
      html << "<h2>The cosine-angle matrix</h2>\n"
           << matrix_figure(++number,
                            "Cosine of the angle between data sets' cosym "
                            "coordinates",
                            cosine, n, which);
      {
        std::map<int, std::pair<json::Array, json::Array>> by;
        double top = 0.0;
        for (double v : clusters.model.reachability)
          if (std::isfinite(v))
            top = std::max(top, v);
        for (std::size_t k = 0; k < clusters.model.ordering.size(); ++k) {
          const std::size_t i = clusters.model.ordering[k];
          const double v = clusters.model.reachability[i];
          auto &[xs, ys] = by[clusters.labels[i]];
          xs.push_back(json::Value(static_cast<double>(k)));
          ys.push_back(json::Value(std::isfinite(v) ? v : 1.1 * top));
        }
        json::Array traces;
        for (auto &[label, xy] : by)
          traces.push_back(json::Value(json::Object{
              {"type", text("bar")},
              {"name",
               text(label == -1 ? std::string("outliers")
                                : "cluster " + std::to_string(label + 1))},
              {"x", json::Value(std::move(xy.first))},
              {"y", json::Value(std::move(xy.second))}}));
        html << "<h2>OPTICS</h2>\n"
             << figure(++number, traces,
                       json::Object{
                           {"title", text("Reachability, in OPTICS's order "
                                          "(infinite drawn at 1.1 x the "
                                          "largest)")},
                           {"xaxis", json::Value(json::Object{
                                         {"title", text("OPTICS's order")}})},
                           {"yaxis", json::Value(json::Object{
                                         {"title", text("reachability")}})},
                           {"bargap", json::Value(0.1)}});
      }
      if (c.dim >= 2) {
        std::map<int, std::array<json::Array, 3>> by;
        for (std::size_t i = 0; i < n; ++i) {
          auto &v = by[clusters.labels[i]];
          v[0].push_back(json::Value(c.components[i * c.dim]));
          v[1].push_back(json::Value(c.components[i * c.dim + 1]));
          v[2].push_back(text("experiment " + std::to_string(which[i])));
        }
        json::Array traces;
        for (auto &[label, v] : by)
          traces.push_back(json::Value(json::Object{
              {"type", text("scatter")},
              {"mode", text("markers")},
              {"name",
               text(label == -1 ? std::string("outliers")
                                : "cluster " + std::to_string(label + 1))},
              {"x", json::Value(std::move(v[0]))},
              {"y", json::Value(std::move(v[1]))},
              {"text", json::Value(std::move(v[2]))}}));
        std::ostringstream xt, yt;
        xt.precision(3);
        yt.precision(3);
        xt << "principal component 1 (" << 100.0 * c.explained_variance_ratio[0]
           << " per cent)";
        yt << "principal component 2 (" << 100.0 * c.explained_variance_ratio[1]
           << " per cent)";
        html << "<h2>The coordinates' principal components</h2>\n"
             << figure(++number, traces,
                       json::Object{
                           {"title", text("Each data set's cosym "
                                          "coordinates, on their first two "
                                          "principal components")},
                           {"xaxis", json::Value(json::Object{
                                         {"title", text(xt.str())}})},
                           {"yaxis", json::Value(json::Object{
                                         {"title", text(yt.str())}})}});
      }
      if (!c.functional.empty()) {
        json::Array xs, ys;
        for (std::size_t k = 0; k < c.functional.size(); ++k) {
          xs.push_back(json::Value(c.dimensions[k]));
          ys.push_back(json::Value(c.functional[k]));
        }
        json::Array traces{
            json::Value(json::Object{{"type", text("scatter")},
                                     {"mode", text("lines+markers")},
                                     {"name", text("functional")},
                                     {"x", json::Value(std::move(xs))},
                                     {"y", json::Value(std::move(ys))}})};
        const std::size_t at = std::min(c.dim, c.functional.size()) - 1;
        traces.push_back(json::Value(json::Object{
            {"type", text("scatter")},
            {"mode", text("markers")},
            {"name", text("chosen")},
            {"marker",
             json::Value(json::Object{{"size", json::Value(14L)},
                                      {"symbol", text("circle-open")}})},
            {"x",
             json::Value(json::Array{json::Value(static_cast<double>(c.dim))})},
            {"y", json::Value(json::Array{json::Value(c.functional[at])})}}));
        html << "<h2>The dimensions</h2>\n"
             << figure(++number, traces,
                       json::Object{
                           {"title", text("The cosym functional by "
                                          "dimension, with outlier "
                                          "rejection; the dimension at its "
                                          "elbow")},
                           {"xaxis", json::Value(json::Object{
                                         {"title", text("dimensions")}})},
                           {"yaxis", json::Value(json::Object{
                                         {"title", text("functional")},
                                         {"type", text("log")}})}});
      }
      html << "</body></html>\n";
      const std::string html_path =
          args.value("--html", "mxi_correlation.html");
      std::ofstream file(html_path);
      file << html.str();
      if (!file)
        throw std::runtime_error("cannot write " + html_path);
      std::printf("\nWrote the report to %s\n", html_path.c_str());
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
