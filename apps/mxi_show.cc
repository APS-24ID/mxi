// mxi_show: the state of processing a file holds (docs/output.md) -- its
// models, its reflections by what has been done to them, and the summary each
// program left in it, in the order they ran -- as dials.show shows models, and
// more.

#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "args.hh"
#include "expt.hh"
#include "processing.hh"
#include "refl.hh"
#include "rflx.hh"
#include "symmetry.hh"

namespace mxi {

namespace {

void usage(const char *program) {
  std::printf(
      "usage: %s FILE [REFL]\n\n"
      "The state of processing FILE holds: a .rflx, or an .expt with its "
      ".refl if given.\n"
      "  the models: each sweep's images and geometry, each crystal's cell "
      "and space group\n"
      "  the reflections: how many, and what has been done to them\n"
      "  the summary: each program's figures and tables, in the order they "
      "ran (a .rflx's)\n",
      program);
}

std::string number(const json::Value &v) {
  char text[64];
  if (v.is_integral())
    std::snprintf(text, sizeof text, "%lld",
                  static_cast<long long>(v.as_number()));
  else
    std::snprintf(text, sizeof text, "%.6g", v.as_number());
  return text;
}

std::string value_text(const json::Value &v) {
  if (v.is_number())
    return number(v);
  if (v.is_string())
    return v.as_string();
  if (v.is_null())
    return "";
  if (v.is_array()) {
    std::string out;
    for (const json::Value &x : v.as_array())
      out += (out.empty() ? "" : " ") + value_text(x);
    return out;
  }
  if (v.is_object()) {
    std::string out;
    for (const auto &[name, x] : v.as_object())
      out += (out.empty() ? "" : "; ") + name + " " + value_text(x);
    return out;
  }
  return "";
}

bool is_table(const json::Value &v) {
  return v.is_object() && v["columns"].is_array();
}

void show_table(const std::string &name, const json::Value &t) {
  const json::Array &columns = t["columns"].as_array();
  std::vector<std::string> heads;
  std::vector<std::vector<std::string>> cells;
  std::size_t rows = 0;
  for (const json::Value &c : columns) {
    heads.push_back(c.as_string());
    std::vector<std::string> col;
    if (t[c.as_string()].is_array())
      for (const json::Value &x : t[c.as_string()].as_array())
        col.push_back(value_text(x));
    rows = std::max(rows, col.size());
    cells.push_back(std::move(col));
  }
  std::vector<std::size_t> width(heads.size());
  for (std::size_t c = 0; c < heads.size(); ++c) {
    width[c] = heads[c].size();
    for (const std::string &s : cells[c])
      width[c] = std::max(width[c], s.size());
  }
  std::printf("    %s\n     ", name.c_str());
  for (std::size_t c = 0; c < heads.size(); ++c)
    std::printf(" %*s", static_cast<int>(width[c]), heads[c].c_str());
  std::printf("\n");
  for (std::size_t r = 0; r < rows; ++r) {
    std::printf("     ");
    for (std::size_t c = 0; c < heads.size(); ++c)
      std::printf(" %*s", static_cast<int>(width[c]),
                  r < cells[c].size() ? cells[c][r].c_str() : "");
    std::printf("\n");
  }
}

void show_models(const ExperimentList &list) {
  // A sweep at a time -- its images and geometry once -- then each experiment
  // on it: a crystal, or the sweep's own experiment, crystal-less.
  std::map<std::size_t, std::vector<std::size_t>> by_sweep;
  for (std::size_t i = 0; i < list.size(); ++i)
    by_sweep[list.sweep(i)].push_back(i);
  std::printf("Experiments: %zu, on %zu sweep%s\n", list.size(),
              by_sweep.size(), by_sweep.size() == 1 ? "" : "s");
  for (const auto &[sweep, members] : by_sweep) {
    const Experiment &e = list[members.front()];
    std::printf("  sweep %zu: images %lld to %lld, %.4f deg each from %.4f deg "
                "(%.1f deg in all); wavelength %.6f A\n",
                sweep, static_cast<long long>(e.scan.first_image),
                static_cast<long long>(e.scan.last_image), e.scan.osc_width,
                e.scan.osc_start,
                e.scan.osc_width * static_cast<double>(e.scan.num_images()),
                e.beam.wavelength);
    if (e.detector.size() > 0) {
      const Panel &p = e.detector[0];
      const Vec3 normal = p.fast.cross(p.slow);
      const double distance = std::fabs(p.origin.dot(normal / normal.norm()));
      std::printf("    detector: %zu panel%s, distance %.3f mm",
                  e.detector.size(), e.detector.size() == 1 ? "" : "s",
                  distance);
      if (const auto hit = e.detector.intersect(e.beam.s0()))
        std::printf(", the beam at pixel %.2f, %.2f", std::get<1>(*hit),
                    std::get<2>(*hit));
      std::printf("\n");
    }
    for (std::size_t i : members) {
      const Experiment &x = list[i];
      if (!x.crystal) {
        std::printf("    experiment %zu: no crystal -- the sweep's own, its "
                    "unindexed spots\n",
                    i);
        continue;
      }
      const UnitCell c = x.crystal->cell();
      std::string group = x.crystal->space_group_hall;
      try {
        group = SpaceGroup::from_hall(x.crystal->space_group_hall).name();
      } catch (const std::exception &) {
      }
      std::printf("    experiment %zu: %.3f %.3f %.3f A, %.3f %.3f %.3f deg; "
                  "%s%s\n",
                  i, c.a, c.b, c.c, c.alpha, c.beta, c.gamma, group.c_str(),
                  x.crystal->A_points.empty() ? "" : "; scan-varying");
      const ProfileBlock &profile =
          list.profiles.size() == list.size() ? list.profiles[i] : list.profile;
      if (profile.present)
        std::printf("      profile model: sigma_b %.4f deg, sigma_m %.4f deg\n",
                    profile.sigma_b, profile.sigma_m);
    }
  }
}

void show_reflections(const Table &t) {
  std::printf("Reflections: %zu\n", t.nrows);
  if (!t.has("flags") || t.nrows == 0)
    return;
  const Column &flags = t.at("flags");
  const bool has_id = t.has("id");
  struct Count {
    std::size_t rows = 0, strong = 0, indexed = 0, sum = 0, prf = 0, scaled = 0,
                outlier = 0;
  };
  std::map<long long, Count> by_id;
  for (std::size_t r = 0; r < t.nrows; ++r) {
    const std::int64_t f = flags.integer(r);
    Count &c =
        by_id[has_id ? static_cast<long long>(t.at("id").integer(r)) : 0];
    ++c.rows;
    c.strong += (f >> 5) & 1;
    c.indexed += (f >> 2) & 1;
    c.sum += (f >> 8) & 1;
    c.prf += (f >> 9) & 1;
    c.scaled += (f >> 26) & 1;
    c.outlier += (f >> 23) & 1;
  }
  std::printf("  %4s %9s %9s %9s %9s %9s %9s %9s\n", "id", "rows", "strong",
              "indexed", "summed", "fitted", "scaled", "outliers");
  for (const auto &[id, c] : by_id)
    std::printf("  %4lld %9zu %9zu %9zu %9zu %9zu %9zu %9zu\n", id, c.rows,
                c.strong, c.indexed, c.sum, c.prf, c.scaled, c.outlier);
}

void show_summary(const json::Value &summary) {
  const std::vector<json::Value> steps = Summary::steps(summary);
  std::printf("Processing: %zu step%s\n", steps.size(),
              steps.size() == 1 ? "" : "s");
  for (const json::Value &step : steps) {
    std::printf("  %s, %s\n", step["program"].as_string().c_str(),
                step["time"].is_string() ? step["time"].as_string().c_str()
                                         : "");
    for (const auto &[name, value] : step.as_object()) {
      if (name == "program" || name == "time" || is_table(value))
        continue;
      std::printf("    %s: %s\n", name.c_str(), value_text(value).c_str());
    }
    for (const auto &[name, value] : step.as_object())
      if (is_table(value))
        show_table(name, value);
  }
}

} // namespace

int run_show(int argc, char **argv) {
  const std::set<std::string> known = {"--help"};
  const std::set<std::string> takes_value = {};
  Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (args.has("--help") || args.positional.empty() ||
      args.positional.size() > 2) {
    usage(argv[0]);
    return args.has("--help") ? 0 : 2;
  }
  try {
    const std::string &first = args.positional[0];
    std::printf("%s\n", first.c_str());
    const ExperimentList list = read_experiments(first);
    show_models(list);
    const std::string refl =
        args.positional.size() == 2 ? args.positional[1] : first;
    if (refl != first || rflx::is_hdf5(first)) {
      try {
        show_reflections(read_reflections(refl));
      } catch (const std::exception &) {
        // An experiment list alone: nothing more to show.
      }
    }
    if (rflx::is_hdf5(first)) {
      if (const auto summary = rflx::read_summary(first))
        show_summary(*summary);
      else
        std::printf("Processing: no summary in the file\n");
    }
    return 0;
  } catch (const std::exception &error) {
    std::fprintf(stderr, "mxi_show: %s\n", error.what());
    return 1;
  }
}

} // namespace mxi

// The program's entry, run by mxi (show) -- apps/mxi.cc.
int mxi_show_main(int argc, char **argv) { return mxi::run_show(argc, argv); }
