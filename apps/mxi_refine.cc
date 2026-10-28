// mxi_refine: refine the experimental model against indexed spot centroids.
//
//   mxi_refine indexed.expt indexed.refl

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <string>

#include "../src/expt.hh"
#include "../src/refl.hh"
#include "../src/shoebox.hh"
#include "args.hh"
#include "combine.hh"
#include "log_mirror.hh"
#include "postrefine.hh"
#include "processing.hh"
#include "refine.hh"
#include "rflx.hh"
#include "timing.hh"
#include "views.hh"
#include <algorithm>
#include <vector>

namespace mxi {

namespace {
void usage() {
  std::printf(
      "usage: mxi_refine INDEXED.expt INDEXED.refl [options]\n"
      "  --no-crystal      hold the crystal fixed\n"
      "  --jacobian-threads N  threads for the Jacobian; 0 is one per core "
      "(0)\n"
      "  --no-detector     hold the detector fixed\n"
      "  --no-beam         hold the beam direction fixed (refined by default;\n"
      "                    on a single sweep it is correlated with the\n"
      "                    detector)\n"
      "  --joint           several sweeps refined in one fit, sharing what\n"
      "                    models they share, where by default each sweep is\n"
      "                    refined on its own\n"
      "  --shared-crystal  several sweeps: one crystal for them all, where by\n"
      "                    default each refines its own (--separate, as was)\n"
      "  --conditional-depth  mean depth given absorption, not eqn (6)\n"
      "  --scan-varying [N]  control points in A across the scan. "
      "Scan-varying\n"
      "                    is the default, one per 10 degrees and at least "
      "five,\n"
      "                    on a scan of 10 degrees or more; N sets the number\n"
      "  --interval DEG    scan-varying control points one per DEG degrees of\n"
      "                    each sweep (10, at least five by default; at least\n"
      "                    four here) -- 20 or 30 for sparse data\n"
      "  --sample-offset   the centre of diffraction off the rotation axis,\n"
      "                    turning with the crystal: two parameters a lattice\n"
      "                    (amplitude and phase), where the crystal is not\n"
      "                    centred or the beam not on it\n"
      "  --crystal-restraint W  each control point's A tied to the static A, "
      "a\n"
      "                    prior of about 1/sqrt(W) of A's size: where the\n"
      "                    reflections hold a point it is theirs, where they "
      "do\n"
      "                    not the crystal stays static (1e6; 0 for none)\n"
      "  --static          a crystal that does not move: one A for the scan\n"
      "  --unit-weights    ignore the centroid variances\n"
      "  --strong-only     build the model from the stronger half only\n"
      "  --numerical       finite differences, not the analytical derivatives\n"
      "                    (the default), for testing\n"
      "  --timing          where the time went, by phase\n"
      "  --normal-threads N  threads for the normal equations; 0 is one per\n"
      "                   core (0). A reduction, so a threaded run differs "
      "from\n"
      "                   a serial one in the last bits; 1 to avoid that\n"
      "  --detector-in-scan-varying  keep refining the detector during the\n"
      "                    scan-varying pass; it is degenerate with the cell\n"
      "  --min-volume V    drop reflections whose rotation angle is not\n"
      "                    determined by the data; 0 keeps them all (0.05)\n"
      "  --z-weight W      scale the weight on the rotation-angle residual\n"
      "  --macrocycles N   (3)\n"
      "  --outlier-sigma S (4; 0 disables rejection)\n"
      "  -o PATH           where to write: one .rflx (refined.rflx), or with\n"
      "                    a .refl the DIALS pair, as --output-expt, -refl\n"
      "  --output-expt P   the DIALS pair's experiment list (refined.expt)\n"
      "  --output-refl P   the DIALS pair's reflections (refined.refl)\n");
}
} // namespace

namespace {
double now_wall() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}
} // namespace

namespace {

// A refinement pass as a table: each macrocycle's kept reflections, the
// outliers rejected after its fit, and its RMSDs. Printed by the program from
// what refine() records; the library used to print these lines itself.
void print_cycles(const mxi::RefineResult &r) {
  std::printf("  %5s %11s %8s %8s %8s %8s\n", "cycle", "reflections",
              "rejected", "RMSD x", "RMSD y", "RMSD z");
  std::printf("  %5s %11s %8s %8s %8s %8s\n", "", "", "", "(px)", "(px)",
              "(images)");
  for (std::size_t c = 0; c < r.cycles.size(); ++c) {
    const mxi::RefineCycle &k = r.cycles[c];
    std::printf("  %5zu %11zu %8zu %8.3f %8.3f %8.3f\n", c + 1, k.reflections,
                k.rejected, k.rmsd_x, k.rmsd_y, k.rmsd_z);
  }
}

double quantile(std::vector<double> v, double q) {
  if (v.empty())
    return 0.0;
  std::sort(v.begin(), v.end());
  const double at = q * static_cast<double>(v.size() - 1);
  const std::size_t lo = static_cast<std::size_t>(at);
  const std::size_t hi = std::min(lo + 1, v.size() - 1);
  return v[lo] + (at - static_cast<double>(lo)) * (v[hi] - v[lo]);
}

// Calculated minus observed over the reflections refined on, as
// dials.refine's summary statistics give them -- in pixels and images, the
// units of this program's RMSDs, where DIALS uses millimetres and degrees. A
// median well away from zero is a systematic the RMSD cannot show.
void print_residuals(const mxi::Table &t,
                     const std::vector<std::size_t> &unpredicted) {
  if (!t.has("xyzcal.px") || !t.has("xyzobs.px.value") || !t.has("flags"))
    return;
  std::vector<bool> skip(t.nrows, false);
  for (std::size_t i : unpredicted)
    if (i < t.nrows)
      skip[i] = true;
  std::size_t left_out = 0;
  const mxi::Column &cal = t.at("xyzcal.px");
  const mxi::Column &obs = t.at("xyzobs.px.value");
  const mxi::Column &flags = t.at("flags");
  std::vector<double> d[3];
  for (std::size_t i = 0; i < t.nrows; ++i) {
    if ((flags.ints[i] & mxi::flag::kUsedInRefinement) == 0)
      continue;
    // No prediction was written for it, so its calculated position is the
    // column's default and not a prediction: a residual from it is fiction.
    if (skip[i]) {
      ++left_out;
      continue;
    }
    for (int k = 0; k < 3; ++k)
      d[k].push_back(cal.reals[i * 3 + k] - obs.reals[i * 3 + k]);
  }
  if (d[0].empty())
    return;
  std::printf(
      "\nCalculated minus observed, over the %zu reflections refined on",
      d[0].size());
  if (left_out > 0)
    std::printf(
        " (%zu more could not be predicted with the final model, and are "
        "left out)",
        left_out);
  std::printf(":\n");
  std::printf("  %-12s %9s %9s %9s %9s %9s\n", "", "min", "Q1", "median", "Q3",
              "max");
  const char *names[3] = {"x (px)", "y (px)", "z (images)"};
  for (int k = 0; k < 3; ++k) {
    std::printf("  %-12s %9.4f %9.4f %9.4f %9.4f %9.4f\n", names[k],
                quantile(d[k], 0.0), quantile(d[k], 0.25), quantile(d[k], 0.5),
                quantile(d[k], 0.75), quantile(d[k], 1.0));
  }
}

} // namespace

// The default weight of --crystal-restraint.
constexpr double kCrystalRestraint = 1e6;

int run_refine(int argc, char **argv);

// Each independent group refined in a run of its own, on files of its own, as
// mxi_index --no-joint and mxi_integrate's sweeps are; then the lists joined
// in the original order, a slice an experiment, and the reflections' ids
// mapped back. Each group's refinement figures kept in the summary.
int refine_apart(const Arguments &args,
                 const std::set<std::string> &takes_value,
                 const std::vector<std::vector<std::size_t>> &groups) {
  const rflx::Outputs outputs = rflx::outputs(args, "refined");
  const std::string base = outputs.rflx.empty() ? outputs.refl : outputs.rflx;
  const json::Value document = read_experiment_document(args.positional[0]);
  const std::size_t n = document["experiment"].as_array().size();
  const Table indexed = read_reflections(args.positional[1]);
  if (!indexed.has("id")) {
    std::fprintf(stderr, "mxi_refine: %s has no id column\n",
                 args.positional[1].c_str());
    return 1;
  }
  std::vector<std::string> made;
  const auto tidy = [&] {
    for (const std::string &f : made)
      std::remove(f.c_str());
  };
  const auto run = [](std::vector<std::string> words) {
    std::vector<char *> argv;
    for (std::string &w : words)
      argv.push_back(w.data());
    argv.push_back(nullptr);
    return run_refine(static_cast<int>(words.size()), argv.data());
  };
  std::printf("Refining %zu experiments sweep by sweep, %zu runs, each "
              "sweep's beam and detector its own (--joint for one refinement "
              "of them all)\n",
              n, groups.size());
  std::vector<json::Value> lists(groups.size());
  std::vector<Table> tables(groups.size());
  json::Array steps;
  const Column &id = indexed.at("id");
  std::size_t unassigned = 0;
  for (std::size_t r = 0; r < indexed.nrows; ++r)
    if (id.integer(r) < 0 || static_cast<std::size_t>(id.integer(r)) >= n)
      ++unassigned;
  for (std::size_t g = 0; g < groups.size(); ++g) {
    const std::vector<std::size_t> &members = groups[g];
    const std::string stem = base + ".group" + std::to_string(g);
    const std::string in_expt = stem + ".in.expt", in_refl = stem + ".in.refl",
                      got = stem + ".rflx";
    made.insert(made.end(), {in_expt, in_refl, got});
    std::vector<json::Value> slices;
    for (std::size_t i : members)
      slices.push_back(slice_experiment_list(document, i));
    json::Value own = join_experiment_lists(slices);
    for (const char *kind :
         {"beam", "detector", "goniometer", "scan", "imageset"})
      share_identical(&own, kind);
    json::dump_file(in_expt, own);
    std::vector<long> local(n, -1);
    for (std::size_t k = 0; k < members.size(); ++k)
      local[members[k]] = static_cast<long>(k);
    std::vector<std::size_t> rows;
    for (std::size_t r = 0; r < indexed.nrows; ++r) {
      const std::int64_t i = id.integer(r);
      if (i >= 0 && static_cast<std::size_t>(i) < n &&
          local[static_cast<std::size_t>(i)] >= 0)
        rows.push_back(r);
    }
    Table part = select_rows_with_shoeboxes(indexed, rows);
    for (std::int64_t &x : part.modify_int_column("id", "int", 1).ints)
      x = local[static_cast<std::size_t>(x)];
    part.identifiers.clear();
    write_reflections(in_refl, part);
    std::vector<std::string> words = {"mxi_refine", in_expt, in_refl};
    for (const auto &[flag, value] : args.options) {
      if (flag == "-o" || flag == "--output-expt" || flag == "--output-refl")
        continue;
      words.push_back(flag);
      if (takes_value.count(flag))
        words.push_back(value);
    }
    words.insert(words.end(), {"-o", got});
    std::printf("\n=== Group %zu of %zu: experiment%s", g + 1, groups.size(),
                members.size() == 1 ? "" : "s");
    for (std::size_t i : members)
      std::printf(" %zu", i);
    std::printf(" ===\n");
    std::fflush(stdout);
    if (run(words) != 0) {
      std::fprintf(stderr, "mxi_refine: group %zu did not refine\n", g);
      tidy();
      return 1;
    }
    lists[g] = *rflx::read_experiments(got);
    tables[g] = read_reflections(got);
    if (const auto s = rflx::read_summary(got)) {
      const std::vector<json::Value> theirs = Summary::steps(*s);
      if (!theirs.empty()) {
        json::Object entry = theirs.back().as_object();
        json::Array which;
        for (std::size_t i : members)
          which.push_back(json::Value(static_cast<long>(i)));
        entry["experiments"] = json::Value(std::move(which));
        steps.push_back(json::Value(std::move(entry)));
      }
    }
  }
  // The lists joined in the original order, a slice an experiment.
  std::vector<std::pair<std::size_t, std::size_t>> where(n); // group, local
  for (std::size_t g = 0; g < groups.size(); ++g)
    for (std::size_t k = 0; k < groups[g].size(); ++k)
      where[groups[g][k]] = {g, k};
  std::vector<json::Value> slices;
  for (std::size_t i = 0; i < n; ++i)
    slices.push_back(
        slice_experiment_list(lists[where[i].first], where[i].second));
  json::Value joined_expt = join_experiment_lists(slices);
  for (const char *kind :
       {"beam", "detector", "goniometer", "scan", "imageset"})
    share_identical(&joined_expt, kind);
  // The tables, each group's ids mapped back to the whole.
  for (std::size_t g = 0; g < groups.size(); ++g)
    for (std::int64_t &x : tables[g].modify_int_column("id", "int", 1).ints)
      if (x >= 0 && static_cast<std::size_t>(x) < groups[g].size())
        x = static_cast<std::int64_t>(groups[g][static_cast<std::size_t>(x)]);
  Table joined = concat_rows_with_shoeboxes(tables);
  joined.identifiers.clear();
  const json::Array &all = joined_expt["experiment"].as_array();
  for (std::size_t k = 0; k < all.size(); ++k)
    if (all[k]["identifier"].is_string())
      joined.identifiers[k] = all[k]["identifier"].as_string();
  Summary step("mxi_refine", args.positional[0]);
  step.set("apart", json::Value(static_cast<long>(groups.size())));
  step.set("groups", json::Value(std::move(steps)));
  const json::Value summary_document = step.document();
  rflx::write_outputs(outputs, &joined_expt, &joined, "mxi_refine",
                      &summary_document);
  tidy();
  if (unassigned > 0)
    std::printf("\n%zu reflections with no experiment left out\n", unassigned);
  std::printf("\n%zu experiments refined in %zu groups, each on its own: %s\n",
              n, groups.size(), rflx::describe(outputs).c_str());
  return 0;
}

int run_refine(int argc, char **argv) {
  const double t_start = now_wall();
  double t_read = 0.0;
  double t_write = 0.0;
  const std::set<std::string> known = {"--static",
                                       "--verbose",
                                       "--no-crystal",
                                       "--no-detector",
                                       "--beam",
                                       "--no-beam",
                                       "--sample-offset",
                                       "--separate",
                                       "--shared-crystal",
                                       "--macrocycles",
                                       "--outlier-sigma",
                                       "-o",
                                       "--output-expt",
                                       "--output-refl",
                                       "--conditional-depth",
                                       "--scan-varying",
                                       "--unit-weights",
                                       "--strong-only",
                                       "--z-weight",
                                       "--analytic",
                                       "--numerical",
                                       "--min-volume",
                                       "--detector-in-scan-varying",
                                       "--jacobian-threads",
                                       "--timing",
                                       "--normal-threads",
                                       "--joint",
                                       "--interval",
                                       "--crystal-restraint"};
  const std::set<std::string> takes_value = {
      "--macrocycles", "--outlier-sigma",    "-o",
      "--output-expt", "--output-refl",      "--z-weight",
      "--min-volume",  "--jacobian-threads", "--normal-threads",
      "--interval",    "--crystal-restraint"};
  Arguments args =
      parse_arguments(argc, argv, known, takes_value, {"--scan-varying"});
  // One .rflx stands for the experiment list and the reflections
  // (docs/rflx.md).
  args.positional = rflx::as_pair(args.positional);
  // 0 means one per core, 1 means none. Exposed because a threading change
  // that cannot be switched off cannot be measured against its absence.
  g_normal_threads =
      static_cast<std::size_t>(args.number("--normal-threads", 0.0));
  g_jacobian_threads =
      static_cast<std::size_t>(args.number("--jacobian-threads", 0.0));
  if (args.help) {
    usage();
    return 0;
  }
  if (!args.ok) {
    std::fprintf(stderr, "mxi_refine: %s\n", args.error.c_str());
    return 2;
  }
  if (args.has("--static") && args.has("--scan-varying")) {
    std::fprintf(
        stderr,
        "mxi_refine: --static and --scan-varying contradict each other\n");
    return 2;
  }
  if (args.positional.size() != 2) {
    std::fprintf(stderr,
                 "mxi_refine: expected an .expt and a .refl, got %zu file "
                 "arguments\n",
                 args.positional.size());
    usage();
    return 2;
  }
  // Each sweep refined on its own (refine_apart) -- even one crystal's
  // sweeps, a shared beam or detector a copy each: better results throughout
  // (Graeme, 18 October) -- unless asked for one refinement of all (--joint)
  // or a crystal shared across sweeps (--shared-crystal), which joins them.
  if (!args.has("--joint") && !args.has("--shared-crystal")) {
    const json::Value document = read_experiment_document(args.positional[0]);
    if (document.is_object() && document["experiment"].is_array() &&
        document["experiment"].as_array().size() > 1) {
      const std::vector<std::vector<std::size_t>> groups =
          independent_groups(document);
      if (groups.size() > 1)
        return refine_apart(args, takes_value, groups);
    }
  }

  RefineOptions options;
  options.crystal = !args.has("--no-crystal");
  options.detector = !args.has("--no-detector");
  // The beam and analytical derivatives by default (Graeme); --beam and
  // --analytic still taken, from older scripts.
  options.beam = !args.has("--no-beam");
  options.sample_offset = args.has("--sample-offset");
  // Several sweeps refine a crystal each by default, as DIALS refines them:
  // indexing put them in one basis, and one matrix cannot fit sweeps no
  // goniometer quite returns between. --separate was the way to ask for this
  // and is still taken.
  options.shared_crystal = args.has("--shared-crystal");
  // The tie of the scan-varying crystal to the static one (src/refine.hh);
  // the default set after measuring it (docs/refinement.md).
  options.crystal_restraint =
      args.number("--crystal-restraint", kCrystalRestraint);
  options.unit_weights = args.has("--unit-weights");
  options.strong_only = args.has("--strong-only");
  options.analytic = !args.has("--numerical");
  options.macrocycles = static_cast<int>(args.number("--macrocycles", 3));
  options.outlier_sigma = args.number("--outlier-sigma", 4.0);
  options.z_weight = args.number("--z-weight", 1.0);
  options.min_volume = args.number("--min-volume", 0.05);
  const bool conditional_depth = args.has("--conditional-depth");
  const rflx::Outputs out = rflx::outputs(args, "refined");

  try {
    const double t_read_start = now_wall();
    // The lattices refined; each sweep's own crystal-less experiment and its
    // unindexed spots carried through, the sweep's shared models refined
    // into them (docs/multi-lattice.md, stage 2).
    const ExperimentList whole = read_experiments(args.positional[0]);
    Table whole_table = read_reflections(args.positional[1]);
    const Lattices views = lattices_of(whole);
    ExperimentList experiments = views.list;
    // The crystal-less rows apart, the rest moved, not copied: one table
    // held, not three (Graeme: symmetry 22 GB, scaling 26, on 36000
    // images).
    const Table crystal_less = crystal_less_rows(whole_table, views);
    Table reflections = lattice_rows(std::move(whole_table), views);
    t_read = now_wall() - t_read_start;
    // Scan-varying by default -- nearly every real crystal moves -- one control
    // point per 10 degrees of each sweep's own scan; static if asked, or on a
    // scan under 10 degrees, too little rotation to tell a moving crystal from
    // noise. Each sweep its own: the first sweep's count for all gave a 120
    // degree sweep the 35 of the 350 degree ones beside it.
    std::vector<std::size_t> points_of;
    for (const Experiment &e : experiments) {
      const double degrees =
          std::abs(e.scan.osc_width) * static_cast<double>(e.scan.num_images());
      std::size_t points = 1;
      if (args.has("--static")) {
        points = 1;
      } else if (args.has("--scan-varying") &&
                 !args.value("--scan-varying", "").empty()) {
        points = static_cast<std::size_t>(args.number("--scan-varying", 1));
      } else if (degrees >= 10.0 && args.has("--interval")) {
        points =
            scan_varying_points(e.scan, args.number("--interval", 10.0), 4);
      } else if (degrees >= 10.0) {
        points = scan_varying_points(e.scan);
      } else {
        std::printf("Static: the scan is %.1f degrees, under the 10 a "
                    "scan-varying crystal "
                    "needs\n",
                    degrees);
      }
      points_of.push_back(points);
    }
    // One count, as one sweep has: printed as before; several, each.
    std::string points_text = std::to_string(points_of.front());
    if (std::any_of(points_of.begin(), points_of.end(),
                    [&](std::size_t p) { return p != points_of.front(); })) {
      points_text.clear();
      for (std::size_t p : points_of)
        points_text += (points_text.empty() ? "" : ", ") + std::to_string(p);
    }
    if (conditional_depth) {
      for (Experiment &e : experiments) {
        for (Panel &p : e.detector.panels)
          p.parallax_conditional = true;
      }
      detail("using the conditional absorption depth\n");
    }
    set_verbose(args.has("--verbose"));
    Summary summary("mxi_refine", args.positional[0]);
    const ExperimentList given_models = experiments;
    std::printf("Refining against %zu reflections from %zu experiment%s\n",
                reflections.nrows, experiments.size(),
                experiments.size() == 1 ? "" : "s");

    const TwoPassRefinement two =
        refine_in_two_passes(experiments, reflections, options, points_of,
                             !args.has("--detector-in-scan-varying"));
    if (two.static_pass.n_used > 0) {
      std::printf("\nScan-static: %zu parameters\n",
                  two.static_pass.n_parameters);
      print_cycles(two.static_pass);
    }
    if (two.varied && two.varying_pass.n_used > 0) {
      std::printf("\nScan-varying, %s control points: %zu parameters%s\n",
                  points_text.c_str(), two.varying_pass.n_parameters,
                  two.detector_held
                      ? "; the detector held where the scan-static pass put "
                        "it (--detector-in-scan-varying to refine it too)"
                      : "");
      print_cycles(two.varying_pass);
    }
    RefineResult result = two.final();
    if (result.n_used == 0) {
      std::fprintf(stderr, "mxi_refine: nothing to refine against\n");
      return 1;
    }
    std::printf(
        "\nRefined on %zu of %zu reflections: %zu rejected as outliers, "
        "%zu with an undetermined rotation angle; %d steps\n",
        result.n_used, reflections.nrows, result.n_rejected,
        result.n_ill_conditioned, result.iterations);
    // Said plainly, because a residual averaged over a different set of
    // reflections is not comparable with anything -- including dials.refine,
    // which applies the same cutoff and reports over what is left.
    std::printf("RMSD over those %zu: %.4f px, %.4f px, %.4f images\n",
                result.n_used, result.rmsd_x, result.rmsd_y, result.rmsd_z);

    // s1, rlp and entering follow the refined model, as they do in DIALS.
    // xyzobs.mm deliberately does not: it is what the spot finder measured
    // through the model as imported, and recomputing it here would silently
    // change the observations refinement was just fitted to.
    set_refinement_flags(result, reflections);
    add_reciprocal_columns(experiments, reflections);
    const std::vector<std::size_t> unpredicted =
        update_predictions(experiments, reflections);
    print_residuals(reflections, unpredicted);
    std::printf("\n");
    // Each crystal: the cell with the scan-static pass's e.s.d.s, as DIALS
    // prints it; scan-varying, how far the cell and the orientation move
    // over the scan -- the check the tutorial makes in dials.report.
    std::vector<json::Value> t_exp, t_cell, t_lo, t_hi, t_turn;
    // Each crystal's cell at every scan point, as refinement left it: the
    // report draws these, not the file's later crystals, which symmetry
    // constrains (Graeme: a cubic cell's angles all 90.000 by then).
    json::Array over_scan;
    for (std::size_t i = 0; i < experiments.size(); ++i) {
      if (!experiments[i].crystal)
        continue;
      const Crystal &x = *experiments[i].crystal;
      const UnitCell c = x.cell();
      const std::array<double, 6> v = {c.a, c.b, c.c, c.alpha, c.beta, c.gamma};
      // No e.s.d.s here: the cell's are mxi_scale's, from two-theta (Graeme).
      const std::array<double, 6> e = {0, 0, 0, 0, 0, 0};
      std::printf("Unit cell%s: %s %s %s A, %s %s %s deg; volume %.0f A^3\n",
                  experiments.size() > 1
                      ? (" [" + std::to_string(i) + "]").c_str()
                      : "",
                  with_esd(v[0], e[0]).c_str(), with_esd(v[1], e[1]).c_str(),
                  with_esd(v[2], e[2]).c_str(), with_esd(v[3], e[3], 3).c_str(),
                  with_esd(v[4], e[4], 3).c_str(),
                  with_esd(v[5], e[5], 3).c_str(), c.volume());
      // Where the centre of diffraction sits off the axis, if refined: its
      // distance, and its direction about the axis from the frame's first
      // perpendicular, as refinement parameterises it.
      const Goniometer &g = experiments[i].goniometer;
      if (g.has_offset()) {
        const Vec3 ax = g.axis.normalized();
        const Vec3 ref =
            std::fabs(ax.z) < 0.9 ? Vec3{0.0, 0.0, 1.0} : Vec3{1.0, 0.0, 0.0};
        const Vec3 pu = ax.cross(ref).normalized(), pv = ax.cross(pu);
        const Vec3 &c0 = g.sample_offset;
        std::printf("  the centre of diffraction %.1f um off the axis, at %.1f "
                    "deg\n",
                    1000.0 * c0.norm(),
                    std::atan2(c0.dot(pv), c0.dot(pu)) * 180.0 / M_PI);
      }
      json::Array jc, je, jlo, jhi;
      for (int k = 0; k < 6; ++k) {
        jc.push_back(json::Value(v[static_cast<std::size_t>(k)]));
        je.push_back(json::Value(e[static_cast<std::size_t>(k)]));
      }
      double turned = 0.0;
      if (x.A_points.size() > 1) {
        std::array<double, 6> lo, hi;
        lo.fill(HUGE_VAL);
        hi.fill(-HUGE_VAL);
        // The model at each degree of the scan (crystal_over_scan), not its
        // control points read as scan positions: those had drawn a 3600-image
        // sweep's cell over its first 37 images, and the ranges below from
        // coefficients the cell never takes in the interior.
        const std::vector<ScanSample> samples =
            crystal_over_scan(x, experiments[i].scan);
        json::Array s_angle, s_image, s_a, s_b, s_c, s_alpha, s_beta, s_gamma;
        for (const ScanSample &sample : samples) {
          Crystal at = x;
          at.A = sample.A;
          const UnitCell p = at.cell();
          s_angle.push_back(json::Value(sample.angle));
          s_image.push_back(json::Value(sample.image));
          s_a.push_back(json::Value(p.a));
          s_b.push_back(json::Value(p.b));
          s_c.push_back(json::Value(p.c));
          s_alpha.push_back(json::Value(p.alpha));
          s_beta.push_back(json::Value(p.beta));
          s_gamma.push_back(json::Value(p.gamma));
          const std::array<double, 6> q = {p.a,     p.b,    p.c,
                                           p.alpha, p.beta, p.gamma};
          for (std::size_t k = 0; k < 6; ++k) {
            lo[k] = std::min(lo[k], q[k]);
            hi[k] = std::max(hi[k], q[k]);
          }
        }
        json::Object scan_cell;
        scan_cell["experiment"] = json::Value(static_cast<long>(i));
        scan_cell["angle"] = json::Value(std::move(s_angle));
        scan_cell["image"] = json::Value(std::move(s_image));
        scan_cell["a"] = json::Value(std::move(s_a));
        scan_cell["b"] = json::Value(std::move(s_b));
        scan_cell["c"] = json::Value(std::move(s_c));
        scan_cell["alpha"] = json::Value(std::move(s_alpha));
        scan_cell["beta"] = json::Value(std::move(s_beta));
        scan_cell["gamma"] = json::Value(std::move(s_gamma));
        over_scan.push_back(json::Value(std::move(scan_cell)));
        // The orientation: A = U B, B the Cholesky factor of A^T A, so U =
        // A B^-1; the turn from the first scan point's U to the last's.
        const auto orientation = [](const Mat3 &A) {
          const Mat3 G = A.transpose() * A;
          Mat3 B{0, 0, 0, 0, 0, 0, 0, 0, 0};
          B(0, 0) = std::sqrt(G(0, 0));
          B(0, 1) = G(0, 1) / B(0, 0);
          B(0, 2) = G(0, 2) / B(0, 0);
          B(1, 1) = std::sqrt(G(1, 1) - B(0, 1) * B(0, 1));
          B(1, 2) = (G(1, 2) - B(0, 1) * B(0, 2)) / B(1, 1);
          B(2, 2) = std::sqrt(G(2, 2) - B(0, 2) * B(0, 2) - B(1, 2) * B(1, 2));
          return A * B.inverse();
        };
        const Mat3 R = orientation(samples.back().A) *
                       orientation(samples.front().A).transpose();
        const double cosine = (R(0, 0) + R(1, 1) + R(2, 2) - 1.0) / 2.0;
        turned =
            std::acos(std::fmax(-1.0, std::fmin(1.0, cosine))) * 180.0 / M_PI;
        std::printf("  over the scan: a %.4f to %.4f, b %.4f to %.4f, c %.4f "
                    "to %.4f A; alpha %.3f to %.3f, beta %.3f to %.3f, gamma "
                    "%.3f to %.3f deg; the orientation turns %.3f deg\n",
                    lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], lo[3], hi[3],
                    lo[4], hi[4], lo[5], hi[5], turned);
        // A crystal does not change shape like this within a sweep: radiation
        // damage and goniometer wobble stay well inside these. Said, as a
        // scale at zero is (a small molecule's cell edge had run 5.32 to
        // 6.51 A, its orientation turning 7.8 deg, in silence).
        double edge = 0.0, angle = 0.0;
        for (std::size_t k = 0; k < 3; ++k)
          edge = std::fmax(edge, (hi[k] - lo[k]) / lo[k] * 100.0);
        for (std::size_t k = 3; k < 6; ++k)
          angle = std::fmax(angle, hi[k] - lo[k]);
        if (edge > 2.0 || angle > 2.0 || turned > 2.0)
          std::printf("  WARNING: experiment %zu's crystal changes more than "
                      "a crystal does within a sweep -- a cell edge by %.1f "
                      "per cent, an angle by %.1f deg, the orientation by %.1f "
                      "deg; the scan-varying model is not held by its data "
                      "(--interval for fewer control points, or a stronger "
                      "--crystal-restraint)\n",
                      i, edge, angle, turned);
        for (std::size_t k = 0; k < 6; ++k) {
          jlo.push_back(json::Value(lo[k]));
          jhi.push_back(json::Value(hi[k]));
        }
      }
      t_exp.push_back(json::Value(static_cast<long>(i)));
      t_cell.push_back(json::Value(std::move(jc)));
      t_lo.push_back(json::Value(std::move(jlo)));
      t_hi.push_back(json::Value(std::move(jhi)));
      t_turn.push_back(json::Value(turned));
    }
    // The geometry as given and as refined.
    // Once a detector: lattices on one sweep share it.
    std::vector<GeometryFigures> said;
    for (std::size_t i = 0; i < experiments.size() && i < given_models.size();
         ++i) {
      const GeometryFigures before = geometry_figures(given_models[i]);
      const GeometryFigures after = geometry_figures(experiments[i]);
      bool seen = false;
      for (const GeometryFigures &g : said)
        seen = seen || (g.distance == after.distance &&
                        g.beam_x == after.beam_x && g.beam_y == after.beam_y);
      if (seen)
        continue;
      said.push_back(after);
      std::printf("Beam centre %.2f, %.2f px and distance %.3f mm as given; "
                  "refined %.2f, %.2f px and %.3f mm\n",
                  before.beam_x, before.beam_y, before.distance, after.beam_x,
                  after.beam_y, after.distance);
      if (i == 0) {
        summary.set("geometry_given", to_json(before));
        summary.set("geometry", to_json(after));
      }
    }
    summary.set("reflections",
                json::Value(static_cast<long>(reflections.nrows)));
    summary.set("refined_on", json::Value(static_cast<long>(result.n_used)));
    summary.set("rejected", json::Value(static_cast<long>(result.n_rejected)));
    summary.set("rmsd_px_px_images",
                json::Value(json::Array{json::Value(result.rmsd_x),
                                        json::Value(result.rmsd_y),
                                        json::Value(result.rmsd_z)}));
    summary.set("scan_varying", json::Value(two.varied ? 1L : 0L));
    if (!over_scan.empty())
      summary.set("cell_over_scan", json::Value(std::move(over_scan)));
    summary.table(
        "crystals",
        {"experiment", "cell", "scan_min", "scan_max", "orientation_turn_deg"},
        {t_exp, t_cell, t_lo, t_hi, t_turn});
    {
      std::vector<json::Value> p_pass, p_cycle, p_n, p_rej, p_x, p_y, p_z;
      const auto add = [&](const char *pass, const RefineResult &r) {
        for (std::size_t c = 0; c < r.cycles.size(); ++c) {
          const RefineCycle &k = r.cycles[c];
          p_pass.push_back(json::Value(std::string(pass)));
          p_cycle.push_back(json::Value(static_cast<long>(c + 1)));
          p_n.push_back(json::Value(static_cast<long>(k.reflections)));
          p_rej.push_back(json::Value(static_cast<long>(k.rejected)));
          p_x.push_back(json::Value(k.rmsd_x));
          p_y.push_back(json::Value(k.rmsd_y));
          p_z.push_back(json::Value(k.rmsd_z));
        }
      };
      add("scan-static", two.static_pass);
      if (two.varied)
        add("scan-varying", two.varying_pass);
      summary.table("cycles",
                    {"pass", "cycle", "reflections", "rejected", "rmsd_x_px",
                     "rmsd_y_px", "rmsd_z_images"},
                    {p_pass, p_cycle, p_n, p_rej, p_x, p_y, p_z});
    }
    const double t_write_start = now_wall();
    const ExperimentList written = whole_again(whole, views, experiments);
    Table written_table =
        rows_again(std::move(reflections), crystal_less, views, written);
    const json::Value document = experiments_to_json(written);
    const json::Value summary_document = summary.document();
    rflx::write_outputs(out, &document, &written_table, "mxi_refine",
                        &summary_document);
    t_write = now_wall() - t_write_start;
    std::printf("Wrote %s\n", rflx::describe(out).c_str());

    if (args.has("--timing")) {
      Timing timing(true, t_start);
      timing.add("reading", t_read);
      timing.add("building the target rows", g_observations_seconds);
      timing.add("the jacobian", g_jacobian_seconds);
      timing.add("the normal equations", g_normal_seconds);
      timing.add("the solve", g_solve_seconds);
      timing.add("the trial residuals", g_residual_seconds);
      timing.add("outlier rejection", g_outlier_seconds);
      timing.add("writing", t_write);
      timing.report(stdout);
    }
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi_refine: %s\n", e.what());
    return 1;
  }
}

} // namespace mxi

// The program's entry, run by mxi (refine) -- apps/mxi.cc.
int mxi_refine_main(int argc, char **argv) {
  // Mirrored to mxi_refine.log in the working directory, as DIALS writes
  // dials.<program>.log; not for a run that only asks for help.
  if (!mxi::only_asks_for_help(argc, argv))
    mxi::mirror_to_log("mxi_refine.log");
  return mxi::run_refine(argc, argv);
}
