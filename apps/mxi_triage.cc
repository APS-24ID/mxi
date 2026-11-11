// mxi triage: what the scaled data say of twinning, as phenix.xtriage says
// it -- a start (src/triage.hh, docs/triage.md). The scaled observations of
// every crystal merged, as mxi correlation merges them, in their space group;
// the statistics that do not depend on twin laws, and those that do; and a
// summary that calls the data twinned only where both kinds say so.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "args.hh"
#include "expt.hh"
#include "log_mirror.hh"
#include "refl.hh"
#include "rflx.hh"
#include "scale.hh"
#include "symmetry.hh"
#include "triage.hh"

namespace mxi {

namespace {

const char *const kUsage =
    "usage: mxi triage SCALED.rflx [options]\n"
    "  -j, --threads N   threads for everything this does in parallel; 0, "
    "every\n"
    "                    core this job may use (0)\n"
    "What the scaled intensities say of twinning, as phenix.xtriage says it:\n"
    "their moments, N(Z) and the L test, independent of twin laws; the twin\n"
    "laws the lattice allows, each one's R and H-test fraction; and a summary\n"
    "that calls the data twinned only where both kinds say so.\n"
    "  --d-max D         the twin tests' low resolution (10 A, as xtriage's)\n"
    "  --d-min D         their high resolution (as xtriage chooses it: the\n"
    "                    finest d to which I/sigma >= 3 in 85 per cent of\n"
    "                    each shell, never coarser than 3.5 A)\n"
    "  --max-delta D     twin laws within D degrees' obliquity (3)\n"
    "  --min-partiality P  observations kept, as mxi correlation's (0.4)\n";

// Every crystal's scaled observations merged, in `group`: scaled, not an
// outlier, partiality at least min_partiality, as mxi correlation takes them
// (src/correlation.cc); merged as cctbx's merge_equivalents, the inverse-
// variance weighted mean, its sigma sqrt(max(V_w / n, 1 / sum w)).
TriageData merged(const Table &t, const SpaceGroup &group,
                  double min_partiality) {
  TriageData out;
  for (const char *c : {"miller_index", "flags", "intensity.scale.value",
                        "intensity.scale.variance"})
    if (!t.has(c))
      throw std::runtime_error(std::string("no ") + c + ": not a scaled file");
  const Column &hkl = t.at("miller_index"), &flags = t.at("flags"),
               &I = t.at("intensity.scale.value"),
               &V = t.at("intensity.scale.variance");
  const Column *partiality =
      t.has("partiality") ? &t.at("partiality") : nullptr;
  std::map<Miller, std::vector<std::pair<double, double>>> groups;
  for (std::size_t r = 0; r < t.nrows; ++r) {
    const std::int64_t f = flags.integer(r);
    if ((f & flag::kScaled) == 0 || (f & flag::kOutlierInScaling) != 0)
      continue;
    if (partiality && !(partiality->reals[r] >= min_partiality))
      continue;
    const double v = V.reals[r];
    if (!(v > 0.0) || !std::isfinite(I.reals[r]))
      continue;
    const Miller m{static_cast<std::int32_t>(hkl.ints[r * 3]),
                   static_cast<std::int32_t>(hkl.ints[r * 3 + 1]),
                   static_cast<std::int32_t>(hkl.ints[r * 3 + 2])};
    groups[group.unique(m)].emplace_back(I.reals[r], std::sqrt(v));
  }
  for (const auto &[m, obs] : groups) {
    double sw = 0, swx = 0, sw2 = 0;
    for (const auto &[x, s] : obs) {
      const double w = 1.0 / (s * s);
      sw += w;
      swx += w * x;
      sw2 += w * w;
    }
    const double mean = swx / sw;
    double sigma = std::sqrt(1.0 / sw);
    if (obs.size() > 1) {
      double swd2 = 0;
      for (const auto &[x, s] : obs)
        swd2 += (x - mean) * (x - mean) / (s * s);
      const double wvar = sw * sw > sw2 ? sw / (sw * sw - sw2) * swd2 : 0.0;
      sigma =
          std::sqrt(std::max(wvar / static_cast<double>(obs.size()), 1.0 / sw));
    }
    out.hkl.push_back(m);
    out.intensity.push_back(mean);
    out.sigma.push_back(sigma);
  }
  return out;
}

void moments_row(const char *name, double a, double a_un, double a_tw, double c,
                 double c_un, double c_tw) {
  std::printf("  %-14s %8.3f %10.3f %9.3f   %8.3f %10.3f %9.3f\n", name, a,
              a_un, a_tw, c, c_un, c_tw);
}

} // namespace

int run_triage(int argc, char **argv) {
  const std::set<std::string> known = {
      "--help", "-h", "--d-max", "--d-min", "--max-delta", "--min-partiality"};
  const std::set<std::string> takes_value = {"--d-max", "--d-min",
                                             "--max-delta", "--min-partiality"};
  Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (!args.ok) {
    std::fprintf(stderr, "mxi triage: %s\n", args.error.c_str());
    std::fputs(kUsage, stderr);
    return 2;
  }
  if (args.has("--help") || args.has("-h")) {
    std::fputs(kUsage, stdout);
    return 0;
  }
  args.positional = rflx::as_pair(args.positional);
  if (args.positional.size() != 2) {
    std::fputs(kUsage, stderr);
    return 2;
  }
  if (!only_asks_for_help(argc, argv))
    mirror_to_log("mxi_triage.log");
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
      std::fprintf(stderr, "mxi triage: no experiment has a crystal\n");
      return 1;
    }
    const SpaceGroup group =
        SpaceGroup::from_hall(first->crystal->space_group_hall);
    Crystal regular = *first->crystal;
    regularise_cell(regular, group);
    const UnitCell cell = regular.cell();
    const TriageData data =
        merged(reflections, group, args.number("--min-partiality", 0.4));
    TriageOptions options;
    options.d_max = args.number("--d-max", options.d_max);
    options.d_min = args.number("--d-min", 0.0);
    options.max_delta = args.number("--max-delta", options.max_delta);
    const Triage t = triage(data, cell, group, options);

    std::printf("%s: %zu unique reflections merged, in %s\n",
                args.positional[0].c_str(), data.hkl.size(),
                group.name().c_str());
    std::printf("  cell %.3f %.3f %.3f %.3f %.3f %.3f\n\n", cell.a, cell.b,
                cell.c, cell.alpha, cell.beta, cell.gamma);
    std::printf("The twin tests: %.2f to %.2f A, %zu reflections (%zu "
                "acentric, %zu centric)\n",
                t.d_max, t.d_min, t.used, t.acentric, t.centric);
    std::printf(
        "  as xtriage's: %s, %zu Wilson outlier%s removed\n\n",
        t.sharpened
            ? "on an absolute scale with the anisotropic fall-off removed"
            : "as merged (too few reflections to fit the fall-off)",
        t.outliers, t.outliers == 1 ? "" : "s");
    std::printf("Moments of the normalised intensities\n");
    std::printf("  %-14s %8s %10s %9s   %8s %10s %9s\n", "", "acentric",
                "untwinned", "perfect", "centric", "untwinned", "perfect");
    const Moments &a = t.acentric_moments, &c = t.centric_moments;
    moments_row("<I^2>/<I>^2", a.i2_over_i_squared, 2.000, 1.500,
                c.i2_over_i_squared, 3.000, 2.000);
    moments_row("<F>^2/<F^2>", a.f_squared_over_f2, 0.785, 0.885,
                c.f_squared_over_f2, 0.637, 0.785);
    moments_row("<|E^2-1|>", a.e2_minus_one, 0.736, 0.541, c.e2_minus_one,
                0.968, 0.736);

    std::printf("\nThe cumulative intensity distribution, N(z)\n");
    std::printf("  %4s %9s %9s %9s %9s\n", "z", "acentric", "Wilson", "centric",
                "Wilson");
    for (std::size_t k = 0; k < t.nz.z.size(); ++k)
      std::printf("  %4.1f %9.3f %9.3f %9.3f %9.3f\n", t.nz.z[k],
                  t.nz.acentric[k], t.nz.acentric_theory[k], t.nz.centric[k],
                  t.nz.centric_theory[k]);
    std::printf("  largest difference: acentric %.3f, centric %.3f\n",
                t.nz.max_deviation_acentric, t.nz.max_deviation_centric);

    std::printf("\nThe L test (Padilla & Yeates 2003): %zu pairs of acentric "
                "neighbours\n",
                t.l_test.pairs);
    std::printf("  <|L|> %.3f (untwinned 0.500, perfect twin 0.375)\n",
                t.l_test.mean_abs_l);
    std::printf("  <L^2> %.3f (untwinned 0.333, perfect twin 0.200)\n",
                t.l_test.mean_l2);
    std::printf("  multivariate Z %.2f (over 3.5, abnormal)\n",
                t.l_test.z_score);

    std::printf("\nTwin laws the lattice allows, within %.1f degrees\n",
                options.max_delta);
    if (t.twin_laws.empty())
      std::printf("  none\n");
    else {
      std::printf("  %-14s %-4s %10s %7s %7s %9s\n", "operator", "type",
                  "R metric", "pairs", "R obs", "H alpha");
      for (const TwinLaw &law : t.twin_laws)
        std::printf("  %-14s %-4s %10.5f %7zu %7.3f %9.3f\n", law.text.c_str(),
                    law.merohedral ? "M" : "PM", law.r_metric, law.pairs,
                    law.r_obs, law.h_alpha);
      std::printf("  (M merohedral, PM pseudo-merohedral; H alpha the H test's "
                  "twin fraction)\n");
    }
    // The summary wrapped to the page, two spaces in.
    std::printf("\nSummary\n");
    std::string line;
    std::size_t at = 0;
    while (at < t.verdict.size()) {
      const std::size_t space = t.verdict.find(' ', at);
      const std::string word = t.verdict.substr(
          at, space == std::string::npos ? std::string::npos : space - at);
      if (!line.empty() && line.size() + 1 + word.size() > 76) {
        std::printf("  %s\n", line.c_str());
        line.clear();
      }
      line += (line.empty() ? "" : " ") + word;
      at = space == std::string::npos ? t.verdict.size() : space + 1;
    }
    if (!line.empty())
      std::printf("  %s\n", line.c_str());
    return 0;
  } catch (const std::exception &error) {
    std::fprintf(stderr, "mxi triage: %s\n", error.what());
    return 1;
  }
}

} // namespace mxi

// The program's entry, run by mxi (triage) -- apps/mxi.cc.
int mxi_triage_main(int argc, char **argv) {
  return mxi::run_triage(argc, argv);
}
