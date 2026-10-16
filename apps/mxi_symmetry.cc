// mxi_symmetry: the Laue group and space group of one sweep's integrated
// data, and the data reindexed into it -- as dials.symmetry does. The lattice's
// metric symmetry by Le Page's method; each element and each subgroup scored
// as Evans (2011); screw axes from the absences.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <map>
#include <numeric>
#include <set>
#include <string>

#include "args.hh"
#include "cosym.hh"
#include "expt.hh"
#include "laue.hh"
#include "log_mirror.hh"
#include "parallel.hh"
#include "processing.hh"
#include "refl.hh"
#include "rflx.hh"
#include "symmetry.hh"
#include "timing.hh"
#include "views.hh"

namespace mxi {

namespace {

// An element as cctbx describes it: its order and its axis, "2 |(1,1,0)" --
// the direction in the cell's lattice basis that the rotation, acting on
// coordinates as h' = h R has it act on indices, leaves where it was: the
// cross product of two independent rows of R - 1, divided by their common
// factor, its first non-zero component positive.
std::string describe_element(const Rotation &r, int order) {
  long m[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      m[i][j] = r[static_cast<std::size_t>(i * 3 + j)] - (i == j ? 1 : 0);
  long axis[3] = {0, 0, 0};
  for (int a = 0; a < 3 && axis[0] == 0 && axis[1] == 0 && axis[2] == 0; ++a)
    for (int b = a + 1; b < 3; ++b) {
      const long u[3] = {m[a][1] * m[b][2] - m[a][2] * m[b][1],
                         m[a][2] * m[b][0] - m[a][0] * m[b][2],
                         m[a][0] * m[b][1] - m[a][1] * m[b][0]};
      if (u[0] != 0 || u[1] != 0 || u[2] != 0) {
        axis[0] = u[0];
        axis[1] = u[1];
        axis[2] = u[2];
        break;
      }
    }
  long g = 0;
  for (long v : axis)
    g = std::gcd(g, std::labs(v));
  if (g > 1)
    for (long &v : axis)
      v /= g;
  for (long v : axis)
    if (v != 0) {
      if (v < 0)
        for (long &w : axis)
          w = -w;
      break;
    }
  char text[64];
  std::snprintf(text, sizeof text, "%d |(%ld,%ld,%ld)", order, axis[0], axis[1],
                axis[2]);
  return text;
}

} // namespace

namespace {

// A rotation as what it does to indices: "h,k,l" to, say, "-k,h,l".
std::string as_hkl(const Rotation &g) {
  const char *names = "hkl";
  std::string out = "h,k,l to ";
  for (int i = 0; i < 3; ++i) {
    std::string term;
    for (int j = 0; j < 3; ++j) {
      const int v = g[3 * j + i]; // (R^T h)_i = sum_j R_ji h_j
      if (v == 0)
        continue;
      if (v < 0)
        term += "-";
      else if (!term.empty())
        term += "+";
      if (std::abs(v) != 1)
        term += std::to_string(std::abs(v));
      term += names[j];
    }
    out += (i ? "," : "") + (term.empty() ? std::string("0") : term);
  }
  return out;
}

void usage() {
  std::printf(
      "usage: mxi_symmetry [options] INTEGRATED_EXPT INTEGRATED_REFL\n"
      "  --no-cosym        several crystals pooled as they are, not first\n"
      "                    brought to one setting by cosym\n"
      "  --threads N       N threads for scoring and naming the subgroups; "
      "every\n"
      "                    core by default, the answer the same on any number\n"
      "  --max-delta D     the lattice's symmetry to D degrees of obliquity "
      "(2)\n"
      "  -o PATH           where to write: one .rflx (symmetrized.rflx), or "
      "with\n"
      "                    a .refl the DIALS pair, as --output-expt, -refl\n"
      "  --output-expt P   the DIALS pair's experiment list "
      "(symmetrized.expt)\n"
      "  --output-refl P   the DIALS pair's reflections (symmetrized.refl)\n"
      "  --timing          where the time goes\n");
}

} // namespace

int run_program(int argc, char **argv) {
  const std::set<std::string> known = {
      "--max-delta",        "--threads",     "--no-cosym",    "-o",
      "--same-orientation", "--output-expt", "--output-refl", "--timing"};
  std::set<std::string> takes_value = known;
  takes_value.erase("--timing");
  takes_value.erase("--no-cosym");
  Arguments args = parse_arguments(argc, argv, known, takes_value);
  // One .rflx stands for the experiment list and the reflections
  // (docs/rflx.md).
  args.positional = rflx::as_pair(args.positional);
  if (args.help) {
    usage();
    return 0;
  }
  if (!args.ok) {
    std::fprintf(stderr, "mxi_symmetry: %s\n", args.error.c_str());
    return 2;
  }
  set_parallel_threads(static_cast<std::size_t>(args.number("--threads", 0.0)));
  if (args.positional.size() != 2) {
    std::fprintf(stderr,
                 "mxi_symmetry: expected an .expt and a .refl, got %zu files\n",
                 args.positional.size());
    usage();
    return 2;
  }
  Timing timing(args.has("--timing"));
  try {
    double mark = Timing::now();
    const auto phase = [&](const char *name) {
      const double t = Timing::now();
      timing.add(name, t - mark);
      mark = t;
    };
    // The lattices; each sweep's own crystal-less experiment carried through
    // (docs/multi-lattice.md, stage 2).
    const ExperimentList whole = read_experiments(args.positional[0]);
    Table whole_table = read_reflections(args.positional[1]);
    const Lattices views = lattices_of(whole);
    ExperimentList experiments = views.list;
    // The crystal-less rows apart, the rest moved, not copied: one table
    // held, not three (Graeme: symmetry 22 GB, scaling 26, on 36000
    // images).
    const Table crystal_less = crystal_less_rows(whole_table, views);
    Table reflections = lattice_rows(std::move(whole_table), views);
    phase("reading");
    // Input in a centred setting -- data already through dials.symmetry or
    // mxi_symmetry, I 2 3 for instance -- taken to its primitive cell first:
    // the lattice's symmetry and everything after work in the basis given, as
    // a primitive one, and the conventional cell of a centred lattice is not
    // one. Taken as it was, a body-centred cubic cell was analysed as a
    // primitive cubic one -- every Patterson group P, P 2 3 recommended for
    // I 2 3 data (Graeme, a cluster of the 36 insulin crystals; item 70).
    // Reflections whose indices are not whole numbers in the primitive basis
    // are forbidden by the centring; any there are dropped, and counted.
    if (experiments.size() > 0 && experiments[0].crystal) {
      char lattice = 'P';
      for (char ch : experiments[0].crystal->space_group_hall)
        if (ch != ' ' && ch != '-') {
          lattice =
              static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
          break;
        }
      static const std::map<char, const char *> to_primitive = {
          {'I', "-1/2a+1/2b+1/2c,1/2a-1/2b+1/2c,1/2a+1/2b-1/2c"},
          {'F', "1/2b+1/2c,1/2a+1/2c,1/2a+1/2b"},
          {'C', "1/2a-1/2b,1/2a+1/2b,c"},
          {'A', "a,1/2b-1/2c,1/2b+1/2c"},
          {'B', "1/2a-1/2c,b,1/2a+1/2c"},
          {'R', "2/3a+1/3b+1/3c,-1/3a+1/3b+1/3c,-1/3a-2/3b+1/3c"}};
      const auto found = to_primitive.find(lattice);
      if (found != to_primitive.end()) {
        const ChangeOfBasis cb = ChangeOfBasis::parse(found->second);
        std::vector<std::size_t> keep;
        if (reflections.has("miller_index")) {
          const Column &hkl = reflections.at("miller_index");
          for (std::size_t r = 0; r < reflections.nrows; ++r) {
            bool whole = true;
            for (int j = 0; j < 3; ++j) {
              double v = 0.0;
              for (int i = 0; i < 3; ++i)
                v += cb.M(i, j) * static_cast<double>(hkl.ints[r * 3 + i]);
              whole = whole && std::fabs(v - std::round(v)) < 1e-6;
            }
            if (whole)
              keep.push_back(r);
          }
        }
        const std::size_t dropped = reflections.nrows - keep.size();
        if (dropped > 0)
          reflections = select_rows(reflections, keep);
        const std::string name =
            SpaceGroup::from_hall(experiments[0].crystal->space_group_hall)
                .name();
        reindex(experiments, reflections, cb, SpaceGroup::from_hall(" P 1"));
        std::printf("The input is in a centred setting, %s: taken to its "
                    "primitive cell first, by %s",
                    name.c_str(), found->second);
        if (dropped > 0)
          std::printf("; %zu reflections forbidden by the centring dropped",
                      dropped);
        std::printf("\n");
      }
    }
    // Several sweeps of one crystal pooled, as dials.symmetry takes them: the
    // lattice is the first's, so every sweep must have a crystal and the same
    // cell -- sweeps indexed apart may not, and may not even share a setting,
    // which joint indexing gives them. Pooled as they are, unscaled: sweeps on
    // much different scales would want each normalised first, not yet done.
    if (experiments.size() == 0)
      throw std::runtime_error("no experiments");
    for (std::size_t i = 0; i < experiments.size(); ++i)
      if (!experiments[i].crystal)
        throw std::runtime_error("experiment " + std::to_string(i) +
                                 " has no crystal: index it first");
    const UnitCell cell = experiments[0].crystal->cell();
    for (std::size_t i = 1; i < experiments.size(); ++i) {
      const UnitCell other = experiments[i].crystal->cell();
      const double worst = std::max(
          {std::fabs(other.a / cell.a - 1.0), std::fabs(other.b / cell.b - 1.0),
           std::fabs(other.c / cell.c - 1.0),
           std::fabs(other.alpha - cell.alpha) / 90.0,
           std::fabs(other.beta - cell.beta) / 90.0,
           std::fabs(other.gamma - cell.gamma) / 90.0});
      if (worst > 0.02)
        throw std::runtime_error(
            "experiment " + std::to_string(i) +
            "'s cell differs from the first's by " +
            std::to_string(static_cast<int>(std::round(100.0 * worst))) +
            " per cent: index the sweeps together, so that they share one "
            "cell and one setting");
    }
    if (experiments.size() > 1)
      std::printf("%zu sweeps, pooled\n", experiments.size());
    std::printf("Cell %.3f %.3f %.3f A, %.3f %.3f %.3f deg\n", cell.a, cell.b,
                cell.c, cell.alpha, cell.beta, cell.gamma);

    // Several crystals: one setting first (docs/multi-lattice.md). Sweeps of
    // one crystal share its setting; two crystals, or two lattices in one
    // sweep, are each in a basis of their own choosing among those their
    // lattice allows -- with I23 in an I432 lattice, half of them the other's
    // twin. So cosym (Gildea & Winter, 2018): each experiment's intensities
    // merged in P1 and normalised, correlated with every other's under every
    // rotation of the lattice, and each reindexed into the first's setting.
    // Crystals indexed together are in one setting already, whatever their
    // refinement did after: their orientations at the goniometer's datum
    // agree. Cosym is for crystals indexed apart. Run on four sweeps of one
    // small molecule crystal -- tetragonal, c/a 1.01, a lattice all but cubic
    // -- it reindexed one by a false cubic rotation, cosines 0.42 (Graeme).
    // So each crystal's U, from A = U B, B the Cholesky factor of A^T A,
    // against the first's: within --same-orientation degrees, all of them,
    // and there is nothing for cosym to do.
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
    const double same_within = args.number("--same-orientation", 5.0);
    double widest = 0.0;
    const Mat3 U0 = orientation(experiments[0].crystal->A);
    for (std::size_t i = 1; i < experiments.size(); ++i) {
      const Mat3 R = orientation(experiments[i].crystal->A) * U0.transpose();
      const double cosine = (R(0, 0) + R(1, 1) + R(2, 2) - 1.0) / 2.0;
      widest =
          std::max(widest, std::acos(std::fmax(-1.0, std::fmin(1.0, cosine))) *
                               180.0 / M_PI);
    }
    const bool one_crystal = widest <= same_within;
    if (experiments.size() > 1 && one_crystal && !args.has("--no-cosym"))
      std::printf("\nOne setting already: the %zu crystals' orientations agree "
                  "within %.3f deg, so no cosym\n",
                  experiments.size(), widest);
    if (experiments.size() > 1 && !one_crystal && !args.has("--no-cosym")) {
      std::vector<P1Intensities> sets;
      const Column &id = reflections.at("id");
      for (std::size_t j = 0; j < experiments.size(); ++j) {
        std::vector<std::size_t> rows;
        for (std::size_t r = 0; r < reflections.nrows; ++r)
          if (id.integer(r) == static_cast<std::int64_t>(j))
            rows.push_back(r);
        Table part = select_rows(reflections, rows);
        Column &own = part.modify_int_column("id", "int", 1);
        std::fill(own.ints.begin(), own.ints.end(), 0);
        ExperimentList one;
        one.experiments.push_back(experiments[j]);
        P1Selection ignored;
        P1Intensities d = merge_in_p1(one, part, &ignored);
        normalise(d, 200);
        sets.push_back(std::move(d));
      }
      const std::vector<Rotation> lattice =
          lattice_symmetry(cell, args.number("--max-delta", 2.0));
      const CosymResult c = cosym(sets, lattice);
      std::printf("\nOne setting for %zu crystals (cosym): the lattice's %zu "
                  "rotations; principal components",
                  experiments.size(), lattice.size());
      for (double share : c.explained)
        std::printf(" %.1f%%", 100.0 * share);
      std::printf(", %zu kept\n", c.dimensions);
      if (args.has("--timing"))
        std::printf(
            "  cosym: correlations %.3f s, eigen-decomposition of %zu x "
            "%zu %.3f s\n",
            c.seconds_correlations, experiments.size() * lattice.size(),
            experiments.size() * lattice.size(), c.seconds_eigen);
      for (std::size_t j = 0; j < experiments.size(); ++j) {
        const Rotation &g = c.operations[c.chosen[j]];
        std::printf("  experiment %zu: %s; cosine %.3f; CC with the first "
                    "%.3f over %zu reflections\n",
                    j, as_hkl(g).c_str(), c.cosine[j], c.cc[j], c.common[j]);
        if (c.chosen[j] == 0)
          continue;
        // Its indices h to R^T h, and its crystal A to A (R^T)^-1, so that
        // A h is unchanged -- every scan point too.
        Mat3 rt;
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            rt.m[3 * a + b] = g[3 * b + a];
        const Mat3 undo = rt.inverse();
        Crystal &crystal = *experiments[j].crystal;
        crystal.A = crystal.A * undo;
        for (Mat3 &at : crystal.A_points)
          at = at * undo;
        Column &miller = reflections.modify_int_column(
            "miller_index", "cctbx::miller::index<>", 3);
        for (std::size_t r = 0; r < reflections.nrows; ++r) {
          if (id.integer(r) != static_cast<std::int64_t>(j))
            continue;
          const Miller h = apply(g, {static_cast<int>(miller.ints[3 * r]),
                                     static_cast<int>(miller.ints[3 * r + 1]),
                                     static_cast<int>(miller.ints[3 * r + 2])});
          for (int k = 0; k < 3; ++k)
            miller.ints[3 * r + static_cast<std::size_t>(k)] = h[k];
        }
      }
      std::printf("\n");
    }

    P1Selection selection;
    const P1Intensities merged =
        merge_in_p1(experiments, reflections, &selection);
    phase("merging in P1, and the resolution limit");
    std::printf("%zu observations; %zu removed with I/sigma below -5\n",
                selection.observations, selection.negative);
    std::printf("Resolution from CC half above 0.6: %.2f A; from <I>/<sigma> "
                "above 4: %.2f A; "
                "the finer, %.2f A, used\n",
                selection.d_min_cc_half, selection.d_min_i_over_sigma,
                selection.d_min);
    P1Intensities normalised = merged;
    WilsonFit wilson_fit;
    const std::size_t wilson = normalise(normalised, 200, &wilson_fit);
    if (wilson_fit.fitted) {
      // B along each reciprocal axis, in A^2: the exponent's b times twice the
      // cell edge squared -- for a cell near orthogonal, what dials.symmetry
      // reports as B_cart's diagonal.
      std::printf(
          "Normalised by an anisotropic Wilson model, %zu reflections "
          "fitted in %d steps: B %.2f, %.2f, %.2f A^2 along a*, b*, c*\n",
          wilson_fit.used, wilson_fit.iterations,
          2.0 * wilson_fit.b[0] * cell.a * cell.a,
          2.0 * wilson_fit.b[1] * cell.b * cell.b,
          2.0 * wilson_fit.b[2] * cell.c * cell.c);
    } else {
      std::printf("Too few reflections for a Wilson model (%zu): normalised "
                  "in shells instead\n",
                  wilson_fit.used);
    }
    phase("normalising");
    std::printf("%zu Wilson outliers removed, E^2 of 16 or more\n", wilson);
    const std::vector<Rotation> lattice =
        lattice_symmetry(cell, args.number("--max-delta", 2.0));
    phase("the lattice's symmetry");
    const LaueScores scores = score_laue_groups(normalised, lattice);
    phase("scoring");
    timing.add("E(CC) and sigma(CC)", scores.t_estimates, 1);
    timing.add("the elements", scores.t_elements, 1);
    timing.add("the subgroups", scores.t_groups, 1);
    std::printf("%zu reflections merged in P1, Friedel mates apart; the "
                "lattice has %zu rotations, "
                "%zu symmetry "
                "elements, %zu subgroups\n",
                merged.size(), lattice.size(), scores.elements.size(),
                scores.groups.size());
    std::printf("E(CC; S) %.3f, from the intensities %.3f and from the "
                "identity %.3f; sigma(CC) "
                "factor %.3f\n",
                scores.cc_true, scores.e_cc_true, scores.cc_identity,
                scores.cc_sig_fac);

    std::printf("\nScoring each symmetry element\n");
    std::printf("  %-4s %-14s %7s %6s %7s %10s\n", "", "operator", "pairs",
                "CC", "Z-CC", "likelihood");
    for (std::size_t k = 0; k < scores.elements.size(); ++k) {
      const ElementScore &e = scores.elements[k];
      std::printf("  %-4c %-14s %7zu %6.3f %7.2f %10.3f %s\n",
                  static_cast<char>('a' + k),
                  describe_element(e.rotation, e.order).c_str(), e.pairs, e.cc,
                  e.z, e.likelihood, e.likelihood > 0.5 ? "***" : "");
    }

    std::printf("\nScoring all possible subgroups\n");
    std::printf("  %-14s %3s %10s %7s %6s %6s %5s %5s  %-20s %s\n",
                "Patterson group", "", "likelihood", "NetZcc", "Zcc+", "Zcc-",
                "CC", "CC-", "reindex", "elements");
    // Every subgroup's reference setting and name, each on its own thread --
    // naming them one after another was most of a small run -- then the table
    // printed in order.
    std::vector<std::optional<Setting>> settings(scores.groups.size());
    for_each_index(scores.groups.size(), [&](std::size_t k) {
      settings[k] = reference_setting(scores.groups[k].rotations, cell);
    });
    std::optional<Setting> best;
    for (std::size_t k = 0; k < scores.groups.size(); ++k) {
      const GroupScore &g = scores.groups[k];
      const std::optional<Setting> &s = settings[k];
      if (k == 0)
        best = s;
      std::string letters;
      for (std::size_t j = 0; j < g.contains.size(); ++j)
        letters += g.contains[j] ? static_cast<char>('a' + j) : '.';
      std::printf(
          "  %-14s %3s %10.3f %7.2f %6.2f %6.2f %5.2f %5.2f  %-20s %s\n",
          s ? s->group.name().c_str() : "?", k == 0 ? "***" : "", g.likelihood,
          g.z_net, g.z_for, g.z_against, g.cc_for, g.cc_against,
          s ? s->cb_text.c_str() : "?", letters.c_str());
    }
    if (!best)
      throw std::runtime_error(
          "the most likely group could not be put in a reference setting");
    phase("naming the subgroups, and printing");

    // Screw axes from the absences, in the chosen group's own setting.
    std::vector<Miller> hkl;
    for (const Miller &h : merged.hkl)
      hkl.push_back(best->cb.apply(h));
    const SpaceGroupChoice choice =
        choose_space_group(hkl, merged.i, merged.sigma, best->group);
    phase("the space group, by absences");
    std::printf(
        "\nSpace groups with Patterson group %s, judged by their absences\n",
        best->group.name().c_str());
    std::printf("  %-14s %8s %10s  %s\n", "space group", "tested", "<I/sigma>",
                "");
    for (const AbsenceTest &t : choice.candidates)
      std::printf("  %-14s %8zu %10.2f  %s\n", t.group.name().c_str(), t.tested,
                  t.mean_i_over_sigma,
                  t.consistent ? "consistent" : "inconsistent");
    std::printf("\nRecommended space group: %s\n",
                choice.chosen.name().c_str());
    for (const std::string &n : choice.indistinguishable)
      std::printf("  indistinguishable by its absences from %s\n", n.c_str());
    std::printf("Reindex operator: %s\n", best->cb_text.c_str());

    reindex(experiments, reflections, best->cb, choice.chosen);
    const UnitCell after = experiments[0].crystal->cell();
    std::printf("Cell %.3f %.3f %.3f A, %.3f %.3f %.3f deg\n", after.a, after.b,
                after.c, after.alpha, after.beta, after.gamma);
    const rflx::Outputs out = rflx::outputs(args, "symmetrized");
    const ExperimentList written = whole_again(whole, views, experiments);
    Table written_table =
        rows_again(std::move(reflections), crystal_less, views, written);
    const json::Value document = experiments_to_json(written);
    // The summary: each element, each subgroup, the choice.
    Summary summary("mxi_symmetry", args.positional[0]);
    {
      std::vector<json::Value> e_name, e_op, e_pairs, e_cc, e_z, e_l;
      for (std::size_t k = 0; k < scores.elements.size(); ++k) {
        const ElementScore &e = scores.elements[k];
        e_name.push_back(
            json::Value(std::string(1, static_cast<char>('a' + k))));
        e_op.push_back(json::Value(describe_element(e.rotation, e.order)));
        e_pairs.push_back(json::Value(static_cast<long>(e.pairs)));
        e_cc.push_back(json::Value(e.cc));
        e_z.push_back(json::Value(e.z));
        e_l.push_back(json::Value(e.likelihood));
      }
      summary.table(
          "elements",
          {"element", "operator", "pairs", "cc", "z_cc", "likelihood"},
          {e_name, e_op, e_pairs, e_cc, e_z, e_l});
      std::vector<json::Value> g_name, g_l, g_net, g_cc, g_against, g_cb;
      for (std::size_t k = 0; k < scores.groups.size(); ++k) {
        const GroupScore &g = scores.groups[k];
        const std::optional<Setting> &st = settings[k];
        g_name.push_back(json::Value(st ? st->group.name() : std::string("?")));
        g_l.push_back(json::Value(g.likelihood));
        g_net.push_back(json::Value(g.z_net));
        g_cc.push_back(json::Value(g.cc_for));
        g_against.push_back(json::Value(g.cc_against));
        g_cb.push_back(json::Value(st ? st->cb_text : std::string("?")));
      }
      summary.table("subgroups",
                    {"patterson_group", "likelihood", "net_z_cc", "cc_for",
                     "cc_against", "reindex"},
                    {g_name, g_l, g_net, g_cc, g_against, g_cb});
      summary.set("laue_group", json::Value(best->group.name()));
      summary.set("space_group", json::Value(choice.chosen.name()));
      json::Array alike;
      for (const std::string &n : choice.indistinguishable)
        alike.push_back(json::Value(n));
      summary.set("indistinguishable_by_absences",
                  json::Value(std::move(alike)));
      summary.set("reindex", json::Value(best->cb_text));
      summary.set("cell",
                  json::Value(json::Array{
                      json::Value(after.a), json::Value(after.b),
                      json::Value(after.c), json::Value(after.alpha),
                      json::Value(after.beta), json::Value(after.gamma)}));
    }
    const json::Value summary_document = summary.document();
    phase("reindexing");
    rflx::write_outputs(out, &document, &written_table, "mxi_symmetry",
                        &summary_document);
    std::printf("\nWrote %s\n", rflx::describe(out).c_str());
    phase("writing");
    timing.report(stdout);
  } catch (const std::exception &error) {
    std::fprintf(stderr, "mxi_symmetry: %s\n", error.what());
    return 1;
  }
  return 0;
}

} // namespace mxi

int main(int argc, char **argv) {
  // Mirrored to mxi_symmetry.log in the working directory, as DIALS writes
  // dials.symmetry.log; not for a run that only asks for help.
  if (!mxi::only_asks_for_help(argc, argv))
    mxi::mirror_to_log("mxi_symmetry.log");
  return mxi::run_program(argc, argv);
}
