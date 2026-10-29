// mxi_convert: between .rflx (docs/rflx.md) and the DIALS pair, .expt and
// .refl, its direction from what it is given.
//
//   mxi_convert indexed.expt indexed.refl     -> indexed.rflx
//   mxi_convert indexed.rflx                  -> indexed.expt, indexed.refl
//   mxi_convert indexed.rflx -o other.rflx    -> a .rflx rewritten
//
// The experiment list goes through as the JSON tree it is -- never made into
// mxi's models -- so what mxi does not model, profiles and scaling models and
// the rest, goes through untouched.

#include <cstdio>
#include <filesystem>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "args.hh"
#include "json.hh"
#include "log_mirror.hh"
#include "refl.hh"
#include "rflx.hh"
#include "views.hh"

namespace {

const char *kUsage =
    "usage: mxi convert [options] FILE...\n"
    "  -j, --threads N   threads for everything this does in parallel; 0, "
    "every\n"
    "                    core this job may use (0)\n"
    "  Between .rflx -- mxi's own files, one HDF5 file holding the experiment\n"
    "  list, the reflections or both (docs/rflx.md) -- and the DIALS pair.\n"
    "  Given a .rflx, writes its .expt and .refl, whichever it holds; given "
    "an\n"
    "  .expt, a .refl or both, writes a .rflx. A .rflx is told by what it is,\n"
    "  HDF5, not by its name.\n"
    "  -o, --output PATH      the .rflx to write (from the pair: the first\n"
    "                         input's name as .rflx); from a .rflx, a .rflx\n"
    "                         rewritten\n"
    "  --output-expt PATH     from a .rflx, the .expt (its name as .expt)\n"
    "  --output-refl PATH     from a .rflx, the .refl (its name as .refl)\n"
    "  --no-compress          a .rflx with no filter at all, for a reader\n"
    "                         without bitshuffle-LZ4 (HDF5 filter 32008),\n"
    "                         which numbers are stored with otherwise\n"
    "  --keep-sweeps          to the pair, keep each sweep's experiment\n"
    "                         with no crystal; by default, once any has a\n"
    "                         crystal, they are left out and their spots'\n"
    "                         id is -1, as DIALS has it\n";

std::string with_extension(const std::string &path, const std::string &ext) {
  return std::filesystem::path(path).replace_extension(ext).string();
}

bool ends_with(const std::string &s, const std::string &tail) {
  return s.size() >= tail.size() &&
         s.compare(s.size() - tail.size(), tail.size(), tail) == 0;
}

} // namespace

// The program's entry, run by mxi (convert) -- apps/mxi.cc.
int mxi_convert_main(int argc, char **argv) {
  using namespace mxi;
  if (!only_asks_for_help(argc, argv))
    mirror_to_log("mxi_convert.log");
  const std::set<std::string> known = {"-o",
                                       "--output",
                                       "--output-expt",
                                       "--output-refl",
                                       "--no-compress",
                                       "--keep-sweeps",
                                       "--drop-sweeps",
                                       "--keep-sweeps",
                                       "--help"};
  const std::set<std::string> takes_value = {"-o", "--output", "--output-expt",
                                             "--output-refl"};
  const Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (args.has("--help")) {
    std::printf("%s", kUsage);
    return 0;
  }
  if (!args.ok || args.positional.empty()) {
    if (!args.ok)
      std::fprintf(stderr, "mxi convert: %s\n", args.error.c_str());
    std::fprintf(stderr, "%s", kUsage);
    return 2;
  }
  const std::string output = args.value("-o", args.value("--output", ""));
  const bool compress = !args.has("--no-compress");
  rflx::set_uncompressed(!compress);
  try {
    std::vector<std::string> hdf5, other;
    for (const std::string &p : args.positional)
      (rflx::is_hdf5(p) ? hdf5 : other).push_back(p);

    if (!hdf5.empty()) {
      // From a .rflx.
      if (hdf5.size() > 1 || !other.empty()) {
        std::fprintf(stderr,
                     "mxi convert: one .rflx at a time, and nothing with it\n");
        return 2;
      }
      const std::string &in = hdf5.front();
      const std::optional<json::Value> experiments = rflx::read_experiments(in);
      std::vector<std::string> skipped;
      const std::optional<Table> reflections =
          rflx::read_reflections(in, &skipped);
      for (const std::string &name : skipped)
        std::printf("Left out the column %s, of a type mxi's tables have no "
                    "place for\n",
                    name.c_str());
      if (!experiments && !reflections) {
        std::fprintf(stderr,
                     "mxi convert: %s holds neither experiments nor "
                     "reflections\n",
                     in.c_str());
        return 1;
      }
      if (!output.empty() && ends_with(output, ".rflx")) {
        // The processing's summary goes with the file, as it was: a
        // conversion is not a step of the processing.
        const std::optional<json::Value> summary = rflx::read_summary(in);
        rflx::write(output, experiments ? &*experiments : nullptr,
                    reflections ? &*reflections : nullptr, "mxi_convert",
                    compress, summary ? &*summary : nullptr);
        std::printf("Wrote %s\n", output.c_str());
        return 0;
      }
      // The pair as DIALS has it: each sweep's crystal-less experiment left
      // out, the ids numbered to match, unindexed spots -1, as dials.index
      // leaves them -- the default (Graeme: with them in, the files would not
      // go straight into DIALS); --keep-sweeps keeps them. Before indexing,
      // when no experiment has a crystal, all are kept: there is nothing else.
      std::optional<Table> table = reflections;
      std::optional<json::Value> list = experiments;
      bool any_crystal = false;
      if (list)
        for (const json::Value &e : (*list)["experiment"].as_array())
          any_crystal = any_crystal || e["crystal"].is_number();
      if (list && any_crystal && !args.has("--keep-sweeps")) {
        const std::size_t before = (*list)["experiment"].as_array().size();
        list = dials_form(*list, table ? &*table : nullptr);
        const std::size_t after = (*list)["experiment"].as_array().size();
        if (after < before)
          std::printf("%zu crystal-less experiment%s left out, their spots' "
                      "id -1 (--keep-sweeps keeps them)\n",
                      before - after, before - after == 1 ? "" : "s");
      }
      // mxi's own scale model, mxi_physical, is not DIALS's physical model and
      // DIALS cannot read it: out of the .expt, the experiments' pointers too.
      if (list && (*list)["scaling_model"].is_array()) {
        bool mxi_model = false;
        for (const json::Value &m : (*list)["scaling_model"].as_array())
          mxi_model = mxi_model || (m["__id__"].is_string() &&
                                    m["__id__"].as_string() == "mxi_physical");
        if (mxi_model) {
          list->as_object()["scaling_model"] = json::Value(json::Array{});
          for (json::Value &e : list->as_object()["experiment"].as_array())
            e.as_object().erase("scaling_model");
          std::printf("mxi's scale model left out: DIALS does not read it\n");
        }
      }
      if (list) {
        const std::string out =
            args.value("--output-expt", with_extension(in, ".expt"));
        json::dump_file(out, *list);
        std::printf("Wrote the experiment list to %s\n", out.c_str());
      }
      if (table) {
        const Table &reflections_out = *table;
        const std::string out =
            args.value("--output-refl", with_extension(in, ".refl"));
        write_reflections(out, reflections_out);
        std::printf("Wrote %zu reflections to %s\n", reflections_out.nrows,
                    out.c_str());
      }
      return 0;
    }

    // From the DIALS pair: an .expt, a .refl or both, told by their contents
    // -- a .refl is msgpack, an .expt JSON -- with their names to break a tie.
    std::optional<json::Value> experiments;
    std::optional<Table> reflections;
    for (const std::string &p : other) {
      const bool looks_expt = ends_with(p, ".expt") || ends_with(p, ".json");
      const bool looks_refl = ends_with(p, ".refl");
      if (looks_expt || (!looks_refl && !experiments)) {
        if (experiments) {
          std::fprintf(stderr, "mxi convert: two experiment lists given\n");
          return 2;
        }
        experiments = json::parse_file(p);
      } else {
        if (reflections) {
          std::fprintf(stderr, "mxi convert: two reflection tables given\n");
          return 2;
        }
        reflections = read_reflections(p);
      }
    }
    const std::string out =
        output.empty() ? with_extension(other.front(), ".rflx") : output;
    rflx::write(out, experiments ? &*experiments : nullptr,
                reflections ? &*reflections : nullptr, "mxi_convert", compress);
    std::printf("Wrote %s%s%s\n", out.c_str(),
                experiments ? ": the experiment list" : ":",
                reflections ? (", " + std::to_string(reflections->nrows) +
                               " reflections")
                                  .c_str()
                            : "");
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mxi convert: %s\n", e.what());
    return 1;
  }
}
