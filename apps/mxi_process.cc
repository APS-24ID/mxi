// mxi process: the whole sequence, import to export, a step at a time -- as
// a script of the eight commands would run it, each its own mxi command on
// the file the one before wrote, its own output and log. Options go where
// they belong: -j to every step, --gpu to spot finding and integration,
// --no-joint to indexing,
// --cell and --space-group to indexing (and the group to scaling),
// --d-min-auto to scaling. A step that fails stops it, said by name, its
// exit code and its log, with the command to run it again once put right
// (Graeme, 5 November: "given the simple script I typically use").
//
// Given a .rflx instead of images, it picks up after the step that wrote it
// and finishes: a cluster mxi correlation wrote is scaled and exported, each
// cluster as the tutorial's cows, pigs and people do by hand (Graeme, 5
// November). The file's own writer, its root attribute "creator", says
// which step that was -- not its summary's last step: a cluster carries the
// summary it was cut from, scaling's last, though it needs scaling again.
//
// Each step is a process of its own -- this program again, by its own path
// (self_path.hh) -- so that nothing differs from running them by hand, and no
// setting of one step's (its threads) is another's.

#include <cstdio>
#include <cstring>
#include <set>
#include <spawn.h>
#include <string>
#include <sys/wait.h>
#include <vector>

#include "args.hh"
#include "log_mirror.hh"
#include "rflx.hh"
#include "self_path.hh"

extern char **environ;

namespace {

const char *const kUsage =
    "usage: mxi process MASTER.nxs... [options]\n"
    "       mxi process FILE.rflx [options]\n"
    "  -j, --threads N   threads for everything this does in parallel; 0, "
    "every\n"
    "                    core this job may use (0)\n"
    "The whole sequence, a step at a time: mxi import, find, index, refine,\n"
    "integrate, symmetry, scale and export, each on the file the one before\n"
    "wrote, as a script of them would run them. Stops at a step that fails,\n"
    "saying which, its log, and the command to run it again.\n"
    "Or given one .rflx, it picks up after the step that wrote it and\n"
    "finishes: a cluster from mxi correlation --write-clusters is scaled and\n"
    "exported; strong.rflx is indexed, and so on.\n"
    "  -j N              to every step\n"
    "  --gpu             spot finding and integration on the GPU\n"
    "  --cell a,b,c,alpha,beta,gamma   the known cell: mxi index --known-cell\n"
    "  --no-joint        several sweeps, each a crystal indexed on its own:\n"
    "                    mxi index --no-joint\n"
    "  --space-group NAME  the known space group: mxi scale --space-group, "
    "and\n"
    "                    with --cell, mxi index --known-space-group\n"
    "  --d-min-auto      mxi scale --d-min-auto: the limit from CC half\n"
    "  --dry-run         the commands, printed; nothing run\n";

struct Step {
  std::string name;              // the mxi command
  std::vector<std::string> args; // after the command: inputs, then options
  std::string output;            // what it writes
};

std::string shown(const Step &s) {
  std::string out = "mxi " + s.name;
  for (const std::string &a : s.args) {
    const bool plain =
        a.find_first_of(" \t'\"$\\") == std::string::npos && !a.empty();
    out += " " + (plain ? a : "'" + a + "'");
  }
  return out;
}

// Runs this program as `mxi STEP ARGS...`, its output to ours; its exit code,
// or 128 + the signal that ended it, or -1 where it could not be started.
int run(const std::string &self, const Step &s) {
  std::vector<std::string> words = {self, s.name};
  words.insert(words.end(), s.args.begin(), s.args.end());
  std::vector<char *> argv;
  for (std::string &w : words)
    argv.push_back(w.data());
  argv.push_back(nullptr);
  std::fflush(stdout);
  std::fflush(stderr);
  pid_t pid = 0;
  if (posix_spawn(&pid, self.c_str(), nullptr, nullptr, argv.data(), environ) !=
      0)
    return -1;
  int status = 0;
  if (waitpid(pid, &status, 0) < 0)
    return -1;
  if (WIFEXITED(status))
    return WEXITSTATUS(status);
  if (WIFSIGNALED(status))
    return 128 + WTERMSIG(status);
  return -1;
}

// The step after the one that wrote a .rflx: by its writer, "creator"; where
// that says nothing, by its summary's last step. "" where neither tells.
std::string next_after(const std::string &path, std::string *by) {
  static const std::vector<std::pair<std::string, std::string>> next = {
      {"mxi_import", "find"},
      {"mxi_find", "index"},
      {"mxi_index", "refine"},
      {"mxi_refine", "integrate"},
      {"mxi_integrate", "symmetry"},
      {"mxi_symmetry", "scale"},
      {"mxi_scale", "export"},
      // A cluster: a subset of the crystals, to be scaled afresh.
      {"mxi_correlation", "scale"}};
  const auto find = [&](const std::string &program) -> std::string {
    for (const auto &[from, step] : next)
      if (from == program)
        return step;
    return "";
  };
  const std::string creator = mxi::rflx::read_creator(path);
  if (const std::string step = find(creator); !step.empty()) {
    *by = creator;
    return step;
  }
  if (const auto summary = mxi::rflx::read_summary(path)) {
    const mxi::json::Value &steps = (*summary)["steps"];
    if (steps.is_array() && !steps.as_array().empty()) {
      const mxi::json::Value &last = steps.as_array().back()["program"];
      if (last.is_string())
        if (const std::string step = find(last.as_string()); !step.empty()) {
          *by = last.as_string();
          return step;
        }
    }
  }
  *by = creator;
  return "";
}

std::string join(const std::vector<std::string> &words) {
  std::string out;
  for (const std::string &w : words)
    out += (out.empty() ? "" : " ") + w;
  return out;
}

std::string ends(const std::string &s, std::size_t n) {
  return s.size() >= n ? s.substr(s.size() - n) : "";
}

} // namespace

// The program's entry, run by mxi (process) -- apps/mxi.cc.
int mxi_process_main(int argc, char **argv) {
  using namespace mxi;
  const std::set<std::string> known = {
      "--gpu",     "--cell",     "--space-group", "--d-min-auto",
      "--dry-run", "--no-joint", "--help",        "-h"};
  const std::set<std::string> takes_value = {"--cell", "--space-group"};
  Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (!args.ok) {
    std::fprintf(stderr, "mxi process: %s\n", args.error.c_str());
    std::fputs(kUsage, stderr);
    return 2;
  }
  if (args.has("--help") || args.has("-h")) {
    std::fputs(kUsage, stdout);
    return 0;
  }
  if (args.positional.empty()) {
    std::fprintf(stderr,
                 "mxi process: the images' master files, one a sweep\n");
    std::fputs(kUsage, stderr);
    return 2;
  }
  const bool dry = args.has("--dry-run");
  if (!dry && !only_asks_for_help(argc, argv))
    mirror_to_log("mxi_process.log");

  std::vector<std::string> threads;
  if (args.has("-j"))
    threads = {"-j", args.value("-j", "0")};
  const bool gpu = args.has("--gpu");
  const std::string cell = args.value("--cell", "");
  const std::string group = args.value("--space-group", "");

  const auto with = [](std::vector<std::string> a,
                       const std::vector<std::string> &b) {
    a.insert(a.end(), b.begin(), b.end());
    return a;
  };
  std::vector<Step> steps;
  steps.push_back({"import", args.positional, "imported.rflx"});
  steps.push_back({"find", {"imported.rflx"}, "strong.rflx"});
  if (gpu)
    steps.back().args.push_back("--gpu");
  steps.push_back({"index", {"strong.rflx"}, "indexed.rflx"});
  if (!cell.empty()) {
    steps.back().args.insert(steps.back().args.end(), {"--known-cell", cell});
    if (!group.empty())
      steps.back().args.insert(steps.back().args.end(),
                               {"--known-space-group", group});
  }
  // Several crystals, a sweep each, indexed apart: as the tutorial's cows,
  // pigs and people, the steps after running as for one.
  if (args.has("--no-joint"))
    steps.back().args.push_back("--no-joint");
  steps.push_back({"refine", {"indexed.rflx"}, "refined.rflx"});
  steps.push_back({"integrate", {"refined.rflx"}, "integrated.rflx"});
  if (gpu)
    steps.back().args.push_back("--gpu");
  steps.push_back({"symmetry", {"integrated.rflx"}, "symmetrized.rflx"});
  steps.push_back({"scale", {"symmetrized.rflx"}, "scaled.rflx"});
  if (!group.empty())
    steps.back().args.insert(steps.back().args.end(), {"--space-group", group});
  if (args.has("--d-min-auto"))
    steps.back().args.push_back("--d-min-auto");
  steps.push_back({"export", {"scaled.rflx"}, "scaled.mtz"});

  // One .rflx: from the step after the one that wrote it, that file its
  // input; the options of the steps done before said to be unused.
  const bool resuming =
      args.positional.size() == 1 && ends(args.positional[0], 5) == ".rflx";
  if (resuming) {
    const std::string &file = args.positional[0];
    std::string by;
    std::string first;
    try {
      first = next_after(file, &by);
    } catch (const std::exception &e) {
      std::fprintf(stderr, "mxi process: %s\n", e.what());
      return 1;
    }
    if (first.empty()) {
      std::fprintf(stderr,
                   "mxi process: %s says neither the step that wrote it nor "
                   "one to follow%s: run the next step by hand\n",
                   file.c_str(),
                   by.empty() ? "" : (" (written by " + by + ")").c_str());
      return 1;
    }
    std::size_t k = 0;
    while (steps[k].name != first)
      ++k;
    std::set<std::string> later;
    for (std::size_t i = k; i < steps.size(); ++i)
      for (const std::string &a : steps[i].args)
        later.insert(a);
    // Said by the options as given -- --cell, not the --known-cell indexing
    // takes -- and --space-group only where scaling, which takes it too, is
    // done as well.
    std::set<std::string> unused;
    for (std::size_t i = 0; i < k; ++i)
      for (const std::string &a : steps[i].args) {
        if (a.rfind("--", 0) != 0 || later.count(a))
          continue;
        if (a == "--known-cell")
          unused.insert("--cell");
        else if (a == "--known-space-group") {
          if (!later.count("--space-group"))
            unused.insert("--space-group");
        } else
          unused.insert(a);
      }
    steps.erase(steps.begin(), steps.begin() + static_cast<long>(k));
    steps.front().args.front() = file;
    // To stderr in a dry run, so that its output is the commands alone.
    FILE *note = dry ? stderr : stdout;
    std::fprintf(note, "mxi process: %s, written by %s: from mxi %s on\n",
                 file.c_str(), by.c_str(), first.c_str());
    for (const std::string &u : unused)
      std::fprintf(note, "mxi process: %s not used: its step is done already\n",
                   u.c_str());
  }
  for (Step &s : steps)
    s.args = with(s.args, threads);

  if (dry) {
    for (const Step &s : steps)
      std::printf("%s\n", shown(s).c_str());
    return 0;
  }

  // The options as given, to say again with mxi process.
  std::string shared;
  if (gpu)
    shared += " --gpu";
  if (!cell.empty())
    shared += " --cell " + cell;
  if (!group.empty())
    shared += " --space-group " + group;
  if (args.has("--no-joint"))
    shared += " --no-joint";
  if (args.has("--d-min-auto"))
    shared += " --d-min-auto";
  if (!threads.empty())
    shared += " -j " + threads[1];
  const std::string self = executable_path(argv[0]);
  const std::size_t n = steps.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Step &s = steps[i];
    std::printf("\n== mxi process: step %zu of %zu: %s\n\n", i + 1, n,
                shown(s).c_str());
    const int code = run(self, s);
    if (code != 0) {
      std::printf("\n");
      if (code < 0)
        std::fprintf(
            stderr,
            "mxi process: stopped at step %zu of %zu, mxi %s: it could "
            "not be started (%s)\n",
            i + 1, n, s.name.c_str(), self.c_str());
      else
        std::fprintf(stderr,
                     "mxi process: stopped at step %zu of %zu, mxi %s, which "
                     "failed (exit %d): its output is above, and in "
                     "mxi_%s.log. Once it is put right, run it again --\n"
                     "  %s\n"
                     "-- and the steps after it, or all that remains:\n"
                     "  mxi process %s%s\n",
                     i + 1, n, s.name.c_str(), code, s.name.c_str(),
                     shown(s).c_str(),
                     s.name == "import" ? join(args.positional).c_str()
                                        : s.args.front().c_str(),
                     shared.c_str());
      return code < 0 ? 1 : code;
    }
  }
  std::printf("\nmxi process: all %zu steps done: the scaled data and their "
              "merging statistics in scaled.rflx (mxi show, mxi report), an "
              "unmerged MTZ in scaled.mtz.\n",
              n);
  return 0;
}
