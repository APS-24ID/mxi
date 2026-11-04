// mxi process: the whole sequence, import to export, a step at a time -- as
// a script of the eight commands would run it, each its own mxi command on
// the file the one before wrote, its own output and log. Options go where
// they belong: -j to every step, --gpu to spot finding and integration,
// --cell and --space-group to indexing (and the group to scaling),
// --d-min-auto to scaling. A step that fails stops it, said by name, its
// exit code and its log, with the command to run it again once put right
// (Graeme, 5 November: "given the simple script I typically use").
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
#include "self_path.hh"

extern char **environ;

namespace {

const char *const kUsage =
    "usage: mxi process MASTER.nxs... [options]\n"
    "  -j, --threads N   threads for everything this does in parallel; 0, "
    "every\n"
    "                    core this job may use (0)\n"
    "The whole sequence, a step at a time: mxi import, find, index, refine,\n"
    "integrate, symmetry, scale and export, each on the file the one before\n"
    "wrote, as a script of them would run them. Stops at a step that fails,\n"
    "saying which, its log, and the command to run it again.\n"
    "  -j N              to every step\n"
    "  --gpu             spot finding and integration on the GPU\n"
    "  --cell a,b,c,alpha,beta,gamma   the known cell: mxi index --known-cell\n"
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

} // namespace

// The program's entry, run by mxi (process) -- apps/mxi.cc.
int mxi_process_main(int argc, char **argv) {
  using namespace mxi;
  const std::set<std::string> known = {
      "--gpu",  "--cell", "--space-group", "--d-min-auto", "--dry-run",
      "--help", "-h"};
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
  for (Step &s : steps)
    s.args = with(s.args, threads);

  if (dry) {
    for (const Step &s : steps)
      std::printf("%s\n", shown(s).c_str());
    return 0;
  }

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
                     "-- and the steps after it.\n",
                     i + 1, n, s.name.c_str(), code, s.name.c_str(),
                     shown(s).c_str());
      return code < 0 ? 1 : code;
    }
  }
  std::printf("\nmxi process: all %zu steps done: the scaled data and their "
              "merging statistics in scaled.rflx (mxi show, mxi report), an "
              "unmerged MTZ in scaled.mtz.\n",
              n);
  return 0;
}
