// mxi: every program of the suite as a subcommand of one --
//
//   mxi import ../data/*.nxs
//   mxi find imported.rflx
//   mxi help
//
// One binary, not 24: each carried its own HDF5 and its own copy of the
// library, 156 MB of programs where one is some 8; one file to install; one
// name, which a program of the same name cannot overwrite (Graeme, 28
// October). Each program is its own source still, its main renamed at build
// time to mxi_<name>_main (CMakeLists.txt), and run here with the arguments
// after its name, argv[0] its name as typed -- "mxi find" -- for its messages.
// Called by another name -- mxi_find, a link to this in the build directory,
// for the tests -- the name is the subcommand.

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unistd.h>
#include <vector>

#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif

#define MXI_PROGRAM(name) int mxi_##name##_main(int argc, char **argv);
MXI_PROGRAM(import)
MXI_PROGRAM(find)
MXI_PROGRAM(index)
MXI_PROGRAM(refine)
MXI_PROGRAM(integrate)
MXI_PROGRAM(symmetry)
MXI_PROGRAM(scale)
MXI_PROGRAM(export)
MXI_PROGRAM(correlation)
MXI_PROGRAM(combine)
MXI_PROGRAM(split)
MXI_PROGRAM(convert)
MXI_PROGRAM(untrusted)
MXI_PROGRAM(show)
MXI_PROGRAM(report)
MXI_PROGRAM(max)
MXI_PROGRAM(profile)
MXI_PROGRAM(foreground)
MXI_PROGRAM(background)
MXI_PROGRAM(forward)
MXI_PROGRAM(grid)
MXI_PROGRAM(residuals)
MXI_PROGRAM(readtest)
MXI_PROGRAM(verify)
#undef MXI_PROGRAM

namespace {

struct Command {
  const char *name;
  int (*run)(int, char **);
  const char *what;
  bool diagnostic;
};

const Command kCommands[] = {
    {"import", mxi_import_main, "NXmx files to an experiment list", false},
    {"find", mxi_find_main, "spots, on the CPU or a GPU", false},
    {"index", mxi_index_main, "a lattice for the spots", false},
    {"refine", mxi_refine_main, "the models against the spots", false},
    {"integrate", mxi_integrate_main, "the reflections' intensities", false},
    {"symmetry", mxi_symmetry_main, "the Laue group and its setting", false},
    {"scale", mxi_scale_main, "scaled and merged intensities", false},
    {"export", mxi_export_main, "MTZ and the like", false},
    {"correlation", mxi_correlation_main, "crystals clustered, as cosym",
     false},
    {"combine", mxi_combine_main, "experiment lists joined", false},
    {"split", mxi_split_main, "an experiment list taken apart", false},
    {"convert", mxi_convert_main, ".rflx to and from .expt and .refl", false},
    {"untrusted", mxi_untrusted_main,
     "the detector's untrusted regions, with no window", false},
    {"view", nullptr, "the images, spots and regions, in napari", false},
    {"ui", nullptr, "the interface: the workflow, its images and its results",
     false},
    {"show", mxi_show_main, "what a file holds", false},
    {"report", mxi_report_main, "an HTML report of the processing", false},
    {"max", mxi_max_main, "a maximum projection of the images", false},
    {"verify", mxi_verify_main,
     "how much of what each lattice predicts was seen", true},
    {"profile", mxi_profile_main, "the profile model, estimated", true},
    {"foreground", mxi_foreground_main,
     "each reflection's foreground, for dials.image_viewer", true},
    {"background", mxi_background_main, "backgrounds fitted to given pixels",
     true},
    {"forward", mxi_forward_main, "spots rendered onto the pixels", true},
    {"grid", mxi_grid_main, "spots in Kabsch space", true},
    {"residuals", mxi_residuals_main, "centroid residuals, decomposed", true},
    {"readtest", mxi_readtest_main, "reading, timed", true},
};

void list(std::FILE *out) {
  std::fprintf(out,
               "mxi %s: processing rotation data, NXmx images to merged "
               "intensities\n\nusage: mxi COMMAND [ARGUMENTS]\n\n",
               MXI_VERSION);
  for (bool diagnostic : {false, true}) {
    std::fprintf(out,
                 diagnostic ? "\nFor looking at a step:\n" : "Commands:\n");
    for (const Command &c : kCommands)
      if (c.diagnostic == diagnostic)
        std::fprintf(out, "  %-12s %s\n", c.name, c.what);
  }
  std::fprintf(out,
               "\nEvery command takes -j N or --threads N: the threads for "
               "everything it\ndoes in parallel; 0, the default, every core "
               "this job may use (on Linux,\nits CPU affinity: a batch "
               "system's allocation, taskset).\n");
  std::fprintf(out, "\nmxi COMMAND --help for each; mxi --version.\n");
}

// Where this program is, symbolic links resolved: the system's answer --
// /proc/self/exe on Linux, _NSGetExecutablePath on macOS -- else argv[0], as a
// path where it has a slash and looked for on the PATH where it has none, as
// the shell found it. Taken from argv[0] alone, as at first, plain "mxi" run
// from the PATH looked in "." on macOS, which has no /proc, and mxi view
// could not find the viewer ("cannot run ./mxi_view", Graeme, 3 November).
std::string executable_path(const char *argv0) {
#if defined(__linux__)
  char buffer[PATH_MAX];
  const ssize_t n = ::readlink("/proc/self/exe", buffer, sizeof buffer - 1);
  if (n > 0) {
    buffer[n] = 0;
    return buffer;
  }
#elif defined(__APPLE__)
  char buffer[PATH_MAX];
  std::uint32_t size = sizeof buffer;
  if (_NSGetExecutablePath(buffer, &size) == 0) {
    char real[PATH_MAX];
    if (::realpath(buffer, real))
      return real;
    return buffer;
  }
#endif
  std::string found = argv0;
  if (found.find('/') == std::string::npos) {
    const char *path = std::getenv("PATH");
    std::string dirs = path ? path : "";
    std::size_t at = 0;
    while (at <= dirs.size()) {
      const std::size_t end = std::min(dirs.find(':', at), dirs.size());
      const std::string dir = end > at ? dirs.substr(at, end - at) : ".";
      const std::string candidate = dir + "/" + argv0;
      if (::access(candidate.c_str(), X_OK) == 0) {
        found = candidate;
        break;
      }
      at = end + 1;
    }
  }
  char real[PATH_MAX];
  return ::realpath(found.c_str(), real) ? std::string(real) : found;
}

// The commands that are Python, beside this binary: mxi view, mxi_view; mxi
// ui, the interface, mxi-ui -- both in the wheel mxi-python too.
const char *python_script(const char *command) {
  if (std::strcmp(command, "view") == 0)
    return "mxi_view";
  if (std::strcmp(command, "ui") == 0)
    return "mxi-ui";
  return nullptr;
}

// Runs the Python command beside this binary, its arguments passed on.
int python_command(const char *command, int argc, char **argv,
                   const char *self) {
  const std::string here = executable_path(self);
  const std::size_t slash = here.rfind('/');
  const std::string script =
      (slash == std::string::npos ? std::string(".") : here.substr(0, slash)) +
      "/" + python_script(command);
  std::vector<char *> args;
  args.push_back(const_cast<char *>(script.c_str()));
  for (int i = 1; i < argc; ++i)
    args.push_back(argv[i]);
  args.push_back(nullptr);
  ::execv(script.c_str(), args.data());
  std::fprintf(stderr, "mxi %s: cannot run %s: %s\n", command, script.c_str(),
               std::strerror(errno));
  return 2;
}

} // namespace

int main(int argc, char **argv) {
  // Called as mxi_find (a link): that is the command, argv as given.
  const char *base = std::strrchr(argv[0], '/');
  base = base ? base + 1 : argv[0];
  if (std::strncmp(base, "mxi_", 4) == 0) {
    for (const Command &c : kCommands)
      if (c.run && std::strcmp(base + 4, c.name) == 0)
        return c.run(argc, argv);
  }
  if (argc < 2 || std::strcmp(argv[1], "help") == 0 ||
      std::strcmp(argv[1], "--help") == 0 || std::strcmp(argv[1], "-h") == 0) {
    list(stdout);
    return 0;
  }
  if (std::strcmp(argv[1], "--version") == 0) {
    std::printf("mxi %s\n", MXI_VERSION);
    return 0;
  }
  for (const Command &c : kCommands) {
    if (std::strcmp(argv[1], c.name) != 0)
      continue;
    if (!c.run)
      return python_command(c.name, argc - 1, argv + 1, argv[0]);
    // argv[0] the command as typed, for its messages: "mxi find".
    std::string typed = std::string("mxi ") + c.name;
    std::vector<char *> args(argv + 1, argv + argc);
    args[0] = typed.data();
    args.push_back(nullptr);
    return c.run(argc - 1, args.data());
  }
  std::fprintf(stderr, "mxi: no command '%s'\n\n", argv[1]);
  list(stderr);
  return 2;
}
