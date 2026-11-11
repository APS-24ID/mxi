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

#include "output.hh"
#include "self_path.hh"

#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif

#define MXI_PROGRAM(name) std::int32_t mxi_##name##_main(int argc, char **argv);
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
MXI_PROGRAM(process)
MXI_PROGRAM(triage)
#undef MXI_PROGRAM

namespace {

struct Command {
  const char *name;
  std::int32_t (*run)(std::int32_t, char **);
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
    {"process", mxi_process_main, "every step above in turn, import to export",
     false},
    {"correlation", mxi_correlation_main, "crystals clustered, as cosym",
     false},
    {"triage", mxi_triage_main,
     "twinning and the intensities' statistics, as xtriage", false},
    {"combine", mxi_combine_main, "experiment lists joined", false},
    {"split", mxi_split_main, "an experiment list taken apart", false},
    {"convert", mxi_convert_main, ".rflx to and from .expt and .refl", false},
    {"untrusted", mxi_untrusted_main,
     "the detector's untrusted regions, with no window", false},
    {"view", nullptr, "the images, spots and regions: mxi ui's viewer", false},
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

void list(mxi::Stream out) {
  mxi::print_to(out,
                "mxi %s: processing rotation data, NXmx images to merged "
                "intensities\n\nusage: mxi COMMAND [ARGUMENTS]\n\n",
                MXI_VERSION);
  for (bool diagnostic : {false, true}) {
    mxi::print_to(out,
                  diagnostic ? "\nFor looking at a step:\n" : "Commands:\n");
    for (const Command &c : kCommands)
      if (c.diagnostic == diagnostic)
        mxi::print_to(out, "  %-12s %s\n", c.name, c.what);
  }
  mxi::print_to(out,
                "\nEvery command takes -j N or --threads N: the threads for "
                "everything it\ndoes in parallel; 0, the default, every core "
                "this job may use (on Linux,\nits CPU affinity: a batch "
                "system's allocation, taskset).\n");
  mxi::print_to(out, "\nmxi COMMAND --help for each; mxi --version.\n");
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
std::int32_t python_command(const char *command, int argc, char **argv,
                            const char *self) {
  const std::string here = mxi::executable_path(self);
  const std::size_t slash = here.rfind('/');
  const std::string script =
      (slash == std::string::npos ? std::string(".") : here.substr(0, slash)) +
      "/" + python_script(command);
  std::vector<char *> args;
  args.push_back(const_cast<char *>(script.c_str()));
  for (std::int32_t i = 1; i < argc; ++i)
    args.push_back(argv[i]);
  args.push_back(nullptr);
  ::execv(script.c_str(), args.data());
  mxi::eprint("mxi %s: cannot run %s: %s\n", command, script.c_str(),
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
    list(mxi::Stream::out);
    return 0;
  }
  if (std::strcmp(argv[1], "--version") == 0) {
    mxi::print("mxi %s\n", MXI_VERSION);
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
  mxi::eprint("mxi: no command '%s'\n\n", argv[1]);
  list(mxi::Stream::err);
  return 2;
}
