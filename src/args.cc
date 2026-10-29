#include "args.hh"

#include "parallel.hh"

#include <cstdio>

#include <cstdlib>

namespace mxi {

double Arguments::number(const std::string &flag, double fallback) const {
  auto it = options.find(flag);
  if (it == options.end() || it->second.empty())
    return fallback;
  char *end = nullptr;
  const double v = std::strtod(it->second.c_str(), &end);
  return end == it->second.c_str() ? fallback : v;
}

Arguments parse_arguments(int argc, char **argv,
                          const std::set<std::string> &known,
                          const std::set<std::string> &takes_value,
                          const std::set<std::string> &optional_value) {
  // --version, for every program: its name and mxi's version, and nothing
  // else done (once only mxi_find answered it).
  for (int i = 1; i < argc; ++i)
    if (std::string(argv[i]) == "--version") {
      std::string name = argv[0];
      const std::size_t slash = name.find_last_of('/');
      if (slash != std::string::npos)
        name = name.substr(slash + 1);
      std::printf("%s %s\n", name.c_str(), MXI_VERSION);
      std::exit(0);
    }
  Arguments out;
  const auto whole_number = [](const std::string &s) {
    if (s.empty())
      return false;
    for (char ch : s)
      if (ch < '0' || ch > '9')
        return false;
    return true;
  };
  for (int i = 1; i < argc; ++i) {
    // -j N and --threads N, one option, taken by every command alike: the
    // threads for everything the command does in parallel, 0 (the default)
    // being every core this job may use (parallel.hh). Commands had asked it
    // differently, or of one part only -- --fft-threads, --jacobian-threads,
    // --normal-threads -- so that a job could not simply be held to a share of
    // a node (APS-24ID/mxi#1).
    std::string arg = argv[i];
    if (arg == "-j" || arg == "--threads") {
      if (i + 1 >= argc || !whole_number(argv[i + 1])) {
        out.ok = false;
        out.error = "option '" + arg +
                    "' needs a number of threads, 0 for "
                    "every core this job may use";
        return out;
      }
      const std::string value = argv[++i];
      out.options["-j"] = value;
      out.options["--threads"] = value;
      out.repeated[arg].push_back(value);
      set_parallel_threads(static_cast<std::size_t>(std::stoul(value)));
      continue;
    }
    if (arg == "-h" || arg == "--help") {
      // Both ways of asking: args.help, and has("--help") as seven programs
      // asked it -- which this never set, so their --help was an error.
      out.help = true;
      out.options["--help"] = "";
      out.options["-h"] = "";
      return out;
    }
    // A bare "-" is a filename by convention, not an option.
    if (arg.size() > 1 && arg[0] == '-') {
      if (!known.count(arg)) {
        out.ok = false;
        out.error = "unknown option '" + arg + "'";
        return out;
      }
      if (takes_value.count(arg)) {
        if (i + 1 >= argc) {
          out.ok = false;
          out.error = "option '" + arg + "' needs a value";
          return out;
        }
        out.options[arg] = argv[i + 1];
        out.repeated[arg].push_back(argv[++i]);
      } else if (optional_value.count(arg)) {
        if (i + 1 < argc && whole_number(argv[i + 1])) {
          out.options[arg] = argv[i + 1];
          out.repeated[arg].push_back(argv[++i]);
        } else {
          out.options[arg] = "";
        }
      } else {
        out.options[arg] = "";
      }
      continue;
    }
    out.positional.push_back(arg);
  }
  return out;
}

} // namespace mxi
