#pragma once

// Everything a program writes for a person, to the terminal and to a log file
// in the working directory, as DIALS writes dials.find_spots.log and the
// rest: mxi::print and mxi::eprint (output.hh), which write both in the same
// call. Standard error goes to the log too, in order: a log without the
// warning or the error that ended a run would mislead whoever reads it later.
//
// Header-only, so the spot finder, which does not link the pipeline's
// library, can use it too.

#include <cerrno>
#include <cstdint>
#include <cstring>
#include <string>

#include "output.hh"

namespace mxi {

//: True for a run that only prints help or a version. Such a run writes no
//: log, since it would otherwise overwrite the log of the last real run.
inline bool only_asks_for_help(int argc, char **argv) {
  for (std::int32_t i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "-h" || a == "--help" || a == "--version")
      return true;
  }
  return false;
}

//: The log opened for the rest of the run: what is printed goes to the
//: terminal and to it.
inline void mirror_to_log(const char *path) {
  if (!open_log(path))
    eprint("cannot write the log %s: %s; carrying on without it\n", path,
           std::strerror(errno));
}

} // namespace mxi
