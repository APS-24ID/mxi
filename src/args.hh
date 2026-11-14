// Splitting a command line into positional arguments and options.
//
// Written because all three programs originally read their file names from
// argv[1] and argv[2] and looked for options from argv[3] onwards, so
//
//     mxi_refine --beam indexed.expt indexed.refl
//
// took "--beam" as the experiment list, "indexed.expt" as the reflections, and
// then complained that "indexed.refl" was an unknown option. Options before
// positionals is what every other command line tool accepts, and the failure
// was confusing rather than merely inconvenient: the error named the wrong
// argument entirely.
//
// It lives here rather than in apps/ so that it can be tested. Argument
// parsing is exactly the kind of code that never gets tested because it feels
// too simple to get wrong.

#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>

namespace mxi {

struct Arguments {
  std::vector<std::string> positional;
  //: Flag to value. A flag that takes no value maps to the empty string.
  //: Given more than once, the last.
  std::map<std::string, std::string> options;
  //: Every value a flag was given, in order: for the options that may be
  //: repeated, as mxi_scale's --exclude.
  std::map<std::string, std::vector<std::string>> repeated;
  bool ok = true;
  std::string error;
  bool help = false;

  bool has(const std::string &flag) const { return options.count(flag) > 0; }
  std::string value(const std::string &flag,
                    const std::string &fallback = "") const {
    auto it = options.find(flag);
    return it == options.end() ? fallback : it->second;
  }
  double number(const std::string &flag, double fallback) const;
  std::vector<std::string> values(const std::string &flag) const {
    auto it = repeated.find(flag);
    return it == repeated.end() ? std::vector<std::string>{} : it->second;
  }
};

// `takes_value` names the flags that consume the following argument. Anything
// beginning with a dash and not in `known` is an error, which is deliberate:
// silently ignoring a misspelled option is how a run comes to use settings
// nobody chose.
//
// `optional_value` names flags whose value may be left out: such a flag takes
// the following argument only if it is a whole number, so that
// `--scan-varying 36` is 36 and `--scan-varying indexed.expt` is the flag alone
// and then a file.
Arguments parse_arguments(int argc, char **argv,
                          const std::set<std::string> &known,
                          const std::set<std::string> &takes_value,
                          const std::set<std::string> &optional_value = {});

//: One option of this run passed on to a run of the command for a part of
//: it -- a sweep, a lattice, a group refined apart -- appended to `words`:
//: with its value where it takes one (`takes_value`, the command's own). The
//: thread count once, as --threads N: parse_arguments holds it under both -j
//: and --threads, which no command's takes_value need list, and each part's
//: run had been given "-j" bare, then "--threads" -- not a number -- and
//: stopped (13 November: mxi integrate on two sweeps, test_integrate_sweeps;
//: so too mxi index --no-joint and mxi refine a group at a time, given -j).
void pass_on(std::vector<std::string> &words, const std::string &flag,
             const std::string &value,
             const std::set<std::string> &takes_value);

} // namespace mxi
