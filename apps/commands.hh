#pragma once

// mxi's commands by name, for what runs one in code: mxi process, each step
// a call where it had been a process of its own (docs/windows.md).

#include <cstdint>

namespace mxi {

//: A command's entry, as main's: argv[0] the command as typed ("mxi find").
using Command = int (*)(int, char **);

//: The command of that name ("find"), or nullptr: those of the table in
//: apps/mxi.cc, the Python ones (view, ui) not among them.
Command command_named(const char *name);

} // namespace mxi
