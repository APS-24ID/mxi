// Where the running program is: mxi's dispatcher finds the Python beside it
// (mxi view, mxi ui) and mxi process runs mxi's commands with it.

#ifndef MXI_SELF_PATH_HH
#define MXI_SELF_PATH_HH

#include <string>

namespace mxi {

//: This program's path, symbolic links resolved: /proc/self/exe on Linux,
//: _NSGetExecutablePath on macOS, else argv0 -- as a path where it has a
//: slash, looked for on the PATH where it has none.
std::string executable_path(const char *argv0);

} // namespace mxi

#endif // MXI_SELF_PATH_HH
