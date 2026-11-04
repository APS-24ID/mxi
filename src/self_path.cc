#include "self_path.hh"

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <unistd.h>

#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif

namespace mxi {

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

} // namespace mxi
