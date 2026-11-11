#pragma once

// What mxi needs of the system that C++ does not give, in one place: each a
// POSIX half (Linux, macOS) and a Windows half (docs/windows.md, the platform
// layer). Header-only, as progress.hh and timing.hh are, for the spot
// finder, which does not link the pipeline's library.
//
//   stderr_is_terminal, stderr_columns -- progress bars (progress.hh)
//   peak_memory_mb                     -- --timing (timing.hh)
//   open_for_reading, read_at          -- chunks read outside HDF5's lock,
//                                         from many threads (nxmx_reader.cc)

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <io.h>
#include <windows.h>
// K32GetProcessMemoryInfo is kernel32's own: no psapi library to link.
#include <psapi.h>
#else
#include <fcntl.h>
#include <sys/ioctl.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

namespace mxi::platform {

//: Whether standard error is a terminal (not a file or a pipe).
inline bool stderr_is_terminal() {
#if defined(_WIN32)
  return _isatty(_fileno(stderr)) != 0;
#else
  return ::isatty(STDERR_FILENO) != 0;
#endif
}

//: The width of standard error's terminal in columns; 0 where it cannot say.
inline std::int32_t stderr_columns() {
#if defined(_WIN32)
  CONSOLE_SCREEN_BUFFER_INFO info;
  if (GetConsoleScreenBufferInfo(GetStdHandle(STD_ERROR_HANDLE), &info) != 0)
    return static_cast<std::int32_t>(info.srWindow.Right - info.srWindow.Left +
                                     1);
  return 0;
#else
  struct winsize w {};
  if (::ioctl(STDERR_FILENO, TIOCGWINSZ, &w) == 0 && w.ws_col > 0)
    return static_cast<std::int32_t>(w.ws_col);
  return 0;
#endif
}

//: The process's peak resident memory so far, in MB; 0 where it cannot say.
//: Only ever rises, so the phase after which it jumps is the one that asked
//: for the memory (Graeme: symmetry 22 GB, scaling 26, on 36000 images).
inline double peak_memory_mb() {
#if defined(_WIN32)
  PROCESS_MEMORY_COUNTERS counters{};
  if (K32GetProcessMemoryInfo(GetCurrentProcess(), &counters,
                              sizeof counters) == 0)
    return 0.0;
  return static_cast<double>(counters.PeakWorkingSetSize) / (1024.0 * 1024.0);
#else
  struct rusage u {};
  if (getrusage(RUSAGE_SELF, &u) != 0)
    return 0.0;
#if defined(__APPLE__)
  return static_cast<double>(u.ru_maxrss) / (1024.0 * 1024.0); // bytes
#else
  return static_cast<double>(u.ru_maxrss) / 1024.0; // kB
#endif
#endif
}

//: A file opened for reading at offsets, by any number of threads at once:
//: a descriptor on POSIX, a handle on Windows; kNoFile where it would not
//: open.
using File = std::intptr_t;
inline constexpr File kNoFile = -1;

//: `path`, UTF-8, opened to be read with read_at.
inline File open_for_reading(const std::string &path) {
#if defined(_WIN32)
  const std::int32_t n =
      MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, nullptr, 0);
  if (n <= 0)
    return kNoFile;
  std::wstring wide(static_cast<std::size_t>(n), L'\0');
  MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, wide.data(), n);
  const HANDLE h = CreateFileW(
      wide.c_str(), GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr,
      OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL | FILE_FLAG_RANDOM_ACCESS, nullptr);
  return h == INVALID_HANDLE_VALUE ? kNoFile : reinterpret_cast<File>(h);
#else
  const std::int32_t fd = ::open(path.c_str(), O_RDONLY);
  return fd < 0 ? kNoFile : static_cast<File>(fd);
#endif
}

//: Up to `size` bytes of `file` at `offset` read into `into`, as POSIX's
//: pread: how many were read, 0 at the end of the file, -1 on an error. Safe
//: from many threads at once on one file, each naming its own offset.
inline std::int64_t read_at(File file, void *into, std::size_t size,
                            std::uint64_t offset) {
#if defined(_WIN32)
  // At most what a DWORD counts in one call; the caller loops for the rest.
  const DWORD ask = static_cast<DWORD>(size > 0x40000000u ? 0x40000000u : size);
  OVERLAPPED at{};
  at.Offset = static_cast<DWORD>(offset & 0xffffffffu);
  at.OffsetHigh = static_cast<DWORD>(offset >> 32);
  DWORD got = 0;
  if (ReadFile(reinterpret_cast<HANDLE>(file), into, ask, &got, &at) == 0)
    return GetLastError() == ERROR_HANDLE_EOF ? 0 : -1;
  return static_cast<std::int64_t>(got);
#else
  return static_cast<std::int64_t>(::pread(
      static_cast<std::int32_t>(file), into, size, static_cast<off_t>(offset)));
#endif
}

} // namespace mxi::platform
