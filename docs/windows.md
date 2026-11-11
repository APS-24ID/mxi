# A Windows build: the plan, and the integer audit

STATUS: investigation and documentation only, 7 November 2026 (Graeme: "Write
up the plan and the long audit ... this should only be an investigation and
documentation exercise"). Nothing here is built yet. `docs/outstanding.md`, 89.

Graeme has no Windows development machine, but can test; the builds would run
on GitHub Actions' Windows runners, as the Linux and macOS release builds run
on theirs.

## The answer

Possible, and a moderate port, not a large one. What mxi depends on builds on
Windows: HDF5, zlib and gemmi through vcpkg, which the release builds already
use; bitshuffle; on the Python side tkinter (with python.org's Python), numpy,
h5py, hdf5plugin and matplotlib, all with Windows wheels. What is POSIX in mxi
is in six places, each with a direct Windows equivalent; and one class of
silent error -- `long` is 32 bits on 64-bit Windows -- wants auditing first,
done below.

## What is POSIX, and its Windows equivalent

Found by searching the source for system headers and calls (src, apps):

| Where | What | On Windows |
| --- | --- | --- |
| `src/self_path.cc` | the executable's own path: `/proc/self/exe`, `_NSGetExecutablePath`, `realpath` | **gone, 7 November**: nothing needs it |
| `apps/mxi_process.cc` | each step started and waited for: `posix_spawn`, `waitpid` | **gone, 7 November**: each step a call |
| `src/progress.hh` | whether a terminal, and its width: `isatty`, `ioctl(TIOCGWINSZ)`; its bar written with `::write` | **done, 7 November**: `platform.hh`'s `stderr_is_terminal`, `stderr_columns`; the bar through C's own `stderr` |
| `src/timing.hh` | peak memory for `--timing`: `getrusage` | **done, 7 November**: `platform.hh`'s `peak_memory_mb` (`K32GetProcessMemoryInfo` on Windows, no psapi to link) |
| `src/log_mirror.hh` | stdout and stderr written to the terminal and a log at once: the C library's `stdout` and `stderr` replaced by streams made with `fopencookie` (glibc) or `funopen` (BSD) | **done, 7 November**: `src/output.hh`, below -- standard C++ only, nothing to port |
| `src/nxmx_reader.cc` | direct reads of NXmx files: `open`, `pread` | **done, 7 November**: `platform.hh`'s `open_for_reading`, `read_at` (`CreateFileW` on the path as UTF-16, `ReadFile` at an offset) |
| `apps/mxi.cc` | `mxi view`, `mxi ui`: `execv` of a Python script with a `#!` line | **gone, 7 November**: `mxi-view` and `mxi-ui` are commands of their own |
| `python/src/mxi_ui/locate.py`, `mxi_tools.py` | `mxi` looked for beside the package by that name | **done, 7 November**: `mxi.exe` where `os.name` is `nt` |
| 10 files | `M_PI` | **done, 7 November**: C++20's `std::numbers::pi`, the same value, no macro |
| `src/processing.hh` | one `__attribute__` | **done, 7 November**: `MXI_PRINTF_LIKE` (output.hh), nothing under MSVC, `gnu_printf` under MinGW |
| `CMakeLists.txt` | `-Wall -Wextra -Wpedantic` on every target | by compiler: `/W4` for MSVC, or clang-cl, which takes them |

Not found: `fork`, `mmap`, `dlopen`, `pthread` (threads are `std::thread`),
`long double`, `H5T_NATIVE_LONG`. `mxi_readtest` (a diagnostic timing reads)
can be left out of a Windows build.

**Text mode.** On Windows a file opened without binary mode has its line
endings translated, which corrupts binary data. Every binary read is
`fopen(..., "rb")`; the four `fopen(..., "w")` and two `std::ofstream` without
`binary` write text (comment-headed tables, HTML), where translation is right;
the binary files mxi writes go through HDF5 and gemmi, which open them so
themselves. Nothing to change.

**The log mirror** (corrected 7 November: first written up here as file
descriptors and a pipe, which it is not). It replaces `stdout` and `stderr`
themselves, so that what a program writes goes to the terminal and its log in
the same call -- deliberately not a pipe copied by a thread, since a crash's
last lines can die in the pipe, and those are the lines that matter. Windows
offers neither way of making such a stream. Two designs:

* **All output through mxi's own two streams** -- `std::fprintf(mxi::out,
  ...)` where it is `printf(...)` now: some 900 calls (715 printf, 193
  fprintf to stdout or stderr; no iostreams), changed mechanically; the same
  on every system, and testable. The larger change and the better one.
* **On Windows only, a pipe copied by a thread**: little code, but the design
  rejected for the crash case, and a second way of doing the same thing.

**Done, 7 November: the first design.** `src/output.hh`: `mxi::print`,
`mxi::eprint`, `mxi::print_to` and `mxi::put` in place of `printf`,
`fprintf(stdout or stderr, ...)` and `fputs` -- 651 calls in src and apps --
each formatting once and writing the terminal and the open log under one
lock, the log flushed each call; stdout flushed at each newline, piped or
not, as the mirror's was (mxi ui shows a step's output as it is written);
logs a stack, `open_log` and `close_log`, for mxi process in code. The four
functions that took a `FILE *` always stdout or stderr (`print_summary`, the
timing table's `report`, find's `usage`, mxi's `list`) take an
`mxi::Stream`. `log_mirror.hh` now only opens the log; no `fopencookie`,
`funopen`, `unistd` or `fcntl`. Not mirrored, as the old mirror did: what
another library writes straight to stderr (HDF5's error stack where not
turned off) -- on the terminal only. The Metal and CUDA sources keep
`printf` (not compiled here) -- their output on the terminal, not in a log,
until converted where they build. A run's log and its terminal output byte
for byte the same (`mxi index` on the planted series); tests of order, nested
logs, long text and threads' whole lines (`tests/test_output.cc`).

## Simpler first: no process starting another

Graeme asked (7 November) how much it would simplify -- the port, and the code
-- if mxi started no other process: `mxi process` running each step in code,
and the interface started only as `mxi-ui` and `mxi-view`.

**`mxi process` in code** (done, 7 November: `apps/commands.hh`'s
`command_named`, each step's log closed as it returns -- `close_logs_to`;
its posix_spawn, waitpid and use of self_path gone; a run that stops at a
failing step, its logs as before). mxi is one binary whose commands are functions
(`mxi_find_main` and the rest), so `mxi process` can call them. Gone:
`posix_spawn`, `waitpid`, `environ`, and its use of `src/self_path.cc`.
Needed, each modest: the log mirror made a scope -- pushed as a step starts,
popped as it ends, its file closed and the output before it restored, where
each step now relies on its process ending; three `std::exit` (`--help` and
`--version`, which mxi process never passes) made returns; the statics of
`apps/mxi_find.cc` (its sweeps' tables and totals, 105 and 1290) into the
call; `--timing`'s peak memory, which a process cannot reset between steps,
said as the run's peak so far or measured otherwise; and every HDF5 handle
closed when a step returns -- a process's exit closed any left open, and a
step's output left open could fail the next step's reading of it. Lost: a
crash's isolation -- a step that crashes takes mxi process with it, without
its message of which step stopped and how to go on; a step that fails and
returns keeps it.

**`mxi-ui` and `mxi-view` only** (done, 7 November: `mxi ui` and `mxi
view`, typed, say where they went; the build installs `mxi-view` beside
`mxi-ui`; `src/self_path.cc` deleted, nothing using it). Gone: mxi's starting of Python
(`python_command` in `apps/mxi.cc`, `execv` of a `#!` script) and the other
use of `src/self_path.cc`, which goes with both changes made. On Windows pip
makes the wheel's entry points small `.exe` launchers itself: nothing to port.
The cost is to users: `mxi view` and `mxi ui` are in the README, the
tutorials, the guides and the tests, and the release archive's interface,
run today from its `bin`, would want installing from the wheel -- a change
for a release's notes.

**What Windows would still need**, both made: the terminal's width
(trivial), peak memory (small), the NXmx reads (small) -- the log mirror
now done -- three places of six, and none of them starting or waiting for a
process, the part hardest to make the same on three systems.

**The platform layer, done 7 November**: `src/platform.hh`, header-only, a
POSIX half and a Windows half -- `stderr_is_terminal`, `stderr_columns`,
`peak_memory_mb`, `open_for_reading` and `read_at`. Outside it, no POSIX
header in the code but the diagnostic `mxi_readtest`, left out of a Windows
build.

**Checked without Windows** (`scripts/check_windows.sh`): every source
compiled for Windows by MinGW-w64, syntax only, warnings on -- all 76, no
warning in mxi's code (with MinGW's C99 stdio, `__USE_MINGW_ANSI_STDIO`: its
default, Microsoft's old runtime, knows no `%zu`; the current one, the planned
build's, does) -- and the platform layer's Windows half built and run under
Wine (`tests/windows/platform_windows.cc`): 2000 reads at offsets from four
threads, every byte right; the end of a file; no file; peak memory; a log
holding what was printed. Wine is not Windows, but the same API is called. And `scripts/check_libcxx.sh`: every source and test against libc++, Apple's
C++ library, for what GCC's own library lets through (126 files, none failing).
Left for the build job: CMake for clang-cl (its flags), vcpkg's Windows
triplet, the workflow's job; HDF5 and paths that are not ASCII.

## The integer audit

64-bit Windows is LLP64: `long` and `unsigned long` are 32 bits, where Linux and
macOS (LP64) make them 64; `long long`, `std::size_t`, `std::ptrdiff_t` and the
`<cstdint>` types are the same on all three. A `long` holding a value past 2^31
compiles on Windows and is silently wrong there. 253 uses of `long` in `src` and
`apps` (46 more in `tests`), each read, 7 November:

### Wrong on Windows: to change before a Windows build

* **`src/rflx.cc` 1015--1026: the shoeboxes' offsets.** A file's shoeboxes are
  held in one flat array of every reflection's voxels, and each reflection's
  are found by `begin() + static_cast<long>(at)`. Two million reflections of a
  thousand voxels are 2e9 voxels, past 2^31: on Windows the offset would wrap,
  and a reflection be given another's pixels, without a word.
* **`src/nxmx_reader.cc` 446: a signed overflow.** The search for HDF5's
  signature, `for (long at = 0; ... at <= (1L << 30); at = at == 0 ? 512 :
  at * 2)`: at 2^30 the body runs and `at * 2` is 2^31 -- on Windows signed
  overflow, undefined. In practice negative, and `fseek` to it fails and ends
  the loop, so no hang (as first written here, and corrected on reading the
  loop whole); reached only for a file of 1 GB or more with no signature at
  any offset tried; but undefined behaviour a compiler may optimise on. Its
  `fseek` takes a `long` too: on Windows it cannot seek past 2 GB --
  `_fseeki64` there, should a seek ever need to.
* **`src/index.cc` 288--371: the max cell's grid of cells.** The cell key
  `(x ny + y) nz + z` is a `long`; in the hashed branch, kept for grids far
  larger than the points, the product can pass 2^31. Measured on 428 349 spots
  of a 128 A cell the grid is some 6e5 cells -- safe there, but nothing in the
  type says so.
* **`apps/mxi_find.cc` 1290: running totals of signal pixels and spots** over
  a whole run, which a large enough run takes past 2^31.

### Safe at real sizes, but the wrong type

* **Iterator offsets**, `begin() + static_cast<long>(...)`: `src/index.cc`
  (1493, 1633, 1637, 1675), `src/laue.cc` 375, `src/correlation.cc` 507 and
  736, `apps/mxi_import.cc` 156, `apps/mxi_process.cc` 279,
  `apps/mxi_grid.cc` 213. An iterator's difference type is `std::ptrdiff_t`,
  64 bits everywhere; `3 i` passes 2^31 only past 700 million rows.
* **Values into the summary**: 81 `json::Value(static_cast<long>(...))` and
  `json::Value(long)`. `json::Value` holds every number as a double, exact to
  2^53; only the cast going in is 32 bits. What is cast is counts --
  reflections, observations, outliers, images, lattices, cluster sizes,
  milliseconds -- below 2^31 in any data seen; a summary of a run with more
  than 2^31 of anything would be wrong.

### Safe: small by their nature

Experiments, sweeps, lattices and clusters numbered (`src/correlation.cc`,
`src/summary_cut.cc`, `src/views.cc`, `src/combine.cc`, `src/expt.cc`,
`apps/mxi_refine.cc`, `apps/mxi_integrate.cc`, `apps/mxi_report.cc`); the
indexing transform's grid, a few hundred a side, and its loops (`src/index.cc`
531--616); Miller indices (`src/index.cc` 1306); symmetry operators' integer
matrices (`apps/mxi_symmetry.cc`); B-spline knots (`src/geometry.cc`,
`src/derivatives.cc`); pixel positions (`src/verify.cc`); a duration in seconds
(`src/progress.hh`); detector modules (`apps/mxi_residuals.cc`); a decimal's
digits (`src/processing.cc`); an image range parsed (`apps/mxi_scale.cc`).

### Formats and parsing

Six `%ld` (`src/processing.cc`, `src/progress.hh`, `apps/mxi_symmetry.cc`,
`apps/mxi_scale.cc` twice, `apps/mxi_residuals.cc`), each with `long`
arguments: consistent now, to change with them (`%lld` with a `long long`, or
`<cinttypes>`' `PRId64`). `std::strtol` in `apps/mxi_split.cc` and `std::stoul`
in `src/args.cc` and `apps/mxi_untrusted.cc` parse small numbers. The Metal
code's `unsigned long` (`src/metal_device.*`, `src/fit_device_metal.cc`,
`src/dext_metal.cc`) is Apple's API, macOS only: left.

## Sized types: a rule for mxi, Windows or not

Graeme asked whether explicitly sized integer types are a good idea anyway,
and then decided (7 November): explicit everywhere, for consistency -- a
reader then sees the range each value is meant to hold, and one rule with no
judgement in it is the easier kept. `int` is 32 bits on all three systems, so
for it the change is of what the code says, not what it does.

* **`std::int64_t` / `std::uint64_t`** for what grows with the data -- rows,
  voxels, pixels, offsets into flat arrays, keys, running totals -- and where
  `long` or `long long` was.
* **`std::int32_t` / `std::uint32_t`** where `int` or `unsigned` was; the
  narrower ones (`std::int16_t`, `std::uint8_t`) where a value or a format is
  that narrow.
* **`std::size_t`** for container sizes and indices, as most of mxi has: it
  says what it is. **`std::ptrdiff_t`** for iterator offsets.
* **The summary's integers**: one `json::Value` constructor for every integer
  type (a template), where a list of `int`, `long`, `long long` and `size_t`
  was a different list on each system.
* **Left as an outside interface fixes them**: `main(int argc, char **argv)`;
  `printf`'s `%*d` widths; HDF5's `hid_t`, `herr_t` and the like; gemmi's,
  Metal's (`unsigned long`) and the C library's own types, `std::lround`'s
  `long` among them, converted where their values enter mxi's.
* **Formats** by `<cinttypes>`' macros (`PRId64`) where a 64-bit value is
  printed or scanned.

**Stage 1 done, 7 November**: every `long`, `long long`, `unsigned long` and
`unsigned long long` in src, apps and tests now `std::int64_t` or
`std::uint64_t` (862, by a script that reads past comments and strings), the
literals' suffixes with them (`std::int64_t{0}`, unsigned ones keeping their
`u`), 92 formats by `<cinttypes>`' macros -- `SCNd64` in the two `sscanf`
calls -- and `json::Value`'s integer constructors one template. What is left
of `long` in the code is Metal's API and one test's `long double` (on Windows
the same as `double`, so its reference there would lose its extra
precision). On Linux and macOS no type changed -- `long` is `std::int64_t`
there -- and the suites pass as before: 378 C++ tests, 745 Python.

**Stage 2 done, 7 November**: `int` as `std::int32_t`, `unsigned` and
`unsigned int` as `std::uint32_t`, `short` as `std::int16_t` -- 1402 in 147
files, the same types on all three systems, so `%d` and `%u` stand and nothing
changes at run time. Kept: `main` and the commands' entry functions, `int
f(int argc, char **argv)` (48); `char` and `unsigned char`, characters and
bytes. **The accelerators converted too**, 7 November (Graeme: before the platform
layer): the CUDA sources (`src/*.cu`, `src/fit_device_cuda.cc`), the Metal
ones (`src/dext_metal.cc`, `src/fit_device_metal.cc`), `src/fft_fftw.cc` and
the GPU spot finder's two tests -- 254 replacements; their three host-side
`fprintf(stderr, ...)` through `mxi::eprint`. Kept: the Metal shader library's
embedded size and the function given it, `unsigned long` as the build
generates the symbol; and `printf` in the GPU tests, as every test prints.
Compiled here with CUDA 12.0's `nvcc` (GCC 11 its host compiler) and FFTW on,
with no warnings; the GPU test skips without a device; FFTW's transform agrees
with the built-in. The Metal sources compiled where they build: Graeme's Mac. 378 C++ tests, 745 Python; no warnings.

Done in stages, each committed alone and each a change of no behaviour on
Linux and macOS: the `long` family first (`long` is `std::int64_t` there,
the same type), then `int` and `unsigned`. On Windows the first stage is what
makes the four cases above right.

## The build

* **Toolchain**: Microsoft's, through vcpkg (`x64-windows-static`), one
  statically linked `mxi.exe` with no DLLs to find, as the Linux and macOS
  archives are static. clang-cl, Clang with Microsoft's ABI, keeps warnings and
  behaviour closest to the other builds and takes the `-W` flags; MSVC's own
  `cl` wants `/W4` and finds different things. clang-cl first.
* **CPU only.** CUDA on Windows can come later; on a hosted runner its toolkit
  is a large install, and the CPU paths are what to prove first. No Metal.
* **The workflow**: a Windows job in the release workflow beside Linux and
  macOS -- vcpkg's dependencies cached, mxi built, ctest and the Python suite
  run (Python 3.12 from actions/setup-python, with tkinter), the window tests
  skipped (`MXI_NO_WINDOW_TESTS`) as on the other runners, and a zip kept:
  `mxi.exe`, `mxi_ui/`, `mxi_view.py`, the licences. Python itself not in the
  zip: `mxi view` and `mxi ui` use the Python found on the system, as on the
  others.
* **Its own branch** until it passes, so that a Windows failure does not hold
  up a Linux or macOS release.

## What Graeme tests

Unzip; `mxi --version`; `mxi process` on a data set, its merged statistics
against a Linux or macOS run of the same (the same within the threads' order,
as between those two); `mxi ui` and `mxi view` opened, the windows CI cannot
open. A path with spaces in it, and one with a non-ASCII character -- HDF5 on
Windows takes paths in the system's code page unless built for UTF-8.

## Phases

1. **The integers**, by the rule above, on Linux and macOS: no Windows needed.
2. **No process starting another** (above), if Graeme chooses it, on Linux
   and macOS; then **the platform layer**: what POSIX remains behind
   `src/platform.{hh,cc}`, Windows' half compiled only there; the log mirror
   by its own design (above).
3. **The build**: CMake and vcpkg for Windows, the workflow job; round trips
   through Graeme until it builds and its tests pass.
4. **Graeme's test** of the zip, and a release with it.
5. Later, if wanted: CUDA on Windows.

## Open

* clang-cl or MSVC's `cl` (clang-cl planned).
* Whether the Windows zip goes in releases from its first passing build, or
  after a release or two on its branch.
* Paths with non-ASCII characters: HDF5's UTF-8 build, or a documented limit.
