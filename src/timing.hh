#pragma once

// Where a run's time goes, for --timing: named phases of wall-clock time, and
// time summed across threads for work done in parallel. One table format for
// every program. Header-only, so that the spot finder, which does not link the
// pipeline's library, prints the same table.

#include <atomic>
#include <chrono>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include <sys/resource.h>

namespace mxi {

//: The process's peak resident memory so far, in MB: kB on Linux, bytes on
//: macOS. Only ever rises, so the phase after which it jumps is the one that
//: asked for the memory (Graeme: symmetry 22 GB, scaling 26, on 36000 images).
inline double peak_memory_mb() {
  struct rusage u {};
  if (getrusage(RUSAGE_SELF, &u) != 0)
    return 0.0;
#ifdef __APPLE__
  return static_cast<double>(u.ru_maxrss) / (1024.0 * 1024.0);
#else
  return static_cast<double>(u.ru_maxrss) / 1024.0;
#endif
}

class Timing {
public:
  //: `on` is whether --timing was asked for; when not, nothing is printed,
  //: and timing costs a clock read per phase.
  explicit Timing(bool on = false) : on_(on), began_(now()) {}
  //: A run that already keeps its own clock: `began` from Timing::now().
  Timing(bool on, double began) : on_(on), began_(began) {}

  static double now() {
    return std::chrono::duration<double>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
  }
  bool on() const { return on_; }

  //: A phase measured elsewhere; depth indents it under the one before.
  void add(const std::string &name, double seconds, int depth = 0) {
    phases_.push_back({name, seconds, depth, peak_memory_mb()});
  }

  //: A scope timed to its end: `{ auto t = timing.scope("reading"); ... }`.
  class Scope {
  public:
    Scope(Timing *t, std::string name, int depth)
        : t_(t), name_(std::move(name)), depth_(depth), start_(now()) {}
    Scope(const Scope &) = delete;
    Scope &operator=(const Scope &) = delete;
    ~Scope() { t_->add(name_, now() - start_, depth_); }

  private:
    Timing *t_;
    std::string name_;
    int depth_;
    double start_;
  };
  Scope scope(const std::string &name, int depth = 0) {
    return Scope(this, name, depth);
  }

  //: A line printed after the table: thread counts, what was read.
  void note(const std::string &text) { notes_.push_back(text); }

  //: The table, each phase in seconds and per cent of the run so far, and the
  //: run's total.
  void report(std::FILE *out = stdout) const {
    if (!on_)
      return;
    const double total = now() - began_;
    std::fprintf(out, "\nTiming\n");
    for (const Phase &p : phases_) {
      const std::string name =
          std::string(static_cast<std::size_t>(2 * p.depth), ' ') + p.name;
      std::fprintf(out, "  %-42s %9.3f s  %5.1f%%  peak %8.0f MB\n",
                   name.c_str(), p.seconds,
                   total > 0.0 ? 100.0 * p.seconds / total : 0.0, p.peak_mb);
    }
    std::fprintf(out, "  %-42s %9.3f s\n", "total", total);
    for (const std::string &n : notes_)
      std::fprintf(out, "  %s\n", n.c_str());
  }

private:
  struct Phase {
    std::string name;
    double seconds;
    int depth;
    double peak_mb; //: the process's peak memory as the phase ended
  };
  bool on_;
  double began_;
  std::vector<Phase> phases_;
  std::vector<std::string> notes_;
};

//: Seconds summed across threads, added to from any of them: the work a
//: parallel phase did, set against threads x wall time, says which part of it
//: the rest waited on.
class ThreadSeconds {
public:
  void add(double seconds) {
    nanoseconds_.fetch_add(static_cast<long long>(seconds * 1e9),
                           std::memory_order_relaxed);
  }
  double seconds() const {
    return static_cast<double>(nanoseconds_.load(std::memory_order_relaxed)) *
           1e-9;
  }

private:
  std::atomic<long long> nanoseconds_{0};
};

} // namespace mxi
