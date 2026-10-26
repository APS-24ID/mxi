// A progress bar, as tqdm draws one, on a terminal:
//
//   integrating  52%|█████████████▌            | 960/1700 images
//   [00:14<00:13, 68.6 images/s]
//
// redrawn in place at most ten times a second, written straight to the
// terminal's file descriptor -- not through stdout or stderr, which the log
// mirror (log_mirror.hh) copies to the program's log, where a line redrawn
// with carriage returns would be noise. Where standard error is not a
// terminal -- a log being read later, a batch job -- a plain line each tenth of
// the way instead, through stdout, which the log keeps (Graeme, 23 October).

#pragma once

#include <unistd.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace mxi {

class Progress {
public:
  //: `total` units of `unit` ("images"), labelled `label` ("integrating").
  Progress(std::string label, double total, std::string unit)
      : label_(std::move(label)), unit_(std::move(unit)),
        total_(total > 0.0 ? total : 1.0), start_(clock::now()),
        terminal_(on_terminal()) {}

  //: Whether a bar is drawn: standard error a terminal, and not a dumb one.
  static bool on_terminal() {
    if (!::isatty(STDERR_FILENO))
      return false;
    const char *term = std::getenv("TERM");
    return term == nullptr || std::strcmp(term, "dumb") != 0;
  }

  //: `done` of the total so far.
  void update(double done) {
    done_ = std::fmin(std::fmax(done, 0.0), total_);
    const double elapsed = seconds();
    if (terminal_) {
      if (elapsed - drawn_at_ < 0.1 && done_ < total_)
        return;
      drawn_at_ = elapsed;
      draw(elapsed, false);
      return;
    }
    const int tenth = static_cast<int>(10.0 * done_ / total_);
    if (tenth > tenth_ && tenth < 10) {
      tenth_ = tenth;
      const double left =
          done_ > 0.0 ? elapsed * (total_ - done_) / done_ : 0.0;
      std::printf("  %3d%%: %.0f of %.0f %s, %.0f s; about %.0f s left\n",
                  tenth * 10, done_, total_, unit_.c_str(), elapsed, left);
      std::fflush(stdout);
    }
  }

  //: The bar at the end, and the line closed: on a terminal only -- the plain
  //: lines need no closing.
  void finish() {
    if (finished_)
      return;
    finished_ = true;
    if (terminal_) {
      done_ = total_;
      draw(seconds(), true);
    }
  }

  ~Progress() { finish(); }

private:
  using clock = std::chrono::steady_clock;

  double seconds() const {
    return std::chrono::duration<double>(clock::now() - start_).count();
  }

  static std::string clock_text(double s) {
    const long t = static_cast<long>(s + 0.5);
    char buffer[32];
    if (t >= 3600)
      std::snprintf(buffer, sizeof buffer, "%ld:%02ld:%02ld", t / 3600,
                    (t / 60) % 60, t % 60);
    else
      std::snprintf(buffer, sizeof buffer, "%02ld:%02ld", t / 60, t % 60);
    return buffer;
  }

  void draw(double elapsed, bool last) const {
    static const char *const kEighths[] = {"",  "▏", "▎", "▍",
                                           "▌", "▋", "▊", "▉"};
    constexpr int kWidth = 30;
    const double fraction = done_ / total_;
    const int eighths = static_cast<int>(fraction * kWidth * 8.0);
    std::string bar;
    for (int i = 0; i < eighths / 8; ++i)
      bar += "█";
    int drawn = eighths / 8;
    if (drawn < kWidth && eighths % 8 != 0) {
      bar += kEighths[eighths % 8];
      ++drawn;
    }
    bar.append(static_cast<std::size_t>(kWidth - drawn), ' ');
    const double rate = elapsed > 0.0 ? done_ / elapsed : 0.0;
    const double left = rate > 0.0 ? (total_ - done_) / rate : 0.0;
    char line[512];
    std::snprintf(line, sizeof line,
                  "\r  %s %3.0f%%|%s| %.0f/%.0f %s [%s<%s, %.1f %s/s]\033[K%s",
                  label_.c_str(), 100.0 * fraction, bar.c_str(), done_, total_,
                  unit_.c_str(), clock_text(elapsed).c_str(),
                  clock_text(left).c_str(), rate, unit_.c_str(),
                  last ? "\n" : "");
    const std::size_t n = std::strlen(line);
    std::size_t written = 0;
    while (written < n) {
      const ssize_t k = ::write(STDERR_FILENO, line + written, n - written);
      if (k <= 0)
        break;
      written += static_cast<std::size_t>(k);
    }
  }

  std::string label_, unit_;
  double total_, done_ = 0.0;
  clock::time_point start_;
  bool terminal_;
  double drawn_at_ = -1.0;
  int tenth_ = 0;
  bool finished_ = false;
};

} // namespace mxi
