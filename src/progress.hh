// A progress bar, as tqdm draws one, on a terminal:
//
//   integrating  52%|█████▋     | 960/1700 images [00:14<00:13, 68.6 images/s]
//
// the bar as wide as the terminal leaves room for, one column spare -- the
// rate and then the times given up on a narrow one -- redrawn in place at most
// ten times a second, written straight to the terminal's file descriptor -- not
// through stdout or stderr, which the log mirror (log_mirror.hh) copies to the
// program's log, where a line redrawn with carriage returns would be noise.
// Where standard error is not a terminal -- a log being read later, a batch job
// -- a plain line each tenth of the way instead, through stdout, which the log
// keeps (Graeme, 23 October).

#pragma once

#include <sys/ioctl.h>
#include <unistd.h>

#include <algorithm>
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

  // The terminal's width, asked each time it is drawn so that a window
  // resized while it runs is followed; COLUMNS, then 80, where it cannot say.
  static int columns() {
    struct winsize w{};
    if (::ioctl(STDERR_FILENO, TIOCGWINSZ, &w) == 0 && w.ws_col > 0)
      return w.ws_col;
    if (const char *c = std::getenv("COLUMNS"))
      if (const int n = std::atoi(c); n > 0)
        return n;
    return 80;
  }

  void draw(double elapsed, bool last) const {
    static const char *const kEighths[] = {"",  "▏", "▎", "▍",
                                           "▌", "▋", "▊", "▉"};
    const double fraction = done_ / total_;
    const double rate = elapsed > 0.0 ? done_ / elapsed : 0.0;
    const double left = rate > 0.0 ? (total_ - done_) / rate : 0.0;
    char head[128];
    std::snprintf(head, sizeof head, "  %s %3.0f%%", label_.c_str(),
                  100.0 * fraction);
    // What follows the bar, shortened until a bar of ten fits: the rate
    // first, then the times; on a terminal narrower still, no bar at all.
    char tails[3][256];
    std::snprintf(tails[0], sizeof tails[0],
                  "| %.0f/%.0f %s [%s<%s, %.1f %s/s]", done_, total_,
                  unit_.c_str(), clock_text(elapsed).c_str(),
                  clock_text(left).c_str(), rate, unit_.c_str());
    std::snprintf(tails[1], sizeof tails[1], "| %.0f/%.0f [%s<%s]", done_,
                  total_, clock_text(elapsed).c_str(),
                  clock_text(left).c_str());
    std::snprintf(tails[2], sizeof tails[2], "| %.0f/%.0f", done_, total_);
    // The bar takes what the terminal leaves, one column spare so the cursor
    // never wraps: a line wider than the terminal had wrapped on 80 columns
    // (Graeme, 23 October). Every character here is one column wide.
    const int cols = columns();
    const char *tail = nullptr;
    int width = 0;
    for (const char *t : tails) {
      const int room = cols - 1 - static_cast<int>(std::strlen(head)) - 1 -
                       static_cast<int>(std::strlen(t));
      if (room >= 10) {
        tail = t;
        width = std::min(60, room);
        break;
      }
    }
    if (tail == nullptr) {
      char line[160];
      std::snprintf(line, sizeof line, "\r%s\033[K%s", head, last ? "\n" : "");
      put(line);
      return;
    }
    const int eighths = static_cast<int>(fraction * width * 8.0);
    std::string bar;
    for (int i = 0; i < eighths / 8; ++i)
      bar += "█";
    int drawn = eighths / 8;
    if (drawn < width && eighths % 8 != 0) {
      bar += kEighths[eighths % 8];
      ++drawn;
    }
    bar.append(static_cast<std::size_t>(width - drawn), ' ');
    char line[1024];
    std::snprintf(line, sizeof line, "\r%s|%s%s\033[K%s", head, bar.c_str(),
                  tail, last ? "\n" : "");
    put(line);
  }

  // To the terminal itself, past the log mirror's copies of stdout and
  // stderr.
  static void put(const char *line) {
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
