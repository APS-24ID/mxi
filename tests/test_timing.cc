#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>

#include "../src/timing.hh"
#include "check.hh"

namespace mxi {

TEST(thread_seconds_sum_what_every_thread_adds) {
  ThreadSeconds total;
  std::vector<std::thread> threads;
  for (std::int32_t t = 0; t < 8; ++t)
    threads.emplace_back([&] {
      for (std::int32_t k = 0; k < 1000; ++k)
        total.add(0.001);
    });
  for (std::thread &t : threads)
    t.join();
  check::close(total.seconds(), 8.0, 1e-6, "8 threads x 1000 x 1 ms");
}

TEST(timing_prints_nothing_unless_asked) {
  // Through a log, as a program's output reaches one (output.hh).
  const auto read = [](const std::string &path) {
    std::FILE *f = std::fopen(path.c_str(), "rb");
    char buffer[512] = {0};
    const std::size_t n = std::fread(buffer, 1, sizeof buffer - 1, f);
    std::fclose(f);
    return std::string(buffer, n);
  };
  const std::string path = check::scratch_path("timing.log");
  check::is_true(open_log(path.c_str()), "a log");
  Timing off(false);
  off.add("a phase", 1.0);
  off.report();
  close_log();
  check::is_true(read(path).empty(), "nothing without --timing");
  check::is_true(open_log(path.c_str()), "a log again");
  Timing on(true);
  on.add("a phase", 1.0);
  on.add("its part", 0.5, 1);
  on.report();
  close_log();
  const std::string text = read(path);
  // Every line inside 80 columns: the names' column as wide as its widest.
  std::size_t at = 0;
  while (at < text.size()) {
    const std::size_t end = std::min(text.find('\n', at), text.size());
    check::is_true(end - at <= 80, "a line of " + std::to_string(end - at) +
                                       ": " + text.substr(at, end - at));
    at = end + 1;
  }
  check::is_true(text.find("\nTiming\n") != std::string::npos,
                 "a table with --timing");
  check::is_true(text.find("    its part") != std::string::npos,
                 "the part indented under it");
  check::is_true(text.find("  total") != std::string::npos,
                 "and the run's total");
}

} // namespace mxi
