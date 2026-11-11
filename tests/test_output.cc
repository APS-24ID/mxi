// output.hh: mxi::print and mxi::eprint, to the terminal and the open log in
// the same call -- the log mirror that Windows can have too (docs/windows.md).

#include <cstdint>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>

#include "../src/output.hh"
#include "check.hh"

namespace mxi {

namespace {

std::string read_all(const std::string &path) {
  std::FILE *f = std::fopen(path.c_str(), "rb");
  std::string text;
  char buffer[4096];
  std::size_t n;
  while ((n = std::fread(buffer, 1, sizeof buffer, f)) > 0)
    text.append(buffer, n);
  std::fclose(f);
  return text;
}

} // namespace

TEST(output_to_the_log_in_the_order_written) {
  const std::string path = check::scratch_path("output_order.log");
  check::is_true(open_log(path.c_str()), "a log");
  print("one %d\n", 1);
  eprint("two %s\n", "to stderr");
  put(Stream::out, "three\n");
  print_to(Stream::err, "four %.1f\n", 4.0);
  close_log();
  check::is_true(read_all(path) == "one 1\ntwo to stderr\nthree\nfour 4.0\n",
                 "out and err, in order: " + read_all(path));
}

TEST(output_logs_nest) {
  const std::string outer = check::scratch_path("output_outer.log"),
                    inner = check::scratch_path("output_inner.log");
  check::is_true(open_log(outer.c_str()), "the outer log");
  print("before\n");
  check::is_true(open_log(inner.c_str()), "the inner log");
  print("during\n");
  close_log();
  print("after\n");
  close_log();
  print("nowhere\n"); // the terminal only: no log open
  check::is_true(read_all(outer) == "before\nafter\n", "the outer, around it");
  check::is_true(read_all(inner) == "during\n", "the inner, only its own");
}

TEST(output_longer_than_its_buffer_whole) {
  const std::string path = check::scratch_path("output_long.log");
  const std::string long_text(5000, 'x');
  check::is_true(open_log(path.c_str()), "a log");
  print("%s|%d\n", long_text.c_str(), 5000);
  close_log();
  check::is_true(read_all(path) == long_text + "|5000\n",
                 "5000 characters and more");
}

TEST(output_from_threads_in_whole_lines) {
  const std::string path = check::scratch_path("output_threads.log");
  check::is_true(open_log(path.c_str()), "a log");
  std::vector<std::thread> threads;
  for (std::int32_t t = 0; t < 8; ++t)
    threads.emplace_back([t] {
      for (std::int32_t k = 0; k < 200; ++k)
        print(
            "thread %d line %d of a line long enough to be split if it could\n",
            t, k);
    });
  for (std::thread &t : threads)
    t.join();
  close_log();
  const std::string text = read_all(path);
  std::size_t lines = 0, at = 0;
  while (at < text.size()) {
    const std::size_t end = text.find('\n', at);
    const std::string line = text.substr(at, end - at);
    check::is_true(line.rfind("thread ", 0) == 0 &&
                       line.find("split if it could") != std::string::npos,
                   "a whole line: " + line);
    ++lines;
    at = end + 1;
  }
  check::equal(static_cast<std::int64_t>(lines), 1600, "every line");
}

} // namespace mxi
