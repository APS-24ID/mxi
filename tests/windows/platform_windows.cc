// The platform layer's Windows half, run: reads at offsets from many
// threads, the end of a file, peak memory, and a log mirrored -- built by
// MinGW-w64 and run under Wine by scripts/check_windows.sh, where neither a
// Windows machine nor its build job is (docs/windows.md).
#include "log_mirror.hh"
#include "output.hh"
#include "platform.hh"
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>
int main() {
  std::int32_t failures = 0;
  const auto check = [&](bool ok, const char *what) {
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok)
      ++failures;
  };
  // A file of known bytes: byte i is i mod 251.
  const std::size_t n = 3000000;
  {
    std::FILE *f = std::fopen("bytes.bin", "wb");
    for (std::size_t i = 0; i < n; ++i)
      std::fputc(static_cast<int>(i % 251), f);
    std::fclose(f);
  }
  const mxi::platform::File file = mxi::platform::open_for_reading("bytes.bin");
  check(file != mxi::platform::kNoFile, "opened");
  std::vector<std::thread> threads;
  std::vector<std::int32_t> bad(4, 0);
  for (std::int32_t t = 0; t < 4; ++t)
    threads.emplace_back([&, t] {
      for (std::int32_t k = 0; k < 500; ++k) {
        const std::uint64_t at =
            (static_cast<std::uint64_t>(t) * 7919u + k * 5003u) % (n - 4096);
        unsigned char buffer[4096];
        const std::int64_t got =
            mxi::platform::read_at(file, buffer, sizeof buffer, at);
        if (got != 4096) {
          ++bad[t];
          continue;
        }
        for (std::size_t i = 0; i < sizeof buffer; ++i)
          if (buffer[i] != (at + i) % 251) {
            ++bad[t];
            break;
          }
      }
    });
  for (std::thread &t : threads)
    t.join();
  check(bad[0] + bad[1] + bad[2] + bad[3] == 0,
        "2000 reads at offsets from 4 threads, every byte right");
  unsigned char tail[100];
  check(mxi::platform::read_at(file, tail, sizeof tail, n - 10) == 10,
        "10 bytes before the end: 10");
  check(mxi::platform::read_at(file, tail, sizeof tail, n) == 0,
        "at the end: 0");
  check(mxi::platform::open_for_reading("no such file") ==
            mxi::platform::kNoFile,
        "no file: kNoFile");
  check(mxi::platform::peak_memory_mb() > 0.0, "peak memory above 0");
  mxi::platform::prepare_terminal(); // no console here: nothing to change
  check(true, "prepare_terminal returns, with or without a console");
  std::printf("     peak %.1f MB, columns %d, terminal %d\n",
              mxi::platform::peak_memory_mb(),
              static_cast<int>(mxi::platform::stderr_columns()),
              mxi::platform::stderr_is_terminal());
  mxi::mirror_to_log("win.log");
  mxi::print("to the log %zu %" PRId64 "\n", static_cast<std::size_t>(42),
             static_cast<std::int64_t>(-7));
  mxi::close_log();
  std::FILE *f = std::fopen("win.log", "rb");
  char text[64] = {0};
  std::fread(text, 1, sizeof text - 1, f);
  std::fclose(f);
  check(std::string(text) == "to the log 42 -7\n",
        "the log holds what was printed, %zu and PRId64 right");
  std::printf("%d failures\n", static_cast<int>(failures));
  return failures;
}
