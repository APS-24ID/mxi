// A test harness in a hundred lines, rather than a dependency.
//
// Registration is by static constructor, so a test file needs no list to be
// added to. Forgetting to register is the commonest way a test silently stops
// running, and this makes it impossible.
#pragma once

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace check {

struct Case {
  std::string name;
  std::function<void()> body;
};

inline std::vector<Case> &registry() {
  static std::vector<Case> cases;
  return cases;
}

struct Failure {
  std::string what;
};

inline void fail(const std::string &message) { throw Failure{message}; }

//: A path for a test's own file, in the system's temporary directory rather
//: than wherever the tests were started, which was once the repository.
inline std::string scratch_path(const std::string &name) {
  return (std::filesystem::temp_directory_path() / ("mxi_test_" + name))
      .string();
}

inline void is_true(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

inline void close(double a, double b, double tolerance,
                  const std::string &message) {
  if (!(std::abs(a - b) <= tolerance)) {
    char buffer[512];
    std::snprintf(buffer, sizeof(buffer), "%s: %.17g vs %.17g (tol %.3g)",
                  message.c_str(), a, b, tolerance);
    fail(buffer);
  }
}

inline void equal(long long a, long long b, const std::string &message) {
  if (a != b) {
    char buffer[512];
    std::snprintf(buffer, sizeof(buffer), "%s: %lld vs %lld", message.c_str(),
                  a, b);
    fail(buffer);
  }
}

//: A test that cannot run here says so and is counted as skipped: neither a
//: pass, which would claim something was checked, nor a failure, which would
//: blame the code for the machine. The four gigabyte test used to catch
//: bad_alloc and assert true, which reported "ok" for a check never made.
struct Skip {
  std::string why;
};
[[noreturn]] inline void skip(const std::string &why) { throw Skip{why}; }

struct Register {
  Register(const char *name, std::function<void()> body) {
    registry().push_back({name, std::move(body)});
  }
};

//: Every test, or with `only` those whose name contains it -- to run one
//: alone, as when chasing a test that hangs. Each is reported as it finishes,
//: flushed, so that a run stopped part way shows how far it got.
inline int run_all(const std::string &only = "") {
  int failed = 0;
  int skipped = 0;
  std::size_t ran = 0;
  for (const Case &c : registry()) {
    if (!only.empty() && c.name.find(only) == std::string::npos)
      continue;
    ++ran;
    try {
      c.body();
      std::printf("  ok    %s\n", c.name.c_str());
      std::fflush(stdout);
    } catch (const Skip &k) {
      std::printf("  skip  %s\n        %s\n", c.name.c_str(), k.why.c_str());
      ++skipped;
    } catch (const Failure &f) {
      std::printf("  FAIL  %s\n        %s\n", c.name.c_str(), f.what.c_str());
      ++failed;
    } catch (const std::exception &e) {
      std::printf("  ERROR %s\n        %s\n", c.name.c_str(), e.what());
      ++failed;
    }
  }
  if (skipped > 0) {
    std::printf("%zu tests, %d failed, %d skipped\n", ran, failed, skipped);
  } else {
    std::printf("%zu tests, %d failed\n", ran, failed);
  }
  return failed == 0 ? 0 : 1;
}

} // namespace check

#define TEST(name)                                                             \
  static void name();                                                          \
  static ::check::Register register_##name(#name, name);                       \
  static void name()
