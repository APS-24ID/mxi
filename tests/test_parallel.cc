#include <atomic>
#include <chrono>
#include <cstddef>
#include <stdexcept>
#include <thread>
#include <vector>

#include "../src/parallel.hh"
#include "check.hh"

namespace mxi {

namespace {

//: A sum over a million terms in fixed blocks, combined in block order.
double blocked_sum() {
  const std::size_t n = 1000000, blocks = kParallelBlocks;
  std::vector<double> part(blocks, 0.0);
  for_each_block(blocks, [&](std::size_t b) {
    double s = 0.0;
    for (std::size_t i = block_begin(n, blocks, b);
         i < block_begin(n, blocks, b + 1); ++i)
      s += 1.0 / (1.0 + static_cast<double>(i));
    part[b] = s;
  });
  double total = 0.0;
  for (double p : part)
    total += p;
  return total;
}

} // namespace

TEST(the_pool_gives_the_same_sum_every_call_and_on_any_count) {
  // Called hundreds of times, as the error model's search calls it: the same
  // bits each time, and the same on one thread as on eight.
  set_parallel_threads(1);
  const double one = blocked_sum();
  set_parallel_threads(8);
  for (int k = 0; k < 300; ++k)
    if (blocked_sum() != one) {
      set_parallel_threads(0);
      check::fail("call " + std::to_string(k) + " differed");
    }
  set_parallel_threads(3); // the pool is made again for a new count
  check::is_true(blocked_sum() == one,
                 "and on three threads, after a change of count");
  set_parallel_threads(0);
}

TEST(parallel_work_started_inside_a_block_runs_serially) {
  // A block that itself calls for_each_block must not wait on the pool it is
  // running in: it runs its inner work serially, and everything is done once.
  set_parallel_threads(4);
  std::atomic<int> count{0};
  for_each_block(8, [&](std::size_t) {
    for_each_block(8, [&](std::size_t) { count.fetch_add(1); });
  });
  set_parallel_threads(0);
  check::equal(static_cast<long long>(count.load()), 64,
               "eight inner blocks in each of eight");
}

TEST(an_exception_in_a_block_reaches_the_caller_and_the_pool_carries_on) {
  set_parallel_threads(4);
  bool caught = false;
  try {
    for_each_block(64, [&](std::size_t b) {
      if (b == 17)
        throw std::runtime_error("block 17");
    });
  } catch (const std::runtime_error &e) {
    caught = std::string(e.what()) == "block 17";
  }
  std::atomic<int> after{0};
  for_each_block(64, [&](std::size_t) { after.fetch_add(1); });
  set_parallel_threads(0);
  check::is_true(caught, "the caller gets the exception");
  check::equal(static_cast<long long>(after.load()), 64,
               "and the pool works after it");
}

} // namespace mxi

namespace mxi {

TEST(parallel_work_started_inside_a_block_finishes_whichever_thread_runs_it) {
  // A block may start parallel work of its own: it runs serially, rather than
  // waiting on the pool it is itself running in. Worker threads were marked
  // so; the calling thread, which takes blocks too, was not -- a nested loop
  // in a block it took locked the pool's one-job mutex it already held, and
  // waited on itself for ever. Only sometimes, as which thread takes which
  // block is the scheduler's; here every block takes 2 ms, so that the caller
  // surely takes some, and the test says it did.
  set_parallel_threads(4);
  const std::thread::id caller = std::this_thread::get_id();
  std::atomic<int> by_caller{0};
  std::vector<long long> sums(16, 0);
  for_each_index(16, [&](std::size_t i) {
    std::this_thread::sleep_for(std::chrono::milliseconds(2));
    if (std::this_thread::get_id() == caller)
      by_caller.fetch_add(1);
    std::vector<long long> inner(8, 0);
    for_each_index(8, [&](std::size_t j) {
      inner[j] = static_cast<long long>(i * 8 + j);
    });
    long long s = 0;
    for (long long v : inner)
      s += v;
    sums[i] = s;
  });
  set_parallel_threads(0);
  check::is_true(by_caller.load() > 0,
                 "the caller took a block with nested work");
  for (std::size_t i = 0; i < 16; ++i)
    check::equal(sums[i], static_cast<long long>(64 * i + 28),
                 "each block's own sum");
}

} // namespace mxi
