// The allocator: where mimalloc is built in (MXI_MIMALLOC, Linux's default),
// malloc, new and the containers' memory are all its -- glibc's had cost
// integration 11 s of 38 on a two-socket Xeon (15 November). Elsewhere, the
// C library's, and nothing to check.

#include <cstdlib>
#include <memory>
#include <vector>

#ifdef MXI_MIMALLOC
#include <mimalloc.h>
#endif

#include "check.hh"

TEST(the_programs_memory_is_mimallocs_where_it_is_built_in) {
#ifdef MXI_MIMALLOC
  void *p = std::malloc(64);
  check::is_true(mi_is_in_heap_region(p), "malloc's memory mimalloc's");
  std::free(p);
  const auto q = std::make_unique<double[]>(1000);
  check::is_true(mi_is_in_heap_region(q.get()), "new's memory mimalloc's");
  std::vector<float> v(100000);
  check::is_true(mi_is_in_heap_region(v.data()),
                 "a vector's memory mimalloc's");
#else
  check::is_true(true, "the C library's malloc: nothing to check");
#endif
}
