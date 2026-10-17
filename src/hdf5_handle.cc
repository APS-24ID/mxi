// HDF5 handles and the lock around HDF5's API: see hdf5_handle.hh.
//
// Every call that may meet another thread's is serialised under this lock --
// the image reader's, and the closes in Hdf5Handle. Nothing is lost by doing it
// unconditionally: a thread-safe HDF5 takes one global lock across its whole
// public API anyway. What runs in parallel is the part that costs:
// decompression and analysis, downstream of the reads.

#include "hdf5_handle.hh"

namespace mxi {

std::recursive_mutex &hdf5_mutex() {
  static std::recursive_mutex mutex;
  return mutex;
}

} // namespace mxi
