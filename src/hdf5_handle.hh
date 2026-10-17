// HDF5 handles for the whole library: one owner of an id and the function that
// closes it, and the one lock every HDF5 call that may meet another thread's is
// taken under.
//
// The image reader, the NXmx import and the .rflx code each had a class of
// their own for this -- Handle, H5 and H -- the reader's alone closing under
// the lock. HDF5 as packaged is not thread-safe: the image reader is called
// from several threads, and a close is a call like any other.

#ifndef MXI_HDF5_HANDLE_HH
#define MXI_HDF5_HANDLE_HH

#include <hdf5.h>

#include <mutex>

namespace mxi {

//: The lock around HDF5's API. Recursive, because the scopes nest: the reader
//: holds it across opening a dataset, and the dataset's temporary handles close
//: inside that.
std::recursive_mutex &hdf5_mutex();

//: An HDF5 id and its close (H5Fclose, H5Dclose, ...), closed once, under the
//: lock, when this goes. Moves; does not copy. An id below zero -- a failed
//: open -- is held but never closed: ok() says which.
class Hdf5Handle {
public:
  Hdf5Handle() = default;
  Hdf5Handle(hid_t id, herr_t (*close)(hid_t)) : id_(id), close_(close) {}
  ~Hdf5Handle() { reset(); }
  Hdf5Handle(Hdf5Handle &&other) noexcept
      : id_(other.id_), close_(other.close_) {
    other.id_ = -1;
  }
  Hdf5Handle &operator=(Hdf5Handle &&other) noexcept {
    if (this != &other) {
      reset();
      id_ = other.id_;
      close_ = other.close_;
      other.id_ = -1;
    }
    return *this;
  }
  Hdf5Handle(const Hdf5Handle &) = delete;
  Hdf5Handle &operator=(const Hdf5Handle &) = delete;

  hid_t get() const { return id_; }
  bool ok() const { return id_ >= 0; }

  void reset() {
    if (id_ >= 0 && close_ != nullptr) {
      const std::lock_guard<std::recursive_mutex> guard(hdf5_mutex());
      close_(id_);
    }
    id_ = -1;
  }

private:
  hid_t id_ = -1;
  herr_t (*close_)(hid_t) = nullptr;
};

} // namespace mxi

#endif // MXI_HDF5_HANDLE_HH
