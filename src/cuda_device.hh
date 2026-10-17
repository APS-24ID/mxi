// The CUDA layer the GPU code shares: whether there is a device, why not, and
// its name. The spot finder's threshold (dext_cuda.cu) and profile fitting
// (fit_cuda.cu) each had these -- the reason word for word twice. Each keeps
// its own error policy for its kernels: the threshold throws, profile fitting
// falls back to the CPU.

#ifndef MXI_CUDA_DEVICE_HH
#define MXI_CUDA_DEVICE_HH

#include <string>

namespace mxi::cuda {

//: Whether a CUDA device is usable.
bool available();

//: Why no CUDA device is usable, in CUDA's own words, with what most often
//: explains it: the driver older than the runtime this was built with, or no
//: device visible to the process -- a login node, a job with no GPU, an empty
//: CUDA_VISIBLE_DEVICES. Empty when a device is usable.
std::string unavailable_reason();

//: The first device's name; empty where there is none.
std::string device_name();

} // namespace mxi::cuda

#endif // MXI_CUDA_DEVICE_HH
