// The shared CUDA layer: see cuda_device.hh.

#include "cuda_device.hh"

#include <cstdint>
#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <string>

namespace mxi::cuda {

bool available() {
  std::int32_t devices = 0;
  return cudaGetDeviceCount(&devices) == cudaSuccess && devices > 0;
}

std::string unavailable_reason() {
  std::int32_t devices = 0;
  const cudaError_t error = cudaGetDeviceCount(&devices);
  if (error == cudaSuccess && devices > 0)
    return std::string();
  std::string why = error != cudaSuccess
                        ? std::string(cudaGetErrorName(error)) + ", " +
                              cudaGetErrorString(error)
                        : std::string("CUDA found no devices");
  std::int32_t driver = 0, runtime = 0;
  cudaDriverGetVersion(&driver);
  cudaRuntimeGetVersion(&runtime);
  char versions[160];
  if (driver == 0)
    std::snprintf(versions, sizeof versions,
                  "; no CUDA driver is installed, and this was built "
                  "with the CUDA %d.%d runtime",
                  runtime / 1000, (runtime % 1000) / 10);
  else
    std::snprintf(versions, sizeof versions,
                  "; the driver supports CUDA %d.%d, and this was built "
                  "with the CUDA %d.%d runtime",
                  driver / 1000, (driver % 1000) / 10, runtime / 1000,
                  (runtime % 1000) / 10);
  why += versions;
  if (const char *visible = std::getenv("CUDA_VISIBLE_DEVICES"))
    why += std::string("; CUDA_VISIBLE_DEVICES is \"") + visible + "\"";
  // What the versions mean, which the bare numbers had left to the reader --
  // on a cluster node, a 13.4 driver beside this build's 12.8 runtime read
  // as the trouble, and was not: CUDA would start for no program there, the
  // GPU not the session's to use (Graeme, 14 November). A driver older than
  // the runtime is the one case the versions explain.
  if (driver > 0 && driver < runtime)
    why += ". The driver is older than this build's runtime: a newer NVIDIA "
           "driver is needed";
  else if (driver > 0)
    why += ". A driver as new as this build's runtime, or newer, is as it "
           "should be: the versions are not the trouble. Where nvidia-smi "
           "works but CUDA will not start, the GPU is not this process's to "
           "use -- on a cluster, run in a job given one (srun --gpus=1, or "
           "your scheduler's like); on a machine of one's own, "
           "/dev/nvidia-uvm missing (nvidia-modprobe -u -c=0, as root)";
  return why;
}

std::string device_name() {
  std::int32_t count = 0;
  cudaDeviceProp prop;
  if (cudaGetDeviceCount(&count) == cudaSuccess && count > 0 &&
      cudaGetDeviceProperties(&prop, 0) == cudaSuccess)
    return prop.name;
  return std::string();
}

} // namespace mxi::cuda
