// The shared CUDA layer: see cuda_device.hh.

#include "cuda_device.hh"

#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <string>

namespace mxi::cuda {

bool available() {
  int devices = 0;
  return cudaGetDeviceCount(&devices) == cudaSuccess && devices > 0;
}

std::string unavailable_reason() {
  int devices = 0;
  const cudaError_t error = cudaGetDeviceCount(&devices);
  if (error == cudaSuccess && devices > 0)
    return std::string();
  std::string why = error != cudaSuccess
                        ? std::string(cudaGetErrorName(error)) + ", " +
                              cudaGetErrorString(error)
                        : std::string("CUDA found no devices");
  int driver = 0, runtime = 0;
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
  return why;
}

std::string device_name() {
  int count = 0;
  cudaDeviceProp prop;
  if (cudaGetDeviceCount(&count) == cudaSuccess && count > 0 &&
      cudaGetDeviceProperties(&prop, 0) == cudaSuccess)
    return prop.name;
  return std::string();
}

} // namespace mxi::cuda
