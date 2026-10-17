// The Metal layer the GPU code shares: the device, an embedded shader library,
// a compute pipeline. The spot finder's threshold (dext_metal.cc) and profile
// fitting (fit_device_metal.cc) each did these for itself; they do them here.
// Errors are said, not thrown: each caller keeps its own policy -- the
// threshold refuses, profile fitting falls back to the CPU.
//
// This file's .cc is where metal-cpp's implementation is compiled
// (NS_PRIVATE_IMPLEMENTATION, MTL_PRIVATE_IMPLEMENTATION), once for any
// program: no other file may define those.

#ifndef MXI_METAL_DEVICE_HH
#define MXI_METAL_DEVICE_HH

#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>

#include <string>

namespace mxi::metal {

//: The system's default device, created once for the process and never
//: released; null where there is none.
MTL::Device *device();

//: A shader library from its embedded bytes (a .metallib), or null with
//: *error saying why.
MTL::Library *load_library(MTL::Device *device, const unsigned char *bytes,
                           unsigned long size, std::string *error);

//: A compute pipeline for the library's kernel `name`, taking at least
//: `min_threads` threads a threadgroup -- the limit is the pipeline's, falling
//: with register pressure, so it is checked rather than assumed -- or null with
//: *error saying why.
MTL::ComputePipelineState *pipeline(MTL::Device *device, MTL::Library *library,
                                    const char *name, NS::UInteger min_threads,
                                    std::string *error);

//: The device's name, as the system gives it.
std::string device_name(MTL::Device *device);

} // namespace mxi::metal

#endif // MXI_METAL_DEVICE_HH
