// The shared Metal layer: see metal_device.hh. metal-cpp's implementation is
// compiled here, and only here.

#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include "metal_device.hh"

#include <dispatch/dispatch.h>

namespace mxi::metal {

namespace {

std::string described(NS::Error *error) {
  if (error != nullptr && error->localizedDescription() != nullptr)
    return error->localizedDescription()->utf8String();
  return "no reason given";
}

} // namespace

MTL::Device *device() {
  static MTL::Device *const it = MTL::CreateSystemDefaultDevice();
  return it;
}

MTL::Library *load_library(MTL::Device *device, const unsigned char *bytes,
                           unsigned long size, std::string *error) {
  if (device == nullptr) {
    if (error)
      *error = "no Metal device";
    return nullptr;
  }
  dispatch_data_t data = dispatch_data_create(bytes, size, nullptr,
                                              DISPATCH_DATA_DESTRUCTOR_DEFAULT);
  NS::Error *why = nullptr;
  MTL::Library *library = device->newLibrary(data, &why);
  dispatch_release(data);
  if (library == nullptr && error)
    *error = described(why);
  return library;
}

MTL::ComputePipelineState *pipeline(MTL::Device *device, MTL::Library *library,
                                    const char *name, NS::UInteger min_threads,
                                    std::string *error) {
  MTL::Function *function =
      library->newFunction(NS::String::string(name, NS::UTF8StringEncoding));
  if (function == nullptr) {
    if (error)
      *error = std::string("no kernel named ") + name;
    return nullptr;
  }
  NS::Error *why = nullptr;
  MTL::ComputePipelineState *state =
      device->newComputePipelineState(function, &why);
  function->release();
  if (state == nullptr) {
    if (error)
      *error = std::string(name) + ": " + described(why);
    return nullptr;
  }
  const NS::UInteger limit = state->maxTotalThreadsPerThreadgroup();
  if (limit < min_threads) {
    if (error)
      *error = std::string(name) + " will only take " + std::to_string(limit) +
               " threads per threadgroup, and needs " +
               std::to_string(min_threads);
    state->release();
    return nullptr;
  }
  return state;
}

std::string device_name(MTL::Device *device) {
  return device != nullptr ? device->name()->utf8String() : std::string();
}

} // namespace mxi::metal
