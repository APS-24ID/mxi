// Profile fitting on a CUDA device: the batch's arrays to fit_cuda_submit and
// fit_cuda_collect, in fit_cuda.cu, which is all nvcc compiles.

#include "fit_batch.hh"

#include <cstdint>

namespace mxi::fitdev {
bool fit_cuda_destination(std::int32_t slot, std::size_t voxels,
                          std::size_t reference_floats, float **data,
                          std::uint8_t **mask, float **reference);
bool fit_cuda_submit(std::int32_t slot, const Setup &setup,
                     const PanelF *panels, std::int32_t n_panels,
                     const BoxF *boxes, std::int32_t n_boxes, const float *data,
                     const std::uint8_t *mask, std::size_t voxels,
                     const float *reference, std::size_t reference_floats,
                     std::size_t corner_floats, bool in_place);
bool fit_cuda_collect(std::int32_t slot, FitF *out);
const char *fit_cuda_name();
std::string fit_cuda_unavailable_reason();
void fit_cuda_times(double out[3]);
} // namespace mxi::fitdev

namespace mxi {

namespace {
std::int32_t g_next_ticket = 0;
std::size_t g_boxes[kFitSlots] = {0, 0};
} // namespace

bool fit_batch_destination(std::size_t voxels, std::size_t reference_floats,
                           FitDestination *out) {
  if (fit_device_name() == nullptr)
    return false;
  return mxi::fitdev::fit_cuda_destination(g_next_ticket % kFitSlots, voxels,
                                           reference_floats, &out->data,
                                           &out->mask, &out->reference);
}

std::int32_t fit_batch_submit(const FitBatch &batch) {
  if (batch.boxes.empty() || fit_device_name() == nullptr)
    return -1;
  const std::int32_t ticket = g_next_ticket;
  if (!mxi::fitdev::fit_cuda_submit(
          ticket % kFitSlots, batch.setup, batch.panels.data(),
          static_cast<std::int32_t>(batch.panels.size()), batch.boxes.data(),
          static_cast<std::int32_t>(batch.boxes.size()), batch.data_at,
          batch.mask_at, batch.voxels, batch.reference_at,
          batch.reference_floats, batch.corner_floats, batch.in_place))
    return -1;
  g_boxes[ticket % kFitSlots] = batch.boxes.size();
  ++g_next_ticket;
  return ticket;
}

bool fit_batch_collect(std::int32_t ticket,
                       std::vector<mxi::fitdev::FitF> *out) {
  out->assign(g_boxes[ticket % kFitSlots], mxi::fitdev::FitF{});
  return mxi::fitdev::fit_cuda_collect(ticket % kFitSlots, out->data());
}

bool fit_batch_device(const FitBatch &batch,
                      std::vector<mxi::fitdev::FitF> *out) {
  if (batch.boxes.empty()) {
    out->clear();
    return true;
  }
  const std::int32_t ticket = fit_batch_submit(batch);
  return ticket >= 0 && fit_batch_collect(ticket, out);
}

const char *fit_device_name() { return mxi::fitdev::fit_cuda_name(); }
std::string fit_device_unavailable_reason() {
  return mxi::fitdev::fit_cuda_unavailable_reason();
}
void fit_device_times(double out[3]) { mxi::fitdev::fit_cuda_times(out); }

} // namespace mxi
