// mxi::series::nxmx -- read a series out of an NXmx HDF5 file.
//
// The master file's /entry/data/data is a virtual dataset spanning a set of
// data files. A virtual dataset cannot be chunk-read, so the VDS mapping is
// unpacked instead: it names its own source files, datasets and offsets, which
// means nothing here has to assume how the writer laid the series out. Frames
// are then read with H5Dread_chunk straight out of the data file that holds
// them, so the compressed bytes reach the threshold without HDF5's filter
// pipeline ever running -- the same path the data took in, in reverse.

#include <fcntl.h>
#include <hdf5.h>
#include <unistd.h>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hdf5_handle.hh"
#include "image_series.hh"

namespace mxi::series {
namespace {

constexpr const char *kImageData = "/entry/data/data";

// HDF5 filter identifiers, as registered with The HDF Group.
constexpr H5Z_filter_t kBitshuffle = 32008;
constexpr H5Z_filter_t kLz4 = 32004;

// The HDF5 C library keeps global state -- identifier tables, free lists, the
// metadata cache -- and is not thread-safe unless it was built that way, which
// vcpkg's default port is not and several distributions are. Giving each thread
// its own file handle is not enough: two threads inside the library at once
// will corrupt something and crash, which is exactly what happens on a build
// without thread safety.
//
// So every call below is serialised, including the closes in Handle. Nothing is
// lost by doing it unconditionally, because a thread-safe HDF5 takes one global
// lock across its whole public API anyway -- the reads were already serialised
// on those builds. What still runs in parallel is the part that costs:
// decompression and analysis, downstream of here.
//
// Recursive because the scopes nest -- read() holds it across open(), and
// open() holds it while its temporary handles close.
using mxi::hdf5_mutex;

using Guard = std::lock_guard<std::recursive_mutex>;

void ok(herr_t status, const char *what) {
  if (status < 0)
    throw std::runtime_error(std::string("HDF5: ") + what + " failed");
}

// HDF5 2.0 added a buffer-size argument to H5Dread_chunk -- in, the buffer's
// capacity; out, what the chunk needed -- and its versioned-API macros point
// the old name at the new function, so the call cannot be written once for
// both. Which shape this HDF5 has is settled by a compile probe in
// CMakeLists.txt, since vcpkg and the distributions are currently on different
// sides of it.
herr_t hdf5_read_chunk(hid_t dataset, const hsize_t *offset,
                       std::uint32_t *filters, void *buffer,
                       std::size_t capacity) {
#ifdef MXI_H5DREAD_CHUNK_TAKES_SIZE
  std::size_t size = capacity;
  return H5Dread_chunk(dataset, H5P_DEFAULT, offset, filters, buffer, &size);
#else
  (void)capacity;
  return H5Dread_chunk(dataset, H5P_DEFAULT, offset, filters, buffer);
#endif
}

hid_t opened(hid_t id, const char *what) {
  if (id < 0)
    throw std::runtime_error(std::string("HDF5: ") + what + " failed");
  return id;
}

// HDF5 hands out integer handles with a different close function for each kind
// of object, so the closer travels with the handle.
using Handle = mxi::Hdf5Handle;

// Where one run of frames lives. One of these per data file, in the usual case
// of a VDS built block by block.
struct Block {
  std::string filename; // resolved against the master file's directory
  std::string dataset;
  std::uint64_t first = 0;  // first frame index in the series
  std::uint64_t last = 0;   // last frame index, inclusive
  std::uint64_t offset = 0; // index of `first` within the data file
};

// The two strings that identify a source dataset, so it can be opened once per
// thread rather than once per frame.
std::string cache_key(const Block &block) {
  return block.filename + "\n" + block.dataset;
}

std::string plist_string(hid_t plist, int index,
                         ssize_t (*get)(hid_t, size_t, char *, size_t)) {
  const ssize_t size = get(plist, static_cast<size_t>(index), nullptr, 0);
  if (size < 0)
    throw std::runtime_error("HDF5: cannot size a virtual mapping name");
  // The size HDF5 reports excludes the terminator but the buffer size it wants
  // includes it, so the string is allocated one longer and trimmed after.
  std::string text(static_cast<std::size_t>(size) + 1, '\0');
  if (get(plist, static_cast<size_t>(index), text.data(), text.size()) < 0) {
    throw std::runtime_error("HDF5: cannot read a virtual mapping name");
  }
  text.resize(static_cast<std::size_t>(size));
  return text;
}

// Which compression the dataset declares. Only the identifier matters: the
// bytes are handed on untouched, so nothing here needs the filter to be
// installed.
std::string_view algorithm_of(hid_t dcpl) {
  const int filters = H5Pget_nfilters(dcpl);
  if (filters < 0)
    throw std::runtime_error("HDF5: cannot count filters");
  for (int i = 0; i < filters; i++) {
    unsigned flags = 0;
    size_t values = 0;
    const H5Z_filter_t filter =
        H5Pget_filter2(dcpl, static_cast<unsigned>(i), &flags, &values, nullptr,
                       0, nullptr, nullptr);
    if (filter == kBitshuffle)
      return "bslz4";
    if (filter == kLz4)
      return "lz4";
  }
  return {};
}

// A data file whose chunks read_chunk can read directly: its descriptor and
// where HDF5's addresses count from (the end of its user block). fd -1 when
// it cannot: the file not on HDF5's default driver, which alone keeps a
// chunk's bytes at its address in the one file named; the file not opened;
// or MXI_READ_CHUNK=hdf5 in the environment, which sends every read through
// H5Dread_chunk, to compare. One descriptor a file for the whole process,
// shared by every reader -- pread names its offset -- since one a reader took
// threads x files of them, past macOS's 256 a process (Graeme). Called with
// the HDF5 lock held.
struct DirectFile {
  int fd = -1;
  std::uint64_t base = 0;
};

struct SourceDataset {
  Handle file;
  Handle dataset;
  DirectFile direct; //: where read_chunk reads its chunks, found once
  std::uint64_t height = 0;
  std::uint64_t width = 0;
  unsigned bit_depth = 0;
  std::string_view algorithm;
};

DirectFile direct_file_of(hid_t dataset) {
  static const bool through_hdf5 = [] {
    const char *how = std::getenv("MXI_READ_CHUNK");
    return how != nullptr && std::string(how) == "hdf5";
  }();
  if (through_hdf5)
    return {};
  static std::map<std::string, DirectFile> files;
  const hid_t file = H5Iget_file_id(dataset);
  if (file < 0)
    return {};
  const ssize_t n = H5Fget_name(file, nullptr, 0);
  std::string name(n > 0 ? static_cast<std::size_t>(n) : 0, '\0');
  if (n > 0)
    H5Fget_name(file, name.data(), static_cast<std::size_t>(n) + 1);
  const auto found = files.find(name);
  if (found != files.end()) {
    H5Fclose(file);
    return found->second;
  }
  DirectFile direct;
  const hid_t fapl = H5Fget_access_plist(file);
  const bool default_driver = fapl >= 0 && H5Pget_driver(fapl) == H5FD_SEC2;
  if (fapl >= 0)
    H5Pclose(fapl);
  if (default_driver && !name.empty()) {
    const hid_t fcpl = H5Fget_create_plist(file);
    hsize_t userblock = 0;
    if (fcpl >= 0) {
      H5Pget_userblock(fcpl, &userblock);
      H5Pclose(fcpl);
    }
    direct.base = userblock;
    direct.fd = ::open(name.c_str(), O_RDONLY);
  }
  H5Fclose(file);
  files.emplace(name, direct);
  return direct;
}

// H5Dread_chunk: the same arguments and the same result -- a chunk's stored
// bytes, still compressed, and the mask of the filters skipped for it -- but
// the bytes read with pread, outside HDF5's global lock, where HDF5 itself
// says they are: H5Dget_chunk_info_by_coord gives the chunk's address and
// size, and on HDF5's default driver those bytes lie there, contiguous, in the
// file holding the dataset (through an external link, the data file, not the
// master). Only that lookup is under the lock. Inside it, 32 threads queued
// behind one -- 1000 thread-seconds of fetching in 67 s of reading. (On one
// thread the two cost the same, warm: 2.2 and 2.6 s on insulin; a 'five times'
// once written here was a cold file cache.)
// Where a direct read is not sound (direct_file_of), H5Dread_chunk itself.
// tests/test_read_chunk.py checks the two give the same spots.
//: Where a chunk is, as H5Dget_chunk_info_by_coord says: given to read_chunk
//: by a caller that has already asked, so that it is not asked twice under
//: the lock (ferritin's reading 6.4 s to 8.7 with the second asking).
struct ChunkPlace {
  haddr_t address = HADDR_UNDEF;
  hsize_t size = 0;
  unsigned mask = 0;
};

herr_t read_chunk(hid_t dataset, const hsize_t *offset, std::uint32_t *filters,
                  void *buffer, std::size_t capacity,
                  const DirectFile *known = nullptr,
                  const ChunkPlace *place = nullptr) {
  DirectFile file;
  unsigned mask = 0;
  haddr_t address = 0;
  hsize_t size = 0;
  if (known != nullptr && known->fd >= 0 && place != nullptr) {
    // Both known: no lock at all.
    file = *known;
    mask = place->mask;
    address = place->address;
    size = place->size;
  } else {
    const Guard guard(hdf5_mutex());
    // The file the dataset is in: given by a caller that knows it, which saves
    // asking HDF5 its name on every read -- under the lock, 16 threads queued
    // on that, ferritin's reading 6.4 s to 16.4 (Graeme).
    file = known != nullptr ? *known : direct_file_of(dataset);
    if (file.fd < 0)
      return hdf5_read_chunk(dataset, offset, filters, buffer, capacity);
    if (H5Dget_chunk_info_by_coord(dataset, offset, &mask, &address, &size) < 0)
      return -1;
  }
  if (address == HADDR_UNDEF || size > capacity)
    return -1;
  auto *into = static_cast<char *>(buffer);
  std::size_t done = 0;
  const auto at = static_cast<off_t>(file.base + address);
  while (done < size) {
    const ssize_t got = ::pread(file.fd, into + done, size - done,
                                at + static_cast<off_t>(done));
    if (got <= 0)
      return -1;
    done += static_cast<std::size_t>(got);
  }
  *filters = mask;
  return 0;
}

class NxmxReader : public Reader {
public:
  NxmxReader(std::vector<Block> blocks,
             std::shared_ptr<const std::vector<std::uint64_t>> masked)
      : blocks_(std::move(blocks)), masked_(std::move(masked)) {}

  bool read(const std::string &key, Frame *frame) override {
    frame->masked = masked_;
    const std::uint64_t index = std::strtoull(key.c_str(), nullptr, 10);
    SourceDataset *source = nullptr;
    hsize_t offset[3] = {0, 0, 0};
    ChunkPlace place;
    {
      const Guard guard(hdf5_mutex());
      const Block &block = locate(index);
      source = &open(block);
      // The chunk holds one frame, so its coordinate is the frame's own index
      // within the data file, and the other two are zero by construction.
      offset[0] = static_cast<hsize_t>(index - block.first + block.offset);
      ok(H5Dget_chunk_info_by_coord(source->dataset.get(), offset, &place.mask,
                                    &place.address, &place.size),
         "H5Dget_chunk_info_by_coord");
      // A chunk that was never allocated is a frame the writer never received.
      // It reads back as the fill value, which is not something to analyse.
      if (place.address == HADDR_UNDEF || place.size == 0)
        return false;
      frame->storage.resize(static_cast<std::size_t>(place.size));
    }
    // Outside the lock: read_chunk takes it only to look the chunk up.
    std::uint32_t filters = 0;
    ok(read_chunk(source->dataset.get(), offset, &filters,
                  frame->storage.data(), frame->storage.size(), &source->direct,
                  &place),
       "read_chunk");
    frame->height = source->height;
    frame->width = source->width;
    frame->bit_depth = source->bit_depth;
    // Bit i of the mask means filter i was skipped for this chunk, which is how
    // a writer stores a payload it could not compress. The algorithm string
    // belongs to the source dataset, which this reader owns and which outlives
    // every frame it fills.
    frame->algorithm = (filters & 1u) ? std::string_view() : source->algorithm;
    frame->data = {
        reinterpret_cast<const std::uint8_t *>(frame->storage.data()),
        frame->storage.size()};
    frame->number = static_cast<std::int64_t>(index);
    return true;
  }

private:
  const Block &locate(std::uint64_t index) const {
    for (const Block &block : blocks_) {
      if (index >= block.first && index <= block.last)
        return block;
    }
    throw std::runtime_error("frame " + std::to_string(index) +
                             " is outside every virtual mapping");
  }

  // Why a dataset could not be opened, when it is an external link to a file
  // that is not there -- the commonest way: a data file renamed or left behind,
  // as an upload that added "_1" to a name already taken did. Empty otherwise.
  // Called with the HDF5 lock held, its error printing off.
  std::string missing_link_target(hid_t file, const std::string &name,
                                  const std::string &file_path) {
    if (H5Lexists(file, name.c_str(), H5P_DEFAULT) <= 0)
      return "";
#if H5_VERSION_GE(1, 12, 0)
    H5L_info2_t link;
    if (H5Lget_info2(file, name.c_str(), &link, H5P_DEFAULT) < 0)
      return "";
#else
    H5L_info_t link;
    if (H5Lget_info(file, name.c_str(), &link, H5P_DEFAULT) < 0)
      return "";
#endif
    if (link.type != H5L_TYPE_EXTERNAL)
      return "";
    std::vector<char> value(link.u.val_size);
    if (H5Lget_val(file, name.c_str(), value.data(), value.size(),
                   H5P_DEFAULT) < 0)
      return "";
    const char *target = nullptr, *object = nullptr;
    unsigned flags = 0;
    if (H5Lunpack_elink_val(value.data(), value.size(), &flags, &target,
                            &object) < 0)
      return "";
    std::filesystem::path where(target);
    if (where.is_relative())
      where = std::filesystem::path(file_path).parent_path() / where;
    std::error_code ignored;
    if (std::filesystem::exists(where, ignored))
      return "";
    return std::string(": it is a link to ") + object + " in " + target +
           ", and there is no " + where.string();
  }

  // Called with the HDF5 lock already held, which the recursive mutex allows.
  SourceDataset &open(const Block &block) {
    const Guard guard(hdf5_mutex());
    const std::string key = cache_key(block);
    const auto found = open_.find(key);
    if (found != open_.end())
      return found->second;

    SourceDataset source;
    source.file = Handle(
        opened(H5Fopen(block.filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT),
               ("opening " + block.filename).c_str()),
        H5Fclose);
    {
      // HDF5's own error stack not printed: what went wrong is said here.
      H5E_auto2_t saved_func = nullptr;
      void *saved_data = nullptr;
      H5Eget_auto2(H5E_DEFAULT, &saved_func, &saved_data);
      H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
      const hid_t id =
          H5Dopen2(source.file.get(), block.dataset.c_str(), H5P_DEFAULT);
      const std::string why =
          id < 0 ? missing_link_target(source.file.get(), block.dataset,
                                       block.filename)
                 : std::string();
      H5Eset_auto2(H5E_DEFAULT, saved_func, saved_data);
      if (id < 0)
        throw std::runtime_error("HDF5: opening " + block.dataset + " failed" +
                                 why);
      source.dataset = Handle(id, H5Dclose);
    }
    source.direct = direct_file_of(source.dataset.get());

    const Handle dcpl(opened(H5Dget_create_plist(source.dataset.get()),
                             "H5Dget_create_plist"),
                      H5Pclose);
    if (H5Pget_layout(dcpl.get()) != H5D_CHUNKED) {
      throw std::runtime_error(block.filename + " " + block.dataset +
                               " is not chunked, so it cannot be chunk-read");
    }
    hsize_t chunk[3] = {0, 0, 0};
    if (H5Pget_chunk(dcpl.get(), 3, chunk) != 3) {
      throw std::runtime_error("expected a three-dimensional chunk in " +
                               block.filename);
    }
    if (chunk[0] != 1) {
      throw std::runtime_error(
          "chunks in " + block.filename + " span " + std::to_string(chunk[0]) +
          " frames; only one frame per chunk can be read this way");
    }
    source.height = chunk[1];
    source.width = chunk[2];
    source.algorithm = algorithm_of(dcpl.get());

    const Handle type(opened(H5Dget_type(source.dataset.get()), "H5Dget_type"),
                      H5Tclose);
    source.bit_depth = static_cast<unsigned>(8 * H5Tget_size(type.get()));

    return open_.emplace(key, std::move(source)).first->second;
  }

  std::vector<Block> blocks_;
  std::shared_ptr<const std::vector<std::uint64_t>> masked_;
  std::map<std::string, SourceDataset> open_;
};

//: Whether a file carries HDF5's signature where a superblock may be: at the
//: start, or after a user block at 512 bytes or any power of two beyond.
bool looks_like_hdf5(const std::string &path) {
  static const unsigned char signature[8] = {0x89, 'H',  'D',  'F',
                                             '\r', '\n', 0x1a, '\n'};
  std::FILE *file = std::fopen(path.c_str(), "rb");
  if (file == nullptr)
    return false;
  bool found = false;
  for (long at = 0; !found && at <= (1L << 30); at = at == 0 ? 512 : at * 2) {
    unsigned char got[8];
    if (std::fseek(file, at, SEEK_SET) != 0 || std::fread(got, 1, 8, file) != 8)
      break;
    found = std::memcmp(got, signature, 8) == 0;
  }
  std::fclose(file);
  return found;
}

//: Whether a file's first character that is not white space is a brace.
bool starts_with_brace(const std::string &path) {
  std::FILE *file = std::fopen(path.c_str(), "rb");
  if (file == nullptr)
    return false;
  int c = 0;
  while ((c = std::fgetc(file)) != EOF && std::isspace(c))
    ;
  std::fclose(file);
  return c == '{';
}

// The detector's own bad pixels, as dials.import reads them: NXdetector's
// pixel_mask, a pixel bad wherever it is not zero; or, in older DECTRIS
// masters, the same under detectorSpecific. None if neither is there, nor if
// it cannot be read -- an external link to a file not beside the master -- nor
// if its shape is not the frame's, said so. Called with the HDF5 lock held.
std::shared_ptr<const std::vector<std::uint64_t>>
read_pixel_mask(hid_t file, std::uint64_t height, std::uint64_t width) {
  // Probing paths that may not be there: HDF5 would print its error stack.
  H5E_auto2_t saved_func = nullptr;
  void *saved_data = nullptr;
  H5Eget_auto2(H5E_DEFAULT, &saved_func, &saved_data);
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  struct Restore {
    H5E_auto2_t func;
    void *data;
    ~Restore() { H5Eset_auto2(H5E_DEFAULT, func, data); }
  } restore{saved_func, saved_data};
  const auto exists = [&](const std::string &path) {
    std::string so_far;
    std::size_t at = 1;
    while (at <= path.size()) {
      const std::size_t slash = path.find('/', at);
      so_far = path.substr(0, slash == std::string::npos ? path.size() : slash);
      if (H5Lexists(file, so_far.c_str(), H5P_DEFAULT) <= 0)
        return false;
      if (slash == std::string::npos)
        break;
      at = slash + 1;
    }
    return true;
  };
  for (const char *path :
       {"/entry/instrument/detector/pixel_mask",
        "/entry/instrument/detector/detectorSpecific/pixel_mask"}) {
    if (!exists(path))
      continue;
    const hid_t dataset = H5Dopen2(file, path, H5P_DEFAULT);
    if (dataset < 0)
      continue;
    const hid_t space = H5Dget_space(dataset);
    hsize_t dims[2] = {0, 0};
    const int rank = space >= 0 ? H5Sget_simple_extent_ndims(space) : -1;
    if (rank == 2)
      H5Sget_simple_extent_dims(space, dims, nullptr);
    if (space >= 0)
      H5Sclose(space);
    if (rank != 2 || dims[0] != height || dims[1] != width) {
      std::fprintf(stderr,
                   "the detector's pixel_mask at %s is not the frame's shape, "
                   "%llu by %llu: not used\n",
                   path, static_cast<unsigned long long>(height),
                   static_cast<unsigned long long>(width));
      H5Dclose(dataset);
      continue;
    }
    std::vector<std::uint32_t> values(height * width);
    const herr_t read = H5Dread(dataset, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, values.data());
    H5Dclose(dataset);
    if (read < 0)
      continue;
    auto masked = std::make_shared<std::vector<std::uint64_t>>();
    for (std::size_t i = 0; i < values.size(); ++i)
      if (values[i] != 0)
        masked->push_back(i);
    return masked;
  }
  return nullptr;
}

class Nxmx : public Series {
public:
  const std::vector<Block> &blocks() const { return blocks_; }
  explicit Nxmx(std::string master) : master_(std::move(master)) {}

  bool try_open(Info *info) override {
    std::error_code ignored;
    if (!std::filesystem::is_regular_file(master_, ignored))
      return false;

    // Not HDF5 at all -- an experiment list given in the images' place, say --
    // said plainly, rather than in HDF5's diagnostic stack, which a file that
    // is HDF5 and fails for another reason still gets.
    if (!looks_like_hdf5(master_))
      throw std::runtime_error(
          master_ + " is not an HDF5 file" +
          (starts_with_brace(master_)
               ? std::string(": it looks like JSON, an experiment list "
                             "perhaps, which is not the images")
               : std::string()));

    const Guard guard(hdf5_mutex());

    const Handle file(
        opened(H5Fopen(master_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT),
               ("opening " + master_).c_str()),
        H5Fclose);
    // The older DECTRIS file writer's master has no /entry/data/data, only
    // /entry/data/data_000001 and on, each a link to a data file's own stack.
    if (H5Lexists(file.get(), kImageData, H5P_DEFAULT) <= 0 &&
        H5Lexists(file.get(), "/entry/data/data_000001", H5P_DEFAULT) > 0) {
      blocks_ = linked_blocks(file.get(), info);
      info->name = std::filesystem::path(master_).stem().string();
      masked_ = read_pixel_mask(file.get(), info->height, info->width);
      info->masked_pixels = masked_ ? masked_->size() : 0;
      return true;
    }
    const Handle data(opened(H5Dopen2(file.get(), kImageData, H5P_DEFAULT),
                             ("opening " + std::string(kImageData)).c_str()),
                      H5Dclose);

    const Handle space(opened(H5Dget_space(data.get()), "H5Dget_space"),
                       H5Sclose);
    if (H5Sget_simple_extent_ndims(space.get()) != 3) {
      throw std::runtime_error(std::string(kImageData) +
                               " is not a stack of images");
    }
    hsize_t dims[3] = {0, 0, 0};
    ok(H5Sget_simple_extent_dims(space.get(), dims, nullptr),
       "H5Sget_simple_extent_dims");

    const Handle dcpl(
        opened(H5Dget_create_plist(data.get()), "H5Dget_create_plist"),
        H5Pclose);
    blocks_ = H5Pget_layout(dcpl.get()) == H5D_VIRTUAL
                  ? unpack_virtual(dcpl.get(), dims[0])
                  : whole_dataset(dims[0]);

    info->name = std::filesystem::path(master_).stem().string();
    info->images = dims[0];
    info->height = dims[1];
    info->width = dims[2];
    masked_ = read_pixel_mask(file.get(), dims[1], dims[2]);
    info->masked_pixels = masked_ ? masked_->size() : 0;
    return true;
  }

  // Every frame is there to be read as soon as the file opens: there is nothing
  // to wait for once the master file exists.
  bool restrict_frames(std::uint64_t first, std::uint64_t last) override {
    first_frame_ = first;
    last_frame_ = last;
    return true;
  }

  std::vector<std::string> ready() override {
    std::vector<std::string> keys;
    if (dispatched_)
      return keys;
    dispatched_ = true;
    keys.reserve(static_cast<std::size_t>(images_of_blocks()));
    for (std::uint64_t i = 0; i < images_of_blocks(); i++) {
      // Only the frames the scan covers, when one was given. The key is the
      // array index, so the range is a filter on it and costs nothing.
      if (i < first_frame_ || i > last_frame_)
        continue;
      char key[24];
      std::snprintf(key, sizeof(key), "%06llu",
                    static_cast<unsigned long long>(i));
      keys.emplace_back(key);
    }
    return keys;
  }

  bool finished() override { return true; }

  std::string describe() const override { return master_; }

  std::unique_ptr<Reader> reader() override {
    return std::unique_ptr<Reader>(new NxmxReader(blocks_, masked_));
  }

private:
  std::uint64_t images_of_blocks() const {
    std::uint64_t images = 0;
    for (const Block &block : blocks_)
      images = std::max(images, block.last + 1);
    return images;
  }

  // Read the mapping out of the dataset creation property list. Each mapping's
  // selection in the virtual dataset gives the range of frames it covers, and
  // the matching selection in the source gives where they start in the data
  // file; that is all that is needed to turn a frame number into a chunk
  // coordinate.
  // Called from try_open with the HDF5 lock held.
  std::vector<Block> unpack_virtual(hid_t dcpl, std::uint64_t images) {
    size_t count = 0;
    ok(H5Pget_virtual_count(dcpl, &count), "H5Pget_virtual_count");
    if (count == 0)
      throw std::runtime_error(std::string(kImageData) + " maps to nothing");

    const std::filesystem::path directory =
        std::filesystem::path(master_).parent_path();

    std::vector<Block> blocks;
    for (size_t i = 0; i < count; i++) {
      Block block;
      const std::string filename =
          plist_string(dcpl, static_cast<int>(i), H5Pget_virtual_filename);
      block.dataset =
          plist_string(dcpl, static_cast<int>(i), H5Pget_virtual_dsetname);

      // "." means the mapping points back into the master file itself, and
      // the master's path is already whatever the caller gave -- it must NOT
      // be resolved against its own directory a second time. Doing so turned
      // "../ins10_1.nxs" into "../../ins10_1.nxs" and failed to open, for any
      // master named by a relative path with a directory in it.
      //
      // Only a genuinely external filename, which HDF5 stores relative to the
      // master, needs the directory prepended.
      std::filesystem::path path;
      if (filename == ".") {
        path = std::filesystem::path(master_);
      } else {
        path = std::filesystem::path(filename);
        if (path.is_relative() && !directory.empty())
          path = directory / path;
      }
      block.filename = path.string();

      const Handle vspace(
          opened(H5Pget_virtual_vspace(dcpl, static_cast<unsigned>(i)),
                 "H5Pget_virtual_vspace"),
          H5Sclose);
      const Handle source(
          opened(H5Pget_virtual_srcspace(dcpl, static_cast<unsigned>(i)),
                 "H5Pget_virtual_srcspace"),
          H5Sclose);

      hsize_t start[3] = {0, 0, 0};
      hsize_t end[3] = {0, 0, 0};
      ok(H5Sget_select_bounds(vspace.get(), start, end),
         "H5Sget_select_bounds on a virtual selection");
      block.first = start[0];
      // An unlimited selection reports an enormous upper bound, so the
      // dataset's own extent is the authority on where the series stops.
      block.last = std::min<std::uint64_t>(end[0], images - 1);

      hsize_t source_start[3] = {0, 0, 0};
      hsize_t source_end[3] = {0, 0, 0};
      ok(H5Sget_select_bounds(source.get(), source_start, source_end),
         "H5Sget_select_bounds on a source selection");
      block.offset = source_start[0];

      blocks.push_back(block);
    }
    std::sort(blocks.begin(), blocks.end(),
              [](const Block &a, const Block &b) { return a.first < b.first; });
    return blocks;
  }

  // A block a link, /entry/data/data_000001 and on, in order: each link's file
  // resolved against the master's directory, as a virtual dataset's sources
  // are, and its frames following the last's. A link whose file is not there
  // is refused by name, rather than the series taken to end before it.
  // Called from try_open with the HDF5 lock held.
  std::vector<Block> linked_blocks(hid_t file, Info *info) {
    const std::filesystem::path directory =
        std::filesystem::path(master_).parent_path();
    std::vector<Block> blocks;
    std::uint64_t next = 0;
    for (int k = 1;; ++k) {
      char name[64];
      std::snprintf(name, sizeof name, "/entry/data/data_%06d", k);
      if (H5Lexists(file, name, H5P_DEFAULT) <= 0)
        break;
        // H5Lget_info2 from HDF5 1.12; before, H5Lget_info, with the same
        // fields this needs.
#if H5_VERSION_GE(1, 12, 0)
      H5L_info2_t link;
      ok(H5Lget_info2(file, name, &link, H5P_DEFAULT), "H5Lget_info2");
#else
      H5L_info_t link;
      ok(H5Lget_info(file, name, &link, H5P_DEFAULT), "H5Lget_info");
#endif
      Block block;
      if (link.type == H5L_TYPE_EXTERNAL) {
        std::vector<char> value(link.u.val_size);
        ok(H5Lget_val(file, name, value.data(), value.size(), H5P_DEFAULT),
           "H5Lget_val");
        const char *target_file = nullptr;
        const char *target_object = nullptr;
        unsigned flags = 0;
        ok(H5Lunpack_elink_val(value.data(), value.size(), &flags, &target_file,
                               &target_object),
           "H5Lunpack_elink_val");
        std::filesystem::path path(target_file);
        if (path.is_relative() && !directory.empty())
          path = directory / path;
        block.filename = path.string();
        block.dataset = target_object;
      } else {
        block.filename = master_;
        block.dataset = name;
      }
      std::error_code ignored;
      if (block.filename != master_ &&
          !std::filesystem::is_regular_file(block.filename, ignored))
        throw std::runtime_error(std::string(name) + " links to " +
                                 block.filename + ", which is not there");
      const Handle data(opened(H5Dopen2(file, name, H5P_DEFAULT),
                               ("opening " + std::string(name)).c_str()),
                        H5Dclose);
      const Handle space(opened(H5Dget_space(data.get()), "H5Dget_space"),
                         H5Sclose);
      if (H5Sget_simple_extent_ndims(space.get()) != 3)
        throw std::runtime_error(std::string(name) +
                                 " is not a stack of images");
      hsize_t dims[3] = {0, 0, 0};
      ok(H5Sget_simple_extent_dims(space.get(), dims, nullptr),
         "H5Sget_simple_extent_dims");
      if (dims[0] == 0)
        continue;
      if (!blocks.empty() &&
          (dims[1] != info->height || dims[2] != info->width))
        throw std::runtime_error(std::string(name) +
                                 " holds frames of another size");
      info->height = dims[1];
      info->width = dims[2];
      block.first = next;
      block.last = next + dims[0] - 1;
      block.offset = 0;
      next += dims[0];
      blocks.push_back(block);
    }
    if (blocks.empty())
      throw std::runtime_error("/entry/data/data_000001 and on hold no frames");
    info->images = next;
    return blocks;
  }

  std::vector<Block> whole_dataset(std::uint64_t images) {
    Block block;
    block.filename = master_;
    block.dataset = kImageData;
    block.first = 0;
    block.last = images == 0 ? 0 : images - 1;
    block.offset = 0;
    return {block};
  }

  std::shared_ptr<const std::vector<std::uint64_t>> masked_;
  std::string master_;
  std::vector<Block> blocks_;
  bool dispatched_ = false;
  std::uint64_t first_frame_ = 0;
  std::uint64_t last_frame_ = std::numeric_limits<std::uint64_t>::max();
};

} // namespace

std::unique_ptr<Series> nxmx(std::string master) {
  return std::unique_ptr<Series>(new Nxmx(std::move(master)));
}

std::vector<SourceBlock> nxmx_blocks(const std::string &master, Info *info) {
  Nxmx series(master);
  if (!series.try_open(info))
    throw std::runtime_error(master + " could not be opened as an NXmx series");
  std::vector<SourceBlock> out;
  for (const Block &b : series.blocks())
    out.push_back({b.filename, b.dataset, b.first, b.last, b.offset});
  return out;
}

} // namespace mxi::series
