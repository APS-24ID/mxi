#include "nxmx_import.hh"

#include <hdf5.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <functional>
#include <map>
#include <random>
#include <stdexcept>
#include <utility>

#include "linalg.hh"

namespace mxi {

namespace {

constexpr double kPi = 3.14159265358979323846;

// The link information H5Literate gives its callback changed type in 1.12.
#if H5_VERSION_GE(1, 12, 0)
using LinkInfo = H5L_info2_t;
#else
using LinkInfo = H5L_info_t;
#endif
constexpr double kHcKeVA = 12.398419843320026; // h c, keV A

//: An HDF5 identifier closed with its owner.
class H5 {
public:
  H5(hid_t id, herr_t (*close)(hid_t)) : id_(id), close_(close) {}
  H5(const H5 &) = delete;
  H5 &operator=(const H5 &) = delete;
  ~H5() {
    if (id_ >= 0)
      close_(id_);
  }
  hid_t get() const { return id_; }
  bool ok() const { return id_ >= 0; }

private:
  hid_t id_;
  herr_t (*close_)(hid_t);
};

bool exists(hid_t file, const std::string &path) {
  if (path.empty() || path == "/")
    return true;
  // Each component in turn, as H5Lexists wants its parents to exist.
  std::string so_far;
  std::size_t at = path[0] == '/' ? 1 : 0;
  if (path[0] == '/')
    so_far = "";
  while (at <= path.size()) {
    const std::size_t next = path.find('/', at);
    const std::string part = path.substr(
        at, next == std::string::npos ? std::string::npos : next - at);
    so_far += "/" + part;
    if (H5Lexists(file, so_far.c_str(), H5P_DEFAULT) <= 0)
      return false;
    if (next == std::string::npos)
      break;
    at = next + 1;
  }
  return H5Oexists_by_name(file, path.c_str(), H5P_DEFAULT) > 0;
}

std::vector<double> read_doubles(hid_t file, const std::string &path) {
  H5 d(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose);
  if (!d.ok())
    throw std::runtime_error("cannot open " + path);
  H5 space(H5Dget_space(d.get()), H5Sclose);
  const hssize_t n = H5Sget_simple_extent_npoints(space.get());
  std::vector<double> out(static_cast<std::size_t>(n > 0 ? n : 0));
  if (n > 0 && H5Dread(d.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                       H5P_DEFAULT, out.data()) < 0)
    throw std::runtime_error("cannot read " + path + " as numbers");
  return out;
}

std::string read_string(hid_t file, const std::string &path) {
  H5 d(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose);
  if (!d.ok())
    throw std::runtime_error("cannot open " + path);
  H5 type(H5Dget_type(d.get()), H5Tclose);
  if (H5Tget_class(type.get()) != H5T_STRING)
    throw std::runtime_error(path + " is not a string");
  std::string out;
  if (H5Tis_variable_str(type.get()) > 0) {
    H5 mem(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(mem.get(), H5T_VARIABLE);
    H5Tset_cset(mem.get(), H5Tget_cset(type.get()));
    H5 space(H5Dget_space(d.get()), H5Sclose);
    if (H5Sget_simple_extent_npoints(space.get()) < 1)
      return out;
    std::vector<char *> values(
        static_cast<std::size_t>(H5Sget_simple_extent_npoints(space.get())),
        nullptr);
    if (H5Dread(d.get(), mem.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                values.data()) >= 0 &&
        values[0])
      out = values[0];
    H5Dvlen_reclaim(mem.get(), space.get(), H5P_DEFAULT, values.data());
  } else {
    const std::size_t size = H5Tget_size(type.get());
    std::vector<char> buffer(size + 1, '\0');
    H5 mem(H5Tcopy(type.get()), H5Tclose);
    if (H5Dread(d.get(), mem.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                buffer.data()) >= 0)
      out.assign(buffer.data(), strnlen(buffer.data(), size));
  }
  return out;
}

std::optional<std::string> attr_string(hid_t object, const char *name) {
  if (H5Aexists(object, name) <= 0)
    return std::nullopt;
  H5 a(H5Aopen(object, name, H5P_DEFAULT), H5Aclose);
  H5 type(H5Aget_type(a.get()), H5Tclose);
  if (H5Tget_class(type.get()) != H5T_STRING)
    return std::nullopt;
  std::string out;
  if (H5Tis_variable_str(type.get()) > 0) {
    H5 mem(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(mem.get(), H5T_VARIABLE);
    H5Tset_cset(mem.get(), H5Tget_cset(type.get()));
    char *value = nullptr;
    if (H5Aread(a.get(), mem.get(), &value) >= 0 && value) {
      out = value;
      H5free_memory(value);
    }
  } else {
    const std::size_t size = H5Tget_size(type.get());
    std::vector<char> buffer(size + 1, '\0');
    if (H5Aread(a.get(), type.get(), buffer.data()) >= 0)
      out.assign(buffer.data(), strnlen(buffer.data(), size));
  }
  return out;
}

std::optional<std::vector<double>> attr_doubles(hid_t object,
                                                const char *name) {
  if (H5Aexists(object, name) <= 0)
    return std::nullopt;
  H5 a(H5Aopen(object, name, H5P_DEFAULT), H5Aclose);
  H5 space(H5Aget_space(a.get()), H5Sclose);
  const hssize_t n = H5Sget_simple_extent_npoints(space.get());
  std::vector<double> out(static_cast<std::size_t>(n > 0 ? n : 0));
  if (n > 0 && H5Aread(a.get(), H5T_NATIVE_DOUBLE, out.data()) < 0)
    return std::nullopt;
  return out;
}

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

//: mm per unit of length; a unit not known is a failure, a missing one mm.
double to_mm(const std::optional<std::string> &units, const std::string &what,
             std::vector<std::string> *notes) {
  if (!units || units->empty()) {
    notes->push_back(what + " has no units: taken as mm");
    return 1.0;
  }
  const std::string u = lower(*units);
  if (u == "m" || u == "metre" || u == "meter" || u == "metres" ||
      u == "meters")
    return 1000.0;
  if (u == "mm" || u == "millimetre" || u == "millimeter" ||
      u == "millimetres" || u == "millimeters")
    return 1.0;
  if (u == "cm")
    return 10.0;
  if (u == "um" || u == "micron" || u == "microns" || u == "micrometre" ||
      u == "micrometer")
    return 1e-3;
  if (u == "nm")
    return 1e-6;
  throw std::runtime_error(what + " has units '" + *units +
                           "', which are not a length known here");
}

//: degrees per unit of angle.
double to_degrees(const std::optional<std::string> &units,
                  const std::string &what, std::vector<std::string> *notes) {
  if (!units || units->empty()) {
    notes->push_back(what + " has no units: taken as degrees");
    return 1.0;
  }
  const std::string u = lower(*units);
  if (u == "deg" || u == "degree" || u == "degrees")
    return 1.0;
  if (u == "rad" || u == "radian" || u == "radians")
    return 180.0 / kPi;
  throw std::runtime_error(what + " has units '" + *units +
                           "', which are not an angle known here");
}

//: One step of a NeXus transformation chain, in McStas mm and degrees.
struct Step {
  std::string path, name, type;
  Vec3 vector, offset;
  Vec3
      raw; //: the vector as the file gives it, which need not be of unit length
  std::vector<double> values; //: mm or degrees, one a frame where it moves
  std::string depends_on;
};

std::string parent_of(const std::string &path) {
  const std::size_t slash = path.rfind('/');
  return slash == 0 || slash == std::string::npos ? "/" : path.substr(0, slash);
}

std::string resolve(const std::string &from_group, const std::string &target) {
  if (target.empty() || target[0] == '/')
    return target;
  std::filesystem::path p = std::filesystem::path(from_group) / target;
  return p.lexically_normal().generic_string();
}

Step read_step(hid_t file, const std::string &path,
               std::vector<std::string> *notes) {
  H5 d(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose);
  if (!d.ok())
    throw std::runtime_error("the transformation " + path + " does not open");
  Step s;
  s.path = path;
  s.name = path.substr(path.rfind('/') + 1);
  s.type = lower(attr_string(d.get(), "transformation_type").value_or(""));
  const auto vec = attr_doubles(d.get(), "vector");
  if (!vec || vec->size() != 3)
    throw std::runtime_error(path + " has no vector");
  s.vector = Vec3{(*vec)[0], (*vec)[1], (*vec)[2]};
  s.raw = s.vector;
  if (s.vector.norm() > 0.0)
    s.vector = s.vector / s.vector.norm();
  const auto off = attr_doubles(d.get(), "offset");
  if (off && off->size() == 3) {
    // An offset without offset_units is in the transformation's own units, as
    // nxmx reads it: Diamond's Eiger masters write offsets so, in metres, and
    // taking them as millimetres put the detector a thousandth of the way out
    // from the beam. A rotation's units are an angle, which no offset can be
    // in, so there the length is still taken as millimetres, with a note.
    auto offset_units = attr_string(d.get(), "offset_units");
    if ((!offset_units || offset_units->empty()) && s.type == "translation") {
      const auto own = attr_string(d.get(), "units");
      if (own && !own->empty()) {
        notes->push_back(path + "'s offset has no offset_units: taken in its " +
                         "units, " + *own + ", as nxmx does");
        offset_units = own;
      }
    }
    const double k = to_mm(offset_units, path + "'s offset", notes);
    s.offset = Vec3{(*off)[0] * k, (*off)[1] * k, (*off)[2] * k};
  }
  s.depends_on = attr_string(d.get(), "depends_on").value_or(".");
  const std::vector<double> raw = read_doubles(file, path);
  const auto units = attr_string(d.get(), "units");
  double k = 1.0;
  if (s.type == "translation")
    k = to_mm(units, path, notes);
  else if (s.type == "rotation")
    k = to_degrees(units, path, notes);
  else
    throw std::runtime_error(path + " is neither a translation nor a rotation");
  for (double v : raw)
    s.values.push_back(v * k);
  if (s.values.empty())
    s.values.push_back(0.0);
  return s;
}

//: The chain from `start`, innermost first, to ".".
std::vector<Step> chain(hid_t file, const std::string &start,
                        std::vector<std::string> *notes) {
  std::vector<Step> out;
  std::string at = start;
  while (!at.empty() && at != ".") {
    if (out.size() > 64)
      throw std::runtime_error("the transformation chain from " + start +
                               " does not end");
    if (!exists(file, at))
      throw std::runtime_error("the transformation " + at + " does not exist");
    out.push_back(read_step(file, at, notes));
    const std::string next = out.back().depends_on;
    at = next == "." ? next : resolve(parent_of(at), next);
  }
  return out;
}

Mat3 rotation_deg(const Vec3 &axis, double degrees) {
  return rotation(axis, degrees * kPi / 180.0);
}

//: A point and the rotation applied to it by a chain, at each step's first
//: value: p -> R p + t, outermost last.
struct Pose {
  Mat3 R = Mat3::identity();
  Vec3 t;
};

Pose pose_of(const std::vector<Step> &steps) {
  Pose p;
  for (const Step &s : steps) { // innermost first: each wraps what came before
    Mat3 r = Mat3::identity();
    Vec3 shift = s.offset;
    if (s.type == "rotation")
      r = rotation_deg(s.vector, s.values[0]);
    else
      shift = shift + s.vector * s.values[0];
    p.R = r * p.R;
    p.t = r * p.t + shift;
  }
  return p;
}

//: McStas, which NeXus uses, to imgCIF, which DIALS uses: a half turn about y.
Vec3 to_imgcif(const Vec3 &v) { return Vec3{-v.x, v.y, -v.z}; }

//: A vector as JSON, -0 written as 0: the frame's half turn negates zeros.
json::Value vec(const Vec3 &v) {
  return json::Array{v.x + 0.0, v.y + 0.0, v.z + 0.0};
}

std::string material_symbol(const std::string &raw) {
  const std::string m = lower(raw);
  if (m == "si" || m == "silicon")
    return "Si";
  if (m == "cdte" || m == "cadmium telluride")
    return "CdTe";
  if (m == "gaas" || m == "gallium arsenide")
    return "GaAs";
  if (m == "ge" || m == "germanium")
    return "Ge";
  return raw;
}

std::string uuid4() {
  std::random_device rd;
  std::mt19937_64 gen((static_cast<std::uint64_t>(rd()) << 32) ^ rd());
  std::uniform_int_distribution<int> hex(0, 15);
  const char *digits = "0123456789abcdef";
  std::string s = "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx";
  for (char &c : s) {
    if (c == 'x')
      c = digits[hex(gen)];
    else if (c == 'y')
      c = digits[8 + hex(gen) % 4];
  }
  return s;
}

std::string now_utc() {
  const std::time_t t = std::time(nullptr);
  char buffer[32];
  std::strftime(buffer, sizeof buffer, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
  return buffer;
}

//: Every group under `root` whose NX_class is `cls`, by path.
std::vector<std::string> groups_of_class(hid_t file, const std::string &root,
                                         const std::string &cls) {
  std::vector<std::string> out;
  struct Ctx {
    hid_t file;
    std::string root, cls;
    std::vector<std::string> *out;
  } ctx{file, root, cls, &out};
  H5 g(H5Gopen2(file, root.c_str(), H5P_DEFAULT), H5Gclose);
  if (!g.ok())
    return out;
  H5Literate(
      g.get(), H5_INDEX_NAME, H5_ITER_INC, nullptr,
      [](hid_t group, const char *name, const LinkInfo *,
         void *data) -> herr_t {
        auto *c = static_cast<Ctx *>(data);
        // A group if it opens as one: the same in every HDF5 version, where
        // asking an object's type is not.
        H5 sub(H5Gopen2(group, name, H5P_DEFAULT), H5Gclose);
        if (sub.ok() &&
            attr_string(sub.get(), "NX_class").value_or("") == c->cls)
          c->out->push_back(c->root + "/" + name);
        return 0;
      },
      &ctx);
  return out;
}

//: A path that is named but does not open: every link to it there, the object
//: not -- an external link into a file that is not beside the master, as a
//: DECTRIS master's saturation_value is into its _meta.h5. Said so, where it
//: leads, rather than taken for absent.
std::optional<std::string> unreachable(hid_t file, const std::string &path) {
  if (path.empty() || path == "/")
    return std::nullopt;
  std::string so_far;
  std::size_t at = path[0] == '/' ? 1 : 0;
  while (at <= path.size()) {
    const std::size_t next = path.find('/', at);
    so_far +=
        "/" + path.substr(at, next == std::string::npos ? std::string::npos
                                                        : next - at);
    if (H5Lexists(file, so_far.c_str(), H5P_DEFAULT) <= 0)
      return std::nullopt; // not named at all
    if (next == std::string::npos)
      break;
    at = next + 1;
  }
  if (H5Oexists_by_name(file, path.c_str(), H5P_DEFAULT) > 0)
    return std::nullopt; // there
#if H5_VERSION_GE(1, 12, 0)
  H5L_info2_t info;
  const herr_t got = H5Lget_info2(file, path.c_str(), &info, H5P_DEFAULT);
#else
  H5L_info_t info;
  const herr_t got = H5Lget_info(file, path.c_str(), &info, H5P_DEFAULT);
#endif
  if (got >= 0 && info.type == H5L_TYPE_EXTERNAL && info.u.val_size > 0) {
    std::vector<char> value(info.u.val_size);
    const char *target_file = nullptr, *target = nullptr;
    unsigned flags = 0;
    if (H5Lget_val(file, path.c_str(), value.data(), value.size(),
                   H5P_DEFAULT) >= 0 &&
        H5Lunpack_elink_val(value.data(), value.size(), &flags, &target_file,
                            &target) >= 0 &&
        target_file && target)
      return path + " links to " + target + " in " + target_file +
             ", which cannot be opened";
  }
  return path + " is named but cannot be opened";
}

std::optional<std::string>
first_existing(hid_t file, std::initializer_list<std::string> paths) {
  for (const std::string &p : paths)
    if (exists(file, p))
      return p;
  return std::nullopt;
}

} // namespace

double attenuation_coefficient(const std::string &material, double wavelength) {
  // Hubbell and Seltzer, NIST mass attenuation coefficients with coherent
  // scattering, energy keV and mu/rho cm^2/g, as cctbx's eltbx tabulates them
  // and dials.import reads them -- the materials dxtbx knows: Si, CdTe, GaAs.
  // Silicon from 2 keV, above its K edge; CdTe and GaAs whole, their edges in
  // the MX range given twice, below and above.
  struct Table {
    const char *name;
    double density; // g/cm^3
    std::vector<std::pair<double, double>> points;
  };
  static const Table silicon{"Si",
                             2.33,
                             {{2.0, 2777.0},
                              {3.0, 978.4},
                              {4.0, 452.9},
                              {5.0, 245.0},
                              {6.0, 147.0},
                              {8.0, 64.68},
                              {10.0, 33.89},
                              {15.0, 10.34},
                              {20.0, 4.464},
                              {30.0, 1.436},
                              {40.0, 0.7012},
                              {50.0, 0.4385},
                              {60.0, 0.3207},
                              {80.0, 0.2228},
                              {100.0, 0.1835}}};
  // cctbx's eltbx table, entry 93: 59 points.
  static const Table cadmium_telluride{
      "CdTe",
      6.2,
      {{1, 7927},        {1.003, 7875},    {1.006, 7824},   {1.006, 8014},
       {1.5, 3291},      {2, 1664},        {3, 614.6},      {3.537, 406.4},
       {3.537, 778.7},   {3.631, 730},     {3.727, 684},    {3.727, 860.1},
       {4, 723},         {4.018, 715.1},   {4.018, 793.4},  {4.177, 722.1},
       {4.341, 656.2},   {4.341, 932.8},   {4.475, 873.9},  {4.612, 813.5},
       {4.612, 943.8},   {4.773, 870.2},   {4.939, 799.9},  {4.939, 865.3},
       {5, 839.2},       {6, 528.6},       {8, 249.2},      {10, 138.1},
       {15, 46.57},      {20, 21.44},      {26.711, 9.834}, {26.711, 29.43},
       {30, 21.82},      {31.814, 18.73},  {31.814, 34.92}, {40, 19.3},
       {50, 10.67},      {60, 6.542},      {80, 3.019},     {100, 1.671},
       {150, 0.6071},    {200, 0.3246},    {300, 0.1628},   {400, 0.1147},
       {500, 0.09291},   {600, 0.08042},   {800, 0.066},    {1000, 0.05742},
       {1250, 0.05043},  {1500, 0.04591},  {2000, 0.0407},  {3000, 0.03649},
       {4000, 0.03525},  {5000, 0.03513},  {6000, 0.03548}, {8000, 0.03687},
       {10000, 0.03857}, {15000, 0.04273}, {20000, 0.04616}}};
  // cctbx's eltbx table, entry 94: 58 points.
  static const Table gallium_arsenide{
      "GaAs",
      5.32,
      {{1, 1917},        {1.05613, 1685},  {1.1154, 1481},   {1.1154, 2772},
       {1.12877, 3180},  {1.1423, 3532},   {1.1423, 4372},   {1.21752, 4130},
       {1.2977, 3657},   {1.2977, 4066},   {1.31034, 3971},  {1.3231, 3879},
       {1.3231, 5652},   {1.34073, 5525},  {1.3586, 5415},   {1.3586, 6266},
       {1.5, 5159},      {1.5265, 4939},   {1.5265, 5278},   {2, 2731},
       {3, 970.2},       {4, 453.9},       {5, 249.5},       {6, 152.4},
       {8, 69.6},        {10, 37.8},       {10.3671, 34.25}, {10.3671, 126},
       {11.0916, 105.9}, {11.8667, 88.99}, {11.8667, 168.5}, {15, 92.2},
       {20, 42.58},      {30, 13.97},      {40, 6.262},      {50, 3.365},
       {60, 2.042},      {80, 0.9587},     {100, 0.5598},    {150, 0.2509},
       {200, 0.1671},    {300, 0.1137},    {400, 0.09371},   {500, 0.08248},
       {600, 0.07484},   {800, 0.06452},   {1000, 0.05751},  {1250, 0.05122},
       {1500, 0.04676},  {2000, 0.04102},  {3000, 0.03538},  {4000, 0.03288},
       {5000, 0.03172},  {6000, 0.03121},  {8000, 0.03117},  {10000, 0.0317},
       {15000, 0.03354}, {20000, 0.03543}}};
  const Table *t = material == "Si"     ? &silicon
                   : material == "CdTe" ? &cadmium_telluride
                   : material == "GaAs" ? &gallium_arsenide
                                        : nullptr;
  if (!t || !(wavelength > 0.0))
    return 0.0;
  const double energy = kHcKeVA / wavelength;
  // cctbx's find_energy_index: the first point above the energy, and the
  // interval before it -- at an edge, given twice, the interval above it; log
  // log between its ends, as mu_rho_at_ev. Outside the table, nothing.
  std::size_t above = 0;
  while (above < t->points.size() && !(energy < t->points[above].first))
    ++above;
  if (above == 0 || above == t->points.size())
    return 0.0;
  const auto [e0, m0] = t->points[above - 1];
  const auto [e1, m1] = t->points[above];
  const double mu_rho =
      std::exp(std::log(m0) + (std::log(m1) - std::log(m0)) *
                                  (std::log(energy) - std::log(e0)) /
                                  (std::log(e1) - std::log(e0)));
  return mu_rho * t->density / 10.0; // 1/cm to 1/mm
}

json::Value import_nxmx(const std::string &master, const ImportOverrides &o,
                        std::vector<std::string> *notes) {
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  H5 file(H5Fopen(master.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
  if (!file.ok())
    throw std::runtime_error("cannot open " + master + " as HDF5");
  const hid_t f = file.get();

  // The beam.
  double wavelength = 0.0;
  if (o.wavelength) {
    wavelength = *o.wavelength;
    notes->push_back("the wavelength given, " + std::to_string(wavelength) +
                     " A");
  } else if (auto p = first_existing(
                 f, {"/entry/instrument/beam/incident_wavelength",
                     "/entry/sample/beam/incident_wavelength"})) {
    const std::vector<double> w = read_doubles(f, *p);
    if (w.empty())
      throw std::runtime_error(*p + " is empty");
    H5 d(H5Dopen2(f, p->c_str(), H5P_DEFAULT), H5Dclose);
    const std::string u =
        lower(attr_string(d.get(), "units").value_or("angstrom"));
    const double k = (u == "nm") ? 10.0 : (u == "m" ? 1e10 : 1.0);
    wavelength = w[0] * k;
    if (w.size() > 1 &&
        std::fabs(w.back() - w.front()) > 1e-9 * std::fabs(w.front()))
      notes->push_back(*p + " varies; the first value used");
  } else if (auto e =
                 first_existing(f, {"/entry/instrument/beam/incident_energy",
                                    "/entry/sample/beam/incident_energy"})) {
    const std::vector<double> energy = read_doubles(f, *e);
    H5 d(H5Dopen2(f, e->c_str(), H5P_DEFAULT), H5Dclose);
    const std::string u = lower(attr_string(d.get(), "units").value_or("eV"));
    const double kev = energy.at(0) * (u == "kev" ? 1.0 : 1e-3);
    wavelength = kHcKeVA / kev;
    notes->push_back("the wavelength from incident_energy");
  } else {
    throw std::runtime_error(
        "no incident_wavelength or incident_energy: give --wavelength");
  }
  json::Object beam{{"__id__", "monochromatic"},
                    {"direction", json::Array{0.0, 0.0, 1.0}},
                    {"wavelength", wavelength},
                    {"divergence", 0.0},
                    {"sigma_divergence", 0.0},
                    {"polarization_normal", json::Array{0.0, 1.0, 0.0}},
                    {"polarization_fraction", 0.999},
                    {"flux", 0.0},
                    {"transmission", 1.0},
                    {"probe", "x-ray"},
                    {"sample_to_source_distance", 0.0}};

  // The detector: a panel per module, its origin and axes through the
  // transformation chains, then in DIALS's frame.
  const std::string detector = "/entry/instrument/detector";
  if (!exists(f, detector))
    throw std::runtime_error("no " + detector);
  std::vector<std::string> modules =
      groups_of_class(f, detector, "NXdetector_module");
  if (modules.empty())
    throw std::runtime_error("no NXdetector_module under " + detector);
  std::string material = "Si";
  if (exists(f, detector + "/sensor_material"))
    material = material_symbol(read_string(f, detector + "/sensor_material"));
  double thickness = 0.0;
  if (exists(f, detector + "/sensor_thickness")) {
    H5 d(H5Dopen2(f, (detector + "/sensor_thickness").c_str(), H5P_DEFAULT),
         H5Dclose);
    thickness = read_doubles(f, detector + "/sensor_thickness").at(0) *
                to_mm(attr_string(d.get(), "units"), "sensor_thickness", notes);
  }
  double mu = o.mu ? *o.mu : attenuation_coefficient(material, wavelength);
  if (!o.mu && mu == 0.0)
    notes->push_back(
        "no attenuation coefficient tabulated for " + material +
        " at this wavelength: mu 0, so no parallax correction; give --mu");
  // The trusted range's top: the lower of two limits. The detector's own --
  // the count a photon counter can still correct for, set by the exposure
  // time, as low as 31881 at 2 ms and higher for longer -- which the master
  // gives as saturation_value or the count cutoff. And the data type's: the
  // two largest values of the image's type are markers, a bad pixel and a
  // tile join, so the largest count is 2^bits - 3. The bit depth from the
  // master, or the first data file's own type; not the virtual dataset's,
  // which a writer may widen (int64 over 32-bit counts).
  double trusted_max = 0.0;
  if (o.trusted_max) {
    trusted_max = *o.trusted_max;
  } else {
    std::optional<double> detector_limit;
    std::optional<int> bits;
    std::string bits_from;
    const std::vector<std::string> limits = {
        detector + "/saturation_value",
        detector + "/detectorSpecific/countrate_correction_count_cutoff"};
    for (const std::string &p : limits) {
      if (exists(f, p)) {
        detector_limit = read_doubles(f, p).at(0);
        break;
      }
      if (auto why = unreachable(f, p))
        notes->push_back(*why);
    }
    for (const char *name : {"/bit_depth_image", "/bit_depth_readout"}) {
      const std::string p = detector + name;
      if (exists(f, p)) {
        bits = static_cast<int>(read_doubles(f, p).at(0));
        bits_from = p;
        break;
      }
      if (auto why = unreachable(f, p))
        notes->push_back(*why);
    }
    if (!bits) {
      // The data files' own type: the first linked one there.
      for (int k = 1; k <= 9 && !bits; ++k) {
        char name[32];
        std::snprintf(name, sizeof name, "/entry/data/data_%06d", k);
        if (!exists(f, name))
          continue;
        H5 d(H5Dopen2(f, name, H5P_DEFAULT), H5Dclose);
        if (!d.ok())
          continue;
        H5 t(H5Dget_type(d.get()), H5Tclose);
        if (t.ok() && H5Tget_class(t.get()) == H5T_INTEGER) {
          bits = static_cast<int>(8 * H5Tget_size(t.get()));
          bits_from = std::string(name) + "'s type";
        }
      }
    }
    const std::optional<double> type_limit =
        bits ? std::optional<double>(std::ldexp(1.0, *bits) - 3.0)
             : std::nullopt;
    if (detector_limit && type_limit) {
      trusted_max = std::min(*detector_limit, *type_limit);
      if (*type_limit < *detector_limit)
        notes->push_back(
            "the detector's count limit, " +
            std::to_string(static_cast<long long>(*detector_limit)) +
            ", is above what " + std::to_string(*bits) +
            "-bit data can hold below its markers: the trusted "
            "range's top is " +
            std::to_string(static_cast<long long>(*type_limit)));
    } else if (detector_limit) {
      trusted_max = *detector_limit;
    } else if (type_limit) {
      trusted_max = *type_limit;
      notes->push_back("no count limit could be read: the trusted range's top "
                       "is " +
                       std::to_string(static_cast<long long>(*type_limit)) +
                       ", the largest " + std::to_string(*bits) +
                       "-bit count below the markers (" + bits_from +
                       "); give --trusted-max for the detector's own");
    } else {
      trusted_max = 2147483647.0;
      notes->push_back(
          "neither the detector's count limit nor the data's bit depth could "
          "be read: the trusted range's top is 2147483647, as dxtbx takes a "
          "file without one, so no count is distrusted for its size -- the "
          "markers are recognised either way. Give --trusted-max, or put the "
          "files the master links to beside it");
    }
  }

  json::Array panels;
  for (const std::string &module : modules) {
    const std::vector<Step> fast_chain =
        chain(f, module + "/fast_pixel_direction", notes);
    const std::vector<Step> slow_chain =
        chain(f, module + "/slow_pixel_direction", notes);
    // The fast and slow directions' own steps give the pixel size and the
    // directions; everything they depend on places the module.
    const Step &fast = fast_chain.front();
    const Step &slow = slow_chain.front();
    const std::vector<Step> rest(fast_chain.begin() + 1, fast_chain.end());
    const Pose place = pose_of(rest);
    Vec3 origin = place.t + place.R * fast.offset;
    Vec3 fast_axis = place.R * fast.vector;
    const Pose slow_place =
        pose_of(std::vector<Step>(slow_chain.begin() + 1, slow_chain.end()));
    Vec3 slow_axis = slow_place.R * slow.vector;
    origin = to_imgcif(origin);
    fast_axis = to_imgcif(fast_axis);
    slow_axis = to_imgcif(slow_axis);
    const std::vector<double> size = read_doubles(f, module + "/data_size");
    const std::vector<double> start =
        exists(f, module + "/data_origin")
            ? read_doubles(f, module + "/data_origin")
            : std::vector<double>{0.0, 0.0};
    if (size.size() < 2)
      throw std::runtime_error(module + "/data_size is not two numbers");
    // NeXus gives slow then fast; DIALS fast then slow.
    const long long nx = static_cast<long long>(size[size.size() - 1]);
    const long long ny = static_cast<long long>(size[size.size() - 2]);
    const long long ox = static_cast<long long>(start[start.size() - 1]);
    const long long oy = static_cast<long long>(start[start.size() - 2]);
    panels.push_back(json::Object{
        {"name", module},
        {"type", "SENSOR_PAD"},
        {"fast_axis", vec(fast_axis)},
        {"slow_axis", vec(slow_axis)},
        {"origin", vec(origin)},
        {"raw_image_offset", json::Array{ox, oy}},
        {"image_size", json::Array{nx, ny}},
        {"pixel_size", json::Array{fast.values[0], slow.values[0]}},
        {"trusted_range", json::Array{0.0, trusted_max}},
        {"thickness", thickness},
        {"material", material},
        {"mu", mu},
        {"identifier", ""},
        {"mask", json::Array{}},
        {"gain", 1.0},
        {"pedestal", 0.0},
        {"px_mm_strategy",
         json::Object{{"type", mu > 0.0 && thickness > 0.0
                                   ? "ParallaxCorrectedPxMmStrategy"
                                   : "SimplePxMmStrategy"}}}});
  }
  if (panels.size() > 1)
    notes->push_back(std::to_string(panels.size()) +
                     " modules, so as many panels: a detector of several "
                     "panels is untested here");

  // The overrides of where the detector is, on the first panel's plane.
  auto &p0 = panels[0].as_object();
  const auto get = [](const json::Value &v) {
    const auto &a = v.as_array();
    return Vec3{a[0].as_number(), a[1].as_number(), a[2].as_number()};
  };
  if (o.distance || o.beam_centre) {
    const Vec3 fast_axis = get(p0["fast_axis"]),
               slow_axis = get(p0["slow_axis"]);
    Vec3 normal = fast_axis.cross(slow_axis);
    normal = normal / normal.norm();
    Vec3 origin = get(p0["origin"]);
    if (origin.dot(normal) < 0)
      normal = normal * -1.0;
    Vec3 shift;
    if (o.distance) {
      const double now = origin.dot(normal);
      shift = normal * (*o.distance - now);
      notes->push_back("the distance given, " + std::to_string(*o.distance) +
                       " mm");
    }
    Vec3 moved = origin + shift;
    if (o.beam_centre) {
      // Where the beam meets the panel's plane, and the origin put so that
      // that point is the pixel given.
      const Vec3 s0 = Vec3{0.0, 0.0, -1.0}; // the beam, towards the detector
      const double t = moved.dot(normal) / s0.dot(normal);
      const Vec3 hit = s0 * t;
      const auto px = p0["pixel_size"].as_array();
      moved = hit - fast_axis * ((*o.beam_centre)[0] * px[0].as_number()) -
              slow_axis * ((*o.beam_centre)[1] * px[1].as_number());
      notes->push_back("the beam centre given, pixel " +
                       std::to_string((*o.beam_centre)[0]) + ", " +
                       std::to_string((*o.beam_centre)[1]));
    }
    const Vec3 delta = moved - origin;
    for (json::Value &panel : panels)
      panel.as_object()["origin"] =
          vec(get(panel.as_object()["origin"]) + delta);
  }
  json::Object hierarchy{
      {"name", ""},
      {"type", ""},
      {"fast_axis", json::Array{1.0, 0.0, 0.0}},
      {"slow_axis", json::Array{0.0, 1.0, 0.0}},
      {"origin", json::Array{0.0, 0.0, 0.0}},
      {"raw_image_offset", json::Array{0, 0}},
      {"image_size", json::Array{0, 0}},
      {"pixel_size", json::Array{0.0, 0.0}},
      {"trusted_range", json::Array{0.0, 0.0}},
      {"thickness", 0.0},
      {"material", ""},
      {"mu", 0.0},
      {"identifier", ""},
      {"mask", json::Array{}},
      {"gain", 1.0},
      {"pedestal", 0.0},
      {"px_mm_strategy", json::Object{{"type", "SimplePxMmStrategy"}}}};
  json::Array children;
  for (std::size_t k = 0; k < panels.size(); ++k)
    children.push_back(json::Object{{"panel", static_cast<long long>(k)}});
  hierarchy["children"] = children;

  // The goniometer: the rotations of the sample's chain, innermost first; the
  // scan axis the one that moves.
  std::string sample_start;
  if (exists(f, "/entry/sample/depends_on"))
    sample_start =
        resolve("/entry/sample", read_string(f, "/entry/sample/depends_on"));
  if (sample_start.empty())
    throw std::runtime_error("/entry/sample has no depends_on: no goniometer");
  json::Array axes, angles, names;
  int scan_axis = -1;
  std::vector<double> scan_values;
  for (const Step &s : chain(f, sample_start, notes)) {
    if (s.type != "rotation")
      continue;
    const bool moves = s.values.size() > 1 &&
                       std::fabs(s.values.back() - s.values.front()) > 1e-12;
    if (moves) {
      if (scan_axis >= 0)
        throw std::runtime_error("more than one goniometer axis moves: " +
                                 s.name + " and the scan axis before it");
      scan_axis = static_cast<int>(axes.size());
      scan_values = s.values;
    }
    // As the file gives it, as dials.import writes it: normalised where used.
    axes.push_back(vec(to_imgcif(s.raw)));
    angles.push_back(moves ? 0.0 : s.values[0]);
    names.push_back(s.name);
  }
  if (scan_axis < 0)
    throw std::runtime_error("no goniometer axis moves: not a rotation scan");
  json::Object goniometer{{"axes", axes},
                          {"angles", angles},
                          {"names", names},
                          {"scan_axis", scan_axis}};

  // The scan.
  const long long frames = static_cast<long long>(scan_values.size());
  long long first = 1, last = frames;
  if (o.image_range) {
    first = (*o.image_range)[0];
    last = (*o.image_range)[1];
    if (first < 1 || last > frames || first > last)
      throw std::runtime_error("--image-range outside the " +
                               std::to_string(frames) + " images");
  }
  double exposure = 0.0;
  if (exists(f, detector + "/count_time"))
    exposure = read_doubles(f, detector + "/count_time").at(0);
  json::Array oscillation, epochs, exposures, indices;
  for (long long i = first - 1; i < last; ++i) {
    oscillation.push_back(scan_values[static_cast<std::size_t>(i)]);
    epochs.push_back(0.0);
    exposures.push_back(exposure);
    indices.push_back(i);
  }
  json::Object scan{{"image_range", json::Array{first, last}},
                    {"batch_offset", 0},
                    {"valid_image_ranges", json::Object{}},
                    {"properties", json::Object{{"epochs", epochs},
                                                {"exposure_time", exposures},
                                                {"oscillation", oscillation}}}};

  const std::string template_path =
      std::filesystem::absolute(master).lexically_normal().string();
  json::Object imageset{{"__id__", "ImageSequence"},
                        {"template", template_path},
                        {"single_file_indices", indices},
                        {"mask", json::Value()},
                        {"gain", json::Value()},
                        {"pedestal", json::Value()},
                        {"dx", json::Value()},
                        {"dy", json::Value()},
                        {"params", json::Object{{"dynamic_shadowing", "Auto"},
                                                {"multi_panel", false}}}};
  json::Object experiment{
      {"__id__", "Experiment"}, {"identifier", uuid4()}, {"beam", 0},
      {"detector", 0},          {"goniometer", 0},       {"scan", 0},
      {"imageset", 0}};
  return json::Object{
      {"__id__", "ExperimentList"},
      {"experiment", json::Array{experiment}},
      {"history", json::Array{now_utc() + "|mxi_import|0.1"}},
      {"imageset", json::Array{imageset}},
      {"beam", json::Array{beam}},
      {"detector",
       json::Array{json::Object{{"panels", panels}, {"hierarchy", hierarchy}}}},
      {"goniometer", json::Array{goniometer}},
      {"scan", json::Array{scan}},
      {"crystal", json::Array{}},
      {"profile", json::Array{}},
      {"scaling_model", json::Array{}}};
}

} // namespace mxi
