// The processing's summary: see processing.hh.

#include "processing.hh"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <stdexcept>

#include "expt.hh"
#include "rflx.hh"

namespace mxi {

namespace {
bool g_verbose = false;
} // namespace

void set_verbose(bool on) { g_verbose = on; }
bool verbose() { return g_verbose; }

void detail(const char *format, ...) {
  if (!g_verbose)
    return;
  std::va_list args;
  va_start(args, format);
  std::vprintf(format, args);
  va_end(args);
}

GeometryFigures geometry_figures(const Experiment &e) {
  GeometryFigures g;
  g.wavelength = e.beam.wavelength;
  if (e.detector.size() == 0)
    return g;
  const Panel &p = e.detector[0];
  const Vec3 normal = p.fast.cross(p.slow);
  g.distance = std::fabs(p.origin.dot(normal / normal.norm()));
  if (const auto hit = e.detector.intersect(e.beam.s0())) {
    g.beam_on_panel = true;
    g.beam_x = std::get<1>(*hit);
    g.beam_y = std::get<2>(*hit);
  }
  return g;
}

std::string with_esd(double value, double esd, int plain) {
  char text[64];
  if (!(esd > 0.0) || !std::isfinite(esd)) {
    std::snprintf(text, sizeof text, "%.*f", plain, value);
    return text;
  }
  int places = static_cast<int>(-std::floor(std::log10(esd)));
  if (esd * std::pow(10.0, places) < 2.0)
    ++places;
  places = std::max(0, std::min(places, 8));
  const long digits = std::lround(esd * std::pow(10.0, places));
  std::snprintf(text, sizeof text, "%.*f(%ld)", places, value, digits);
  return text;
}

json::Value to_json(const GeometryFigures &g) {
  json::Object o;
  o["wavelength"] = json::Value(g.wavelength);
  o["distance"] = json::Value(g.distance);
  if (g.beam_on_panel)
    o["beam_centre_px"] =
        json::Value(json::Array{json::Value(g.beam_x), json::Value(g.beam_y)});
  return json::Value(std::move(o));
}

Summary::Summary(const std::string &program, const std::string &input) {
  if (!input.empty()) {
    try {
      if (rflx::is_hdf5(input))
        if (const auto s = rflx::read_summary(input))
          before_ = steps(*s);
    } catch (const std::exception &) {
      // A summary not to be read is not a run not to be made.
    }
  }
  own_["program"] = json::Value(program);
  char stamp[32];
  const std::time_t t = std::time(nullptr);
  std::strftime(stamp, sizeof stamp, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
  own_["time"] = json::Value(std::string(stamp));
}

void Summary::set(const std::string &name, json::Value value) {
  own_[name] = std::move(value);
}

void Summary::table(const std::string &name,
                    const std::vector<std::string> &columns,
                    const std::vector<std::vector<json::Value>> &values) {
  if (values.size() != columns.size())
    throw std::runtime_error("summary table " + name +
                             ": a column a name wanted");
  json::Object t;
  json::Array names;
  for (std::size_t c = 0; c < columns.size(); ++c) {
    names.push_back(json::Value(columns[c]));
    t[columns[c]] =
        json::Value(json::Array(values[c].begin(), values[c].end()));
  }
  t["columns"] = json::Value(std::move(names));
  own_[name] = json::Value(std::move(t));
}

json::Value Summary::document() const {
  json::Array all(before_.begin(), before_.end());
  all.push_back(json::Value(own_));
  json::Object d;
  d["steps"] = json::Value(std::move(all));
  return json::Value(std::move(d));
}

std::vector<json::Value> Summary::steps(const json::Value &summary) {
  std::vector<json::Value> out;
  if (summary.is_object() && summary["steps"].is_array())
    for (const json::Value &s : summary["steps"].as_array())
      out.push_back(s);
  return out;
}

} // namespace mxi
