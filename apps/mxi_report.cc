// mxi_report: an HTML report of how a file's data were processed and what came
// of it (docs/output.md) -- each step the processing's summary records, its
// figures typeset and its tables drawn, after dials.scale's report. Plotly is
// loaded from its CDN, as dials.scale's and mxeq's reports load it: the report
// is one file that can be mailed, and draws wherever there is a network.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "args.hh"
#include "expt.hh"
#include "processing.hh"
#include "refl.hh"
#include "rflx.hh"
#include "symmetry.hh"

namespace mxi {

namespace {

const char *kPlotly = "https://cdn.plot.ly/plotly-2.27.0.min.js";

void usage(const char *program) {
  std::printf(
      "usage: %s FILE.rflx [-o REPORT.html]\n\n"
      "An HTML report of the processing FILE records (docs/output.md): each\n"
      "step's figures and tables, and plots of them -- spots per image,\n"
      "refinement's RMSDs and the cell over the scan, integration and the\n"
      "merging statistics by resolution, the scale and <I/sigma> by image.\n"
      "  -o REPORT.html   where to write it (mxi_report.html)\n",
      program);
}

std::string escape(const std::string &s) {
  std::string out;
  for (char c : s) {
    switch (c) {
    case '&':
      out += "&amp;";
      break;
    case '<':
      out += "&lt;";
      break;
    case '>':
      out += "&gt;";
      break;
    case '"':
      out += "&quot;";
      break;
    default:
      out += c;
    }
  }
  return out;
}

std::string number(double v, int places = -1) {
  char text[64];
  if (!std::isfinite(v))
    return "&ndash;";
  if (places >= 0)
    std::snprintf(text, sizeof text, "%.*f", places, v);
  else if (v == std::floor(v) && std::fabs(v) < 1e12)
    std::snprintf(text, sizeof text, "%.0f", v);
  else
    std::snprintf(text, sizeof text, "%.4g", v);
  return text;
}

std::string text_of(const json::Value &v) {
  if (v.is_number())
    return number(v.as_number());
  if (v.is_string())
    return escape(v.as_string());
  if (v.is_array()) {
    std::string out;
    for (const json::Value &x : v.as_array())
      out += (out.empty() ? "" : " ") + text_of(x);
    return out;
  }
  if (v.is_object()) {
    std::string out;
    for (const auto &[k, x] : v.as_object())
      out += (out.empty() ? "" : "; ") + escape(k) + " " + text_of(x);
    return out;
  }
  return "";
}

bool is_table(const json::Value &v) {
  return v.is_object() && v["columns"].is_array();
}

// A table's column as numbers.
json::Array column(const json::Value &table, const std::string &name) {
  if (table.is_object() && table[name].is_array())
    return table[name].as_array();
  return {};
}

std::string pretty(const std::string &name) {
  static const std::map<std::string, std::string> names = {
      {"d_max", "d<sub>max</sub> (&Aring;)"},
      {"d_min", "d<sub>min</sub> (&Aring;)"},
      {"cc_half", "CC<sub>&frac12;</sub>"},
      {"cc_anom", "CC<sub>anom</sub>"},
      {"i_over_sigma", "&lang;I/&sigma;&rang;"},
      {"isigma_sum", "&lang;I/&sigma;&rang; summed"},
      {"isigma_prf", "&lang;I/&sigma;&rang; fitted"},
      {"rmerge", "R<sub>merge</sub>"},
      {"rmeas", "R<sub>meas</sub>"},
      {"rpim", "R<sub>pim</sub>"},
      {"mean_i", "&lang;I&rang;"},
      {"cc_prf", "CC profile"},
      {"ibg", "background"},
      {"rmsd_xy_px", "RMSD xy (px)"},
      {"rmsd_x_px", "RMSD x (px)"},
      {"rmsd_y_px", "RMSD y (px)"},
      {"rmsd_z_images", "RMSD z (images)"},
      {"sigma_b_deg", "&sigma;<sub>b</sub> (&deg;)"},
      {"sigma_m_deg", "&sigma;<sub>m</sub> (&deg;)"},
      {"z_cc", "Z-CC"},
      {"cc", "CC"},
  };
  const auto it = names.find(name);
  if (it != names.end())
    return it->second;
  std::string out = escape(name);
  std::replace(out.begin(), out.end(), '_', ' ');
  return out;
}

// A summary table as an HTML table, its columns in their order.
std::string html_table(const json::Value &t) {
  std::ostringstream o;
  const json::Array &names = t["columns"].as_array();
  std::size_t rows = 0;
  for (const json::Value &n : names)
    rows = std::max(rows, column(t, n.as_string()).size());
  o << "<table><thead><tr>";
  for (const json::Value &n : names)
    o << "<th>" << pretty(n.as_string()) << "</th>";
  o << "</tr></thead><tbody>";
  for (std::size_t r = 0; r < rows; ++r) {
    o << "<tr>";
    for (const json::Value &n : names) {
      const json::Array c = column(t, n.as_string());
      o << "<td>" << (r < c.size() ? text_of(c[r]) : "") << "</td>";
    }
    o << "</tr>";
  }
  o << "</tbody></table>";
  return o.str();
}

// Plotly draws its text itself and knows only a few of HTML's entities (and
// its tags <sub>, <sup>, <b>, <br>): the rest as the characters, the page
// being UTF-8.
std::string unentity(std::string s) {
  static const std::vector<std::pair<std::string, std::string>> map = {
      {"&lang;", "\u27e8"},   {"&rang;", "\u27e9"},  {"&sigma;", "\u03c3"},
      {"&frac12;", "\u00bd"}, {"&Aring;", "\u00c5"}, {"&alpha;", "\u03b1"},
      {"&beta;", "\u03b2"},   {"&gamma;", "\u03b3"}, {"&deg;", "\u00b0"}};
  for (const auto &[from, to] : map)
    for (std::size_t at = s.find(from); at != std::string::npos;
         at = s.find(from, at + to.size()))
      s.replace(at, from.size(), to);
  return s;
}

// A plotly figure.
class Figures {
public:
  //: `y_span`: the least the y axis shows, centred on the data -- a cell
  //: constant over the scan is a flat line, not its rounding drawn large
  //: (Graeme: the angles, all 90.000, looked wild). `x_axis`: more of the x
  //: axis's layout, ticks for one.
  std::string plot(const std::string &title, const std::string &x_title,
                   const std::string &y_title, const json::Array &traces,
                   bool reverse_x = false, double y_span = 0.0,
                   const json::Object &x_axis = {}, bool y_log = false,
                   const std::string &y2_title = "",
                   const json::Object &y_axis = {},
                   const json::Object &y2_axis = {}) {
    const std::string id = "plot" + std::to_string(++count_);
    json::Object layout;
    layout["title"] = json::Value(unentity(title));
    json::Object x, y;
    x["title"] = json::Value(unentity(x_title));
    if (reverse_x)
      x["autorange"] = json::Value(std::string("reversed"));
    for (const auto &[k, v] : x_axis)
      x[k] = v;
    y["title"] = json::Value(unentity(y_title));
    // R factors on a log scale: they span decades, and one wild shell had
    // flattened the rest (Graeme).
    if (y_log)
      y["type"] = json::Value(std::string("log"));
    for (const auto &[k, v] : y_axis)
      y[k] = v;
    // A second axis on the right, for the traces marked y2: multiplicity
    // beside completeness.
    if (!y2_title.empty()) {
      json::Object y2{{"title", json::Value(unentity(y2_title))},
                      {"overlaying", json::Value(std::string("y"))},
                      {"side", json::Value(std::string("right"))}};
      for (const auto &[k, v] : y2_axis)
        y2[k] = v;
      layout["yaxis2"] = json::Value(std::move(y2));
    }
    if (y_span > 0.0) {
      double lo = HUGE_VAL, hi = -HUGE_VAL;
      for (const json::Value &t : traces)
        if (t["y"].is_array())
          for (const json::Value &v : t["y"].as_array())
            if (v.is_number()) {
              lo = std::min(lo, v.as_number());
              hi = std::max(hi, v.as_number());
            }
      if (hi >= lo && hi - lo < y_span) {
        const double mid = 0.5 * (hi + lo);
        y["range"] = json::Value(json::Array{json::Value(mid - 0.5 * y_span),
                                             json::Value(mid + 0.5 * y_span)});
      }
    }
    layout["xaxis"] = json::Value(std::move(x));
    layout["yaxis"] = json::Value(std::move(y));
    layout["margin"] = json::Value(json::Object{{"t", json::Value(40L)}});
    layout["height"] = json::Value(360L);
    std::ostringstream o;
    o << "<div class=\"plot\" id=\"" << id
      << "\"></div><script>Plotly.newPlot('" << id << "', "
      << json::dump(json::Value(traces), 0) << ", "
      << json::dump(json::Value(std::move(layout)), 0)
      << ", {responsive: true});</script>";
    return o.str();
  }

private:
  int count_ = 0;
};

json::Value trace(const json::Array &x, const json::Array &y,
                  const std::string &name,
                  const std::string &mode = "lines+markers",
                  const std::string &type = "scatter") {
  json::Object t;
  t["x"] = json::Value(x);
  t["y"] = json::Value(y);
  t["name"] = json::Value(unentity(name));
  t["type"] = json::Value(type);
  if (type == "scatter")
    t["mode"] = json::Value(mode);
  return json::Value(std::move(t));
}

// A trace whose points at the axis's limits are drawn whole, not cut by it.
json::Value unclipped(json::Value t) {
  json::Object o = t.as_object();
  o["cliponaxis"] = json::Value(false);
  return json::Value(std::move(o));
}

// A trace drawn against the right-hand axis.
json::Value on_right(json::Value t) {
  json::Object o = t.as_object();
  o["yaxis"] = json::Value(std::string("y2"));
  return json::Value(std::move(o));
}

json::Array counting(std::size_t n, long first) {
  json::Array a;
  for (std::size_t i = 0; i < n; ++i)
    a.push_back(json::Value(first + static_cast<long>(i)));
  return a;
}

std::string figures_list(const json::Value &step) {
  std::ostringstream o;
  bool any = false;
  for (const auto &[name, value] : step.as_object()) {
    // Not what a plot or a table shows: long or nested lists -- each sweep's
    // spots per image had come out as a wall of numbers.
    if (name == "program" || name == "time" || is_table(value))
      continue;
    if (value.is_array()) {
      bool nested = value.as_array().size() > 12;
      for (const json::Value &v : value.as_array())
        nested = nested || v.is_array() || v.is_object();
      if (nested)
        continue;
    }
    if (!any)
      o << "<dl>";
    any = true;
    o << "<dt>" << pretty(name) << "</dt><dd>" << text_of(value) << "</dd>";
  }
  if (any)
    o << "</dl>";
  return o.str();
}

std::string cell_text(const json::Value &cell, const json::Value &esd) {
  if (!cell.is_array() || cell.as_array().size() != 6)
    return text_of(cell);
  const json::Array &c = cell.as_array();
  std::string out;
  for (std::size_t k = 0; k < 6; ++k) {
    const double e = esd.is_array() && esd.as_array().size() == 6
                         ? esd.as_array()[k].as_number()
                         : 0.0;
    out += with_esd(c[k].as_number(), e, k < 3 ? 4 : 3);
    out += k == 2 ? " &Aring;, " : (k == 5 ? "&deg;" : " ");
  }
  return out;
}

// --- The steps -------------------------------------------------------------

std::string section_find(const json::Value &s, Figures &f) {
  std::ostringstream o;
  if (is_table(s["each_sweep"]))
    o << "<h4>Each sweep</h4>" << html_table(s["each_sweep"]);
  // Spots per image, a plot a sweep, two to a row (Graeme: four lines on one
  // plot read as none). A single list, as files before October's 13th have,
  // a plot alone.
  const json::Value &per = s["spots_per_image"];
  if (per.is_array() && !per.as_array().empty()) {
    const json::Array &list = per.as_array();
    const bool several = list.front().is_array();
    const std::size_t n = several ? list.size() : 1;
    o << "<div class=\"pair\">";
    for (std::size_t k = 0; k < n; ++k) {
      const json::Array counts = several ? list[k].as_array() : list;
      const json::Value &firsts = s["first_image"];
      const long first =
          firsts.is_array() && k < firsts.as_array().size()
              ? static_cast<long>(firsts.as_array()[k].as_number())
              : (firsts.is_number() ? static_cast<long>(firsts.as_number())
                                    : 1);
      o << f.plot(
          n > 1 ? "Spots per image, sweep " + std::to_string(k)
                : std::string("Spots per image"),
          "image", "spots",
          {trace(counting(counts.size(), first), counts, "spots", "lines")});
    }
    o << "</div>";
  }
  // By resolution: drawn against 1/d^2, its ticks labelled in d.
  const json::Value &r = s["spots_by_resolution"];
  if (is_table(r)) {
    const json::Array dmax = column(r, "d_max"), dmin = column(r, "d_min"),
                      n = column(r, "spots");
    json::Array x;
    double top = 0.0;
    for (std::size_t k = 0; k < dmin.size(); ++k) {
      const double a = dmax[k].as_number() > 0.0
                           ? 1.0 / (dmax[k].as_number() * dmax[k].as_number())
                           : 0.0;
      const double b = 1.0 / (dmin[k].as_number() * dmin[k].as_number());
      x.push_back(json::Value(0.5 * (a + b)));
      top = std::max(top, b);
    }
    json::Array ticks, labels;
    for (double d : {10.0, 5.0, 4.0, 3.0, 2.5, 2.0, 1.8, 1.6, 1.5, 1.4, 1.3,
                     1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5})
      if (1.0 / (d * d) <= top) {
        ticks.push_back(json::Value(1.0 / (d * d)));
        char text[16];
        std::snprintf(text, sizeof text, "%g", d);
        labels.push_back(json::Value(std::string(text)));
      }
    json::Object axis;
    axis["tickvals"] = json::Value(ticks);
    axis["ticktext"] = json::Value(labels);
    o << f.plot("Spots by resolution, every sweep",
                "d (&Aring;), on a scale of 1/d<sup>2</sup>", "spots",
                {trace(x, n, "spots", "markers", "bar")}, false, 0.0, axis);
  }
  return o.str();
}

std::string section_index(const json::Value &s, Figures &f) {
  std::ostringstream o;
  if (s["cell"].is_array())
    o << "<p class=\"result\">Cell " << cell_text(s["cell"], json::Value())
      << "; indexed " << number(100.0 * s["indexed_fraction"].as_number(), 1)
      << "% of " << text_of(s["reflections"]) << "</p>";
  const json::Value &c = s["candidates"];
  if (is_table(c)) {
    o << "<h4>Candidate solutions</h4><table><thead><tr><th>cell</th><th>"
         "volume (&Aring;<sup>3</sup>)</th><th>indexed</th><th></th></tr></"
         "thead><tbody>";
    const json::Array cells = column(c, "cell"), volume = column(c, "volume"),
                      indexed = column(c, "indexed"),
                      chosen = column(c, "chosen");
    for (std::size_t k = 0; k < cells.size(); ++k)
      o << "<tr><td>" << cell_text(cells[k], json::Value()) << "</td><td>"
        << number(volume[k].as_number(), 0) << "</td><td>"
        << text_of(indexed[k]) << "</td><td>"
        << (k < chosen.size() && chosen[k].as_number() > 0 ? "chosen" : "")
        << "</td></tr>";
    o << "</tbody></table>";
  }
  const json::Value &k = s["cycles"];
  if (is_table(k)) {
    const json::Array x = column(k, "cycle");
    o << f.plot("RMSDs by macrocycle", "cycle", "RMSD",
                {trace(x, column(k, "rmsd_x_px"), "x (px)"),
                 trace(x, column(k, "rmsd_y_px"), "y (px)"),
                 trace(x, column(k, "rmsd_z_images"), "z (images)")});
  }
  return o.str();
}

std::string section_refine(const json::Value &s, Figures &f,
                           const ExperimentList *list) {
  std::ostringstream o;
  const json::Value &c = s["crystals"];
  if (is_table(c)) {
    const json::Array cells = column(c, "cell"), esd = column(c, "cell_esd"),
                      exp = column(c, "experiment"),
                      turn = column(c, "orientation_turn_deg");
    for (std::size_t k = 0; k < cells.size(); ++k)
      o << "<p class=\"result\">Cell"
        << (cells.size() > 1 ? " " + text_of(exp[k]) : "") << ": "
        << cell_text(cells[k], k < esd.size() ? esd[k] : json::Value())
        << (k < turn.size() && turn[k].as_number() > 0
                ? "; the orientation turns " + number(turn[k].as_number(), 3) +
                      "&deg; over the scan"
                : "")
        << "</p>";
  }
  const json::Value &k = s["cycles"];
  if (is_table(k)) {
    const json::Array x = counting(column(k, "cycle").size(), 1);
    o << f.plot("RMSDs, cycle by cycle (scan-static, then scan-varying)",
                "cycle", "RMSD",
                {trace(x, column(k, "rmsd_x_px"), "x (px)"),
                 trace(x, column(k, "rmsd_y_px"), "y (px)"),
                 trace(x, column(k, "rmsd_z_images"), "z (images)")});
  }
  // The cell over the scan, a crystal at a time, lengths beside angles: as
  // refinement left it, which its step records -- the file's own crystals
  // have since been constrained by symmetry, a cubic cell's angles 90.000.
  // Files written before October's 13th have only the file's.
  struct Series {
    std::string tag;
    json::Array x, a, b, c, alpha, beta, gamma;
  };
  std::vector<Series> series;
  const json::Value &over = s["cell_over_scan"];
  if (over.is_array()) {
    for (const json::Value &one : over.as_array())
      series.push_back(
          {"experiment " + text_of(one["experiment"]), column(one, "image"),
           column(one, "a"), column(one, "b"), column(one, "c"),
           column(one, "alpha"), column(one, "beta"), column(one, "gamma")});
  } else if (list) {
    for (std::size_t i = 0; i < list->size(); ++i) {
      const Experiment &e = (*list)[i];
      if (!e.crystal || e.crystal->A_points.size() < 2)
        continue;
      Series one;
      one.tag = "experiment " + std::to_string(i) + ", as the file has it";
      for (std::size_t p = 0; p < e.crystal->A_points.size(); ++p) {
        Crystal at = *e.crystal;
        at.A = e.crystal->A_points[p];
        const UnitCell u = at.cell();
        one.x.push_back(json::Value(static_cast<double>(e.scan.first_image) +
                                    static_cast<double>(p) - 0.5));
        one.a.push_back(json::Value(u.a));
        one.b.push_back(json::Value(u.b));
        one.c.push_back(json::Value(u.c));
        one.alpha.push_back(json::Value(u.alpha));
        one.beta.push_back(json::Value(u.beta));
        one.gamma.push_back(json::Value(u.gamma));
      }
      series.push_back(std::move(one));
    }
  }
  if (!series.empty()) {
    o << "<h4>The cell over the scan</h4><div class=\"pair\">";
    for (const Series &one : series) {
      o << f.plot("Cell lengths, " + one.tag, "image", "&Aring;",
                  {trace(one.x, one.a, "a", "lines"),
                   trace(one.x, one.b, "b", "lines"),
                   trace(one.x, one.c, "c", "lines")},
                  false, 0.01);
      o << f.plot("Cell angles, " + one.tag, "image", "degrees",
                  {trace(one.x, one.alpha, "&alpha;", "lines"),
                   trace(one.x, one.beta, "&beta;", "lines"),
                   trace(one.x, one.gamma, "&gamma;", "lines")},
                  false, 0.1);
    }
    o << "</div>";
  }
  return o.str();
}

std::string section_integrate(const json::Value &s, Figures &f) {
  std::ostringstream o;
  const json::Value &p = s["profile_model"];
  if (is_table(p))
    o << "<h4>Profile model</h4>" << html_table(p);
  // Two plots of one summary by resolution: <I/sigma> fitted and summed;
  // profile correlation and background.
  const auto pair_of = [&](const json::Value &r, const std::string &whose) {
    const json::Array d = column(r, "d_min");
    std::string out;
    out += f.plot("&lang;I/&sigma;&rang; by resolution" + whose,
                  "d<sub>min</sub> (&Aring;)", "&lang;I/&sigma;&rang;",
                  {trace(d, column(r, "isigma_prf"), "profile fitted"),
                   trace(d, column(r, "isigma_sum"), "summed")},
                  true);
    out += f.plot("Profile correlation and background" + whose,
                  "d<sub>min</sub> (&Aring;)", "",
                  {trace(d, column(r, "cc_prf"), "profile CC"),
                   trace(d, column(r, "ibg"), "background (counts/px)")},
                  true);
    return out;
  };
  // Each sweep's own, a run each and joined (Graeme: four runs, one plot):
  // its two plots and its overall row; then all together.
  const json::Value &each = s["each_sweep"];
  if (each.is_array() && !each.as_array().empty()) {
    o << "<div class=\"pair\">";
    for (const json::Value &sw : each.as_array())
      if (is_table(sw["by_resolution"]))
        o << pair_of(sw["by_resolution"], ", sweep " + text_of(sw["sweep"]));
    o << "</div>";
    const json::Value &first = each.as_array().front()["overall"];
    if (is_table(first)) {
      o << "<h4>Each sweep, overall</h4><table><thead><tr><th>sweep</th>";
      for (const json::Value &n : first["columns"].as_array())
        o << "<th>" << pretty(n.as_string()) << "</th>";
      o << "</tr></thead><tbody>";
      for (const json::Value &sw : each.as_array()) {
        const json::Value &ov = sw["overall"];
        if (!is_table(ov))
          continue;
        o << "<tr><td>" << text_of(sw["sweep"]) << "</td>";
        for (const json::Value &n : ov["columns"].as_array()) {
          const json::Array c = column(ov, n.as_string());
          o << "<td>" << (c.empty() ? "" : text_of(c[0])) << "</td>";
        }
        o << "</tr>";
      }
      o << "</tbody></table>";
    }
  }
  const json::Value &r = s["by_resolution"];
  if (is_table(r)) {
    const bool several = each.is_array() && !each.as_array().empty();
    o << "<div class=\"pair\">" << pair_of(r, several ? ", every sweep" : "")
      << "</div>";
    o << "<h4>By resolution" << (several ? ", every sweep" : "") << "</h4>"
      << html_table(r);
  }
  if (is_table(s["overall"]))
    o << "<h4>Overall</h4>" << html_table(s["overall"]);
  return o.str();
}

std::string section_symmetry(const json::Value &s) {
  std::ostringstream o;
  if (s["space_group"].is_string())
    o << "<p class=\"result\">Laue group " << text_of(s["laue_group"])
      << "; space group <b>" << text_of(s["space_group"]) << "</b>"
      << (s["indistinguishable_by_absences"].is_array() &&
                  !s["indistinguishable_by_absences"].as_array().empty()
              ? " (" + text_of(s["indistinguishable_by_absences"]) +
                    " indistinguishable by absences)"
              : "")
      << "; reindex " << text_of(s["reindex"]) << "</p>";
  if (is_table(s["elements"]))
    o << "<h4>Each symmetry element</h4>" << html_table(s["elements"]);
  if (is_table(s["subgroups"]))
    o << "<h4>The subgroups</h4>" << html_table(s["subgroups"]);
  return o.str();
}

// The merging summary as DIALS lays it out: statistics down, shells across.
std::string merging_summary(const json::Value &t) {
  static const std::vector<std::pair<std::string, int>> rows = {
      {"d_min", 2},
      {"d_max", 2},
      {"completeness", 1},
      {"multiplicity", 1},
      {"i_over_sigma", 1},
      {"rmerge", 3},
      {"rmeas", 3},
      {"rpim", 3},
      {"cc_half", 3},
      {"anom_completeness", 1},
      {"anom_multiplicity", 1},
      {"cc_anom", 3},
      {"observations", 0},
      {"unique", 0}};
  const json::Array shells = column(t, "shell");
  std::ostringstream o;
  o << "<table class=\"summary\"><thead><tr><th></th>";
  for (const json::Value &s : shells)
    o << "<th>" << text_of(s) << "</th>";
  o << "</tr></thead><tbody>";
  for (const auto &[name, places] : rows) {
    const json::Array c = column(t, name);
    if (c.empty())
      continue;
    o << "<tr><th>" << pretty(name) << "</th>";
    for (const json::Value &v : c)
      o << "<td>" << number(v.as_number(), places) << "</td>";
    o << "</tr>";
  }
  o << "</tbody></table>";
  return o.str();
}

std::string section_scale(const json::Value &s, Figures &f,
                          const ExperimentList *list, const Table *table) {
  std::ostringstream o;
  if (is_table(s["merging_summary"]))
    o << merging_summary(s["merging_summary"]);
  if (s["two_theta_cell"].is_array())
    o << "<p class=\"result\">Cell from two-theta: "
      << cell_text(s["two_theta_cell"], s["two_theta_cell_esd"]) << "; "
      << text_of(s["two_theta_reflections"]) << " reflections, "
      << text_of(s["two_theta_system"]) << ", RMSD "
      << number(s["two_theta_rmsd_deg"].as_number(), 4) << "&deg;</p>";
  o << "<p class=\"result\">Error model a = "
    << number(s["error_model_a"].as_number(), 4)
    << ", b = " << number(s["error_model_b"].as_number(), 4) << "; "
    << text_of(s["outliers"]) << " outliers";
  if (s["resolution_cc_half"].is_number())
    o << "; resolution from CC<sub>&frac12;</sub> "
      << number(s["resolution_cc_half"].as_number(), 2) << " &Aring;";
  o << "</p>";
  const json::Value &m = s["merging_by_resolution"];
  if (is_table(m)) {
    const json::Array d = column(m, "d_min");
    o << "<div class=\"pair\">";
    o << f.plot("CC<sub>&frac12;</sub> by resolution",
                "d<sub>min</sub> (&Aring;)", "CC<sub>&frac12;</sub>",
                {trace(d, column(m, "cc_half"), "CC&frac12;"),
                 trace(d, column(m, "cc_anom"), "CC anom")},
                true);
    o << f.plot("&lang;I/&sigma;&rang; by resolution",
                "d<sub>min</sub> (&Aring;)", "&lang;I/&sigma;&rang;",
                {trace(d, column(m, "i_over_sigma"), "&lang;I/&sigma;&rang;")},
                true);
    o << f.plot(
        "Completeness and multiplicity by resolution",
        "d<sub>min</sub> (&Aring;)", "completeness (%)",
        {unclipped(trace(d, column(m, "completeness"), "completeness (%)")),
         on_right(trace(d, column(m, "multiplicity"), "multiplicity"))},
        true, 0.0, {}, false, "multiplicity",
        // Completeness 0 to 100 per cent, multiplicity from zero:
        // easier read than plotly's own ranges (Graeme).
        {{"range",
          json::Value(json::Array{json::Value(0.0), json::Value(100.0)})}},
        {{"rangemode", json::Value(std::string("tozero"))}});
    o << f.plot("R factors by resolution", "d<sub>min</sub> (&Aring;)",
                "R (log scale)",
                {trace(d, column(m, "rmerge"), "R merge"),
                 trace(d, column(m, "rmeas"), "R meas"),
                 trace(d, column(m, "rpim"), "R pim")},
                true, 0.0, {}, true);
    o << "</div><h4>Merging statistics by resolution</h4>" << html_table(m);
  }
  // By image, from the reflections: each sweep's mean scale and <I/sigma>
  // over fifty blocks of its images.
  if (list && table && table->has("inverse_scale_factor") &&
      table->has("intensity.scale.value") &&
      table->has("intensity.scale.variance") && table->has("xyzobs.px.value") &&
      table->has("id") && table->has("flags")) {
    const Column &g = table->at("inverse_scale_factor"),
                 &I = table->at("intensity.scale.value"),
                 &V = table->at("intensity.scale.variance"),
                 &z = table->at("xyzobs.px.value"), &id = table->at("id"),
                 &flags = table->at("flags");
    json::Array scale_traces, isig_traces;
    for (std::size_t e = 0; e < list->size(); ++e) {
      const Experiment &x = (*list)[e];
      if (!x.crystal)
        continue;
      const double first = static_cast<double>(x.scan.first_image);
      const double n = static_cast<double>(x.scan.num_images());
      const std::size_t blocks = 50;
      std::vector<double> sg(blocks, 0.0), si(blocks, 0.0);
      std::vector<std::size_t> count(blocks, 0);
      for (std::size_t r = 0; r < table->nrows; ++r) {
        if (static_cast<std::size_t>(id.integer(r)) != e ||
            ((flags.integer(r) >> 26) & 1) == 0 || !(V.reals[r] > 0.0))
          continue;
        // Images into the sweep: z counts from the scan's first image less one,
        // not from zero -- a sweep starting at image 151 had fallen into one
        // block.
        const double frame = z.reals[r * 3 + 2] - (first - 1.0);
        const std::size_t b = std::min<std::size_t>(
            blocks - 1,
            static_cast<std::size_t>(std::fmax(0.0, frame) / n * blocks));
        sg[b] += g.reals[r];
        si[b] += I.reals[r] / std::sqrt(V.reals[r]);
        ++count[b];
      }
      json::Array xs, gs, is;
      for (std::size_t b = 0; b < blocks; ++b) {
        if (count[b] == 0)
          continue;
        xs.push_back(json::Value(first + (static_cast<double>(b) + 0.5) * n /
                                             static_cast<double>(blocks)));
        gs.push_back(json::Value(sg[b] / static_cast<double>(count[b])));
        is.push_back(json::Value(si[b] / static_cast<double>(count[b])));
      }
      const std::string tag = "experiment " + std::to_string(e);
      scale_traces.push_back(trace(xs, gs, tag));
      isig_traces.push_back(trace(xs, is, tag));
    }
    if (!scale_traces.empty()) {
      o << f.plot("Inverse scale factor by image", "image", "inverse scale",
                  scale_traces);
      o << f.plot("&lang;I/&sigma;&rang; by image", "image",
                  "&lang;I/&sigma;&rang;", isig_traces);
    }
  }
  return o.str();
}

std::string section_generic(const json::Value &s) {
  std::ostringstream o;
  for (const auto &[name, value] : s.as_object())
    if (is_table(value))
      o << "<h4>" << pretty(name) << "</h4>" << html_table(value);
  return o.str();
}

const char *kStyle = R"(
body { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
       margin: 2em auto; max-width: 1100px; color: #222; padding: 0 1em; }
h1 { font-weight: 400; border-bottom: 2px solid #3a6ea5; padding-bottom: .3em; }
h2 { font-weight: 400; color: #3a6ea5; margin-top: 2em;
     border-bottom: 1px solid #ccc; padding-bottom: .2em; }
h4 { margin-bottom: .3em; }
table { border-collapse: collapse; font-size: 0.85em; margin: .5em 0 1em;
        font-variant-numeric: tabular-nums; }
th, td { padding: .25em .6em; text-align: right; border-bottom: 1px solid #eee; }
thead th { border-bottom: 2px solid #999; }
table.summary th { text-align: left; }
dl { display: grid; grid-template-columns: max-content auto; gap: .2em 1em; }
dt { font-weight: 600; } dd { margin: 0; }
p.result { font-size: 1.05em; background: #f3f6fa; padding: .5em .8em;
           border-left: 3px solid #3a6ea5; }
.when { color: #777; font-size: .85em; }
.plot { width: 100%; min-height: 360px; }
div.pair { display: grid; grid-template-columns: 1fr 1fr; gap: 1em; }
@media (max-width: 800px) { div.pair { grid-template-columns: 1fr; } }
)";

} // namespace

int run_program(int argc, char **argv) {
  const std::set<std::string> known = {"-o", "--help"};
  const std::set<std::string> takes_value = {"-o"};
  Arguments args = parse_arguments(argc, argv, known, takes_value);
  if (args.has("--help") || args.positional.size() != 1) {
    usage(argv[0]);
    return args.has("--help") ? 0 : 2;
  }
  const std::string input = args.positional[0];
  const std::string output = args.value("-o", "mxi_report.html");
  try {
    if (!rflx::is_hdf5(input))
      throw std::runtime_error(input + " is not a .rflx");
    const auto summary = rflx::read_summary(input);
    if (!summary)
      throw std::runtime_error(input +
                               " has no processing summary: written before "
                               "mxi kept one, or by mxi_convert from a pair");
    std::optional<ExperimentList> list;
    try {
      list = read_experiments(input);
    } catch (const std::exception &) {
    }
    std::optional<Table> table;
    try {
      table = read_reflections(input);
    } catch (const std::exception &) {
    }

    const std::vector<json::Value> steps = Summary::steps(*summary);
    Figures figures;
    std::ostringstream body;
    body << "<h1>mxi processing report: " << escape(input) << "</h1>";
    // How the data were processed: the steps in order.
    body << "<h2>The processing</h2><table><thead><tr><th>step</th>"
            "<th>program</th><th>when</th></tr></thead><tbody>";
    for (std::size_t k = 0; k < steps.size(); ++k)
      body << "<tr><td>" << k + 1 << "</td><td style=\"text-align:left\">"
           << text_of(steps[k]["program"]) << "</td><td>"
           << text_of(steps[k]["time"]) << "</td></tr>";
    body << "</tbody></table>";
    if (list) {
      body << "<h4>Experiments</h4><table><thead><tr><th>experiment</th>"
              "<th>images</th><th>wavelength (&Aring;)</th><th>crystal</th>"
              "</tr></thead><tbody>";
      for (std::size_t i = 0; i < list->size(); ++i) {
        const Experiment &e = (*list)[i];
        std::string crystal = "the sweep's own, its unindexed spots";
        if (e.crystal) {
          const UnitCell c = e.crystal->cell();
          std::string group = e.crystal->space_group_hall;
          try {
            group = SpaceGroup::from_hall(group).name();
          } catch (const std::exception &) {
          }
          char text[160];
          std::snprintf(text, sizeof text,
                        "%.3f %.3f %.3f &Aring;, %.2f %.2f %.2f&deg;; %s", c.a,
                        c.b, c.c, c.alpha, c.beta, c.gamma, group.c_str());
          crystal = text;
        }
        body << "<tr><td>" << i << "</td><td>" << e.scan.first_image
             << "&ndash;" << e.scan.last_image << "</td><td>"
             << number(e.beam.wavelength, 5)
             << "</td><td style=\"text-align:left\">" << crystal
             << "</td></tr>";
      }
      body << "</tbody></table>";
    }
    static const std::map<std::string, std::string> titles = {
        {"mxi_import", "Import"},
        {"mxi_find", "Spot finding"},
        {"mxi_index", "Indexing"},
        {"mxi_refine", "Refinement"},
        {"mxi_integrate", "Integration"},
        {"mxi_symmetry", "Symmetry"},
        {"mxi_scale", "Scaling and merging"}};
    // The last step first is what most readers want: scaling's, if there is
    // one, then each step in the order it ran.
    for (const json::Value &s : steps) {
      const std::string program = s["program"].as_string();
      const auto t = titles.find(program);
      body << "<h2>" << (t != titles.end() ? t->second : escape(program))
           << " <span class=\"when\">" << escape(program) << ", "
           << text_of(s["time"]) << "</span></h2>";
      body << figures_list(s);
      if (program == "mxi_find")
        body << section_find(s, figures);
      else if (program == "mxi_index")
        body << section_index(s, figures);
      else if (program == "mxi_refine")
        body << section_refine(s, figures, list ? &*list : nullptr);
      else if (program == "mxi_integrate")
        body << section_integrate(s, figures);
      else if (program == "mxi_symmetry")
        body << section_symmetry(s);
      else if (program == "mxi_scale")
        body << section_scale(s, figures, list ? &*list : nullptr,
                              table ? &*table : nullptr);
      else
        body << section_generic(s);
    }

    std::ofstream out(output);
    if (!out)
      throw std::runtime_error("cannot write " + output);
    out << "<!DOCTYPE html><html><head><meta charset=\"utf-8\"><title>mxi "
           "report: "
        << escape(input) << "</title><script src=\"" << kPlotly
        << "\"></script><style>" << kStyle << "</style></head><body>"
        << body.str() << "</body></html>\n";
    std::printf("Wrote %s: %zu steps of processing\n", output.c_str(),
                steps.size());
    return 0;
  } catch (const std::exception &error) {
    std::fprintf(stderr, "mxi_report: %s\n", error.what());
    return 1;
  }
}

} // namespace mxi

int main(int argc, char **argv) { return mxi::run_program(argc, argv); }
