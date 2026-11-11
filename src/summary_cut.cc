#include "summary_cut.hh"

#include <algorithm>
#include <cstdint>
#include <map>
#include <set>
#include <string>

namespace mxi {

namespace {

// Each experiment's sweep: its image set, numbered in the order image sets
// first appear among the experiments (counted over those marked).
std::vector<std::int64_t> sweeps_of(const json::Value &document,
                                    const std::vector<bool> *only) {
  std::vector<std::int64_t> out;
  std::map<std::int64_t, std::int64_t> seen;
  const json::Value &all = document["experiment"];
  for (std::size_t i = 0; i < all.size(); ++i) {
    const json::Value &set = all[i]["imageset"];
    if (!set.is_number() || (only && !(*only)[i])) {
      out.push_back(-1);
      continue;
    }
    const std::int64_t key = static_cast<std::int64_t>(set.as_number());
    auto [at, fresh] =
        seen.emplace(key, static_cast<std::int64_t>(seen.size()));
    (void)fresh;
    out.push_back(at->second);
  }
  return out;
}

bool is_table(const json::Value &v) {
  return v.is_object() && v["columns"].is_array();
}

std::int64_t number_of(const json::Value &v) {
  return v.is_number() ? static_cast<std::int64_t>(v.as_number()) : -1;
}

// A table's rows where keep(row) is true; column `renumber`, if named,
// mapped through `to`.
json::Value rows_of(const json::Value &table, const std::vector<bool> &keep,
                    const std::string &renumber,
                    const std::map<std::int64_t, std::int64_t> &to) {
  json::Object out;
  out["columns"] = table["columns"];
  for (const json::Value &name : table["columns"].as_array()) {
    const std::string column = name.as_string();
    json::Array kept;
    const json::Value &values = table[column];
    for (std::size_t r = 0; r < values.size() && r < keep.size(); ++r)
      if (keep[r]) {
        if (column == renumber && to.count(number_of(values[r])))
          kept.push_back(json::Value(to.at(number_of(values[r]))));
        else
          kept.push_back(values[r]);
      }
    out[column] = json::Value(std::move(kept));
  }
  return json::Value(std::move(out));
}

} // namespace

json::Value cut_summary(const json::Value &summary, const json::Value &document,
                        const std::vector<bool> &keep, const json::Value *own) {
  // Old to new: experiments by rank among those kept; sweeps by first
  // appearance among them.
  const std::vector<std::int64_t> parent_sweep = sweeps_of(document, nullptr);
  const std::vector<std::int64_t> kept_sweep = sweeps_of(document, &keep);
  std::map<std::int64_t, std::int64_t> experiment, sweep;
  std::int64_t n = 0;
  for (std::size_t i = 0; i < keep.size(); ++i)
    if (keep[i]) {
      experiment[static_cast<std::int64_t>(i)] = n++;
      if (i < parent_sweep.size() && parent_sweep[i] >= 0)
        sweep[parent_sweep[i]] = kept_sweep[i];
    }
  std::int64_t parent_sweeps = 0;
  for (std::int64_t s : parent_sweep)
    parent_sweeps = std::max(parent_sweeps, s + 1);
  // A list aligned to the parent's sweeps, its entries those kept, in the
  // new order.
  const auto by_position = [&](const json::Array &list) {
    std::vector<json::Value> placed(sweep.size());
    for (const auto &[from, to] : sweep)
      if (from < static_cast<std::int64_t>(list.size()))
        placed[static_cast<std::size_t>(to)] =
            list[static_cast<std::size_t>(from)];
    return json::Value(json::Array(placed.begin(), placed.end()));
  };

  json::Array steps;
  for (const json::Value &step : summary["steps"].as_array()) {
    const std::string program =
        step["program"].is_string() ? step["program"].as_string() : "";
    if (program == "mxi_scale")
      continue; // every crystal's statistics: scaled afresh
    if (program == "mxi_symmetry" || program == "mxi_correlation" ||
        !step.is_object()) {
      steps.push_back(step);
      continue;
    }
    json::Object out;
    json::Array left_out;
    for (const auto &[key, value] : step.as_object()) {
      if (key == "program" || key == "time") {
        out[key] = value;
      } else if (is_table(value) && value["sweep"].is_array()) {
        std::vector<bool> rows;
        for (const json::Value &s : value["sweep"].as_array())
          rows.push_back(sweep.count(number_of(s)) > 0);
        out[key] = rows_of(value, rows, "sweep", sweep);
      } else if (key == "sweeps" && is_table(value) &&
                 value[value["columns"][0].as_string()].size() ==
                     static_cast<std::size_t>(parent_sweeps)) {
        // Import's table of sweeps, a row a sweep in order.
        json::Object reordered;
        reordered["columns"] = value["columns"];
        for (const json::Value &name : value["columns"].as_array())
          reordered[name.as_string()] =
              by_position(value[name.as_string()].as_array());
        out[key] = json::Value(std::move(reordered));
      } else if ((key == "spots_per_image" || key == "first_image") &&
                 value.is_array() &&
                 value.size() == static_cast<std::size_t>(parent_sweeps)) {
        out[key] = by_position(value.as_array());
      } else if (value.is_array() && value.size() > 0 && value[0].is_object() &&
                 (value[0].contains("sweep") ||
                  value[0].contains("experiment") ||
                  value[0].contains("experiments"))) {
        json::Array kept;
        for (const json::Value &entry : value.as_array()) {
          json::Object e = entry.as_object();
          if (entry.contains("sweep")) {
            const std::int64_t s = number_of(entry["sweep"]);
            if (!sweep.count(s))
              continue;
            e["sweep"] = json::Value(sweep.at(s));
          } else if (entry.contains("experiment")) {
            const std::int64_t x = number_of(entry["experiment"]);
            if (!experiment.count(x))
              continue;
            e["experiment"] = json::Value(experiment.at(x));
          } else {
            json::Array which;
            for (const json::Value &x : entry["experiments"].as_array())
              if (experiment.count(number_of(x)))
                which.push_back(json::Value(experiment.at(number_of(x))));
            if (which.empty())
              continue;
            e["experiments"] = json::Value(std::move(which));
          }
          kept.push_back(json::Value(std::move(e)));
        }
        out[key] = json::Value(std::move(kept));
      } else {
        left_out.push_back(json::Value(key));
      }
    }
    // Spot finding's totals, the sums of the sweeps kept.
    if (program == "mxi_find" && is_table(out["each_sweep"])) {
      for (const char *total : {"spots", "signal_pixels"}) {
        const json::Value &column = out["each_sweep"][total];
        if (!column.is_array())
          continue;
        std::int64_t sum = 0;
        for (const json::Value &v : column.as_array())
          sum += number_of(v) > 0 ? number_of(v) : 0;
        out[total] = json::Value(sum);
        json::Array still;
        for (const json::Value &k : left_out)
          if (k.as_string() != total)
            still.push_back(k);
        left_out = std::move(still);
      }
    }
    if (!left_out.empty())
      out["whole_run_left_out"] = json::Value(std::move(left_out));
    steps.push_back(json::Value(std::move(out)));
  }
  if (own)
    steps.push_back(*own);
  json::Object result;
  result["steps"] = json::Value(std::move(steps));
  return json::Value(std::move(result));
}

} // namespace mxi
