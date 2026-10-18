// Experiment lists split apart and joined together: see combine.hh.

#include "combine.hh"

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

#include "shoebox.hh"
#include "views.hh"

namespace mxi {

namespace {

// The arrays of models an experiment refers to by index.
const char *const kModels[] = {"beam",          "crystal",  "detector",
                               "goniometer",    "imageset", "profile",
                               "scaling_model", "scan"};

std::size_t array_size(const json::Value &document, const char *key) {
  return document.is_object() && document[key].is_array()
             ? document[key].as_array().size()
             : 0;
}

// The table with only `names` of its decoded columns, in its row order.
Table with_columns(const Table &table, const std::vector<std::string> &names) {
  Table out;
  out.nrows = table.nrows;
  out.version = table.version;
  out.identifiers = table.identifiers;
  for (const std::string &name : names) {
    const Column &from = table.at(name);
    Column &to = from.integral ? out.int_column(name, from.type, from.width)
                               : out.real_column(name, from.type, from.width);
    to.ints = from.ints;
    to.reals = from.reals;
  }
  return out;
}

} // namespace

ExperimentsAndReflections select_experiments(const json::Value &document,
                                             const Table &reflections,
                                             const std::vector<bool> &keep) {
  const json::Array &all = document["experiment"].as_array();
  std::vector<bool> chosen(all.size(), false);
  std::set<std::pair<long, long>> sweeps; // image set, scan
  for (std::size_t e = 0; e < all.size() && e < keep.size(); ++e)
    if (keep[e]) {
      chosen[e] = true;
      const json::Value &x = all[e];
      if (x["imageset"].is_number() && x["scan"].is_number())
        sweeps.insert({static_cast<long>(x["imageset"].as_number()),
                       static_cast<long>(x["scan"].as_number())});
    }
  for (std::size_t e = 0; e < all.size(); ++e)
    if (!all[e]["crystal"].is_number() && all[e]["imageset"].is_number() &&
        all[e]["scan"].is_number() &&
        sweeps.count({static_cast<long>(all[e]["imageset"].as_number()),
                      static_cast<long>(all[e]["scan"].as_number())}))
      chosen[e] = true;
  std::vector<std::size_t> rows;
  if (reflections.has("id")) {
    const Column &id = reflections.at("id");
    for (std::size_t r = 0; r < reflections.nrows; ++r) {
      const std::int64_t e = id.integer(r);
      if (e >= 0 && static_cast<std::size_t>(e) < chosen.size() &&
          chosen[static_cast<std::size_t>(e)])
        rows.push_back(r);
    }
  }
  ExperimentsAndReflections out;
  out.reflections = select_rows_with_shoeboxes(reflections, rows);
  out.document = keep_experiments(document, chosen, &out.reflections);
  return out;
}

Combined combine_experiments(const std::vector<json::Value> &documents,
                             const std::vector<Table> &tables) {
  if (documents.empty() || documents.size() != tables.size())
    throw std::runtime_error("combine: an experiment list and its reflections "
                             "for each input");
  // The columns every table has, of the same type and width; shoeboxes and any
  // other undecoded column are left out (they cannot be checked by row).
  std::vector<std::string> common;
  std::set<std::string> dropped;
  for (const std::string &name : tables.front().names()) {
    const Column &first = tables.front().at(name);
    bool everywhere = true;
    for (const Table &t : tables)
      if (!t.has(name) || t.at(name).type != first.type ||
          t.at(name).width != first.width)
        everywhere = false;
    (everywhere ? common.push_back(name) : (void)dropped.insert(name));
  }
  for (const Table &t : tables) {
    for (const std::string &name : t.names())
      if (std::find(common.begin(), common.end(), name) == common.end())
        dropped.insert(name);
    for (const auto &entry : t.opaque())
      if (entry.first != "shoebox")
        dropped.insert(entry.first);
  }
  // The shoeboxes kept where every input has them -- integration's profile
  // model wants them -- and left out, said, where only some do.
  std::size_t with_boxes = 0;
  for (const Table &t : tables)
    with_boxes += t.opaque().count("shoebox");
  const bool keep_boxes = with_boxes == tables.size();
  if (with_boxes > 0 && !keep_boxes)
    dropped.insert("shoebox");

  json::Object out = documents.front().as_object();
  std::map<std::string, json::Array> models;
  for (const char *key : kModels)
    models[key] = json::Array();
  json::Array experiments;
  std::map<std::string, std::size_t> seen; // identifier, input
  std::vector<Table> renumbered;
  for (std::size_t k = 0; k < documents.size(); ++k) {
    const json::Value &document = documents[k];
    if (!document.is_object() || !document["experiment"].is_array())
      throw std::runtime_error("combine: input " + std::to_string(k) +
                               " has no experiment list");
    std::map<std::string, long> offset;
    for (const char *key : kModels) {
      offset[key] = static_cast<long>(models[key].size());
      for (std::size_t i = 0; i < array_size(document, key); ++i)
        models[key].push_back(document[key].as_array()[i]);
    }
    const long first = static_cast<long>(experiments.size());
    for (const json::Value &e : document["experiment"].as_array()) {
      json::Object x = e.as_object();
      for (const char *key : kModels)
        if (x.count(key) && x[key].is_number())
          x[key] =
              json::Value(static_cast<long>(x[key].as_number()) + offset[key]);
      if (x.count("identifier") && x["identifier"].is_string()) {
        const std::string id = x["identifier"].as_string();
        if (!id.empty()) {
          const auto [where, fresh] = seen.emplace(id, k);
          if (!fresh)
            throw std::runtime_error(
                "combine: the same experiment twice -- identifier " + id +
                " in inputs " + std::to_string(where->second) + " and " +
                std::to_string(k));
        }
      }
      experiments.push_back(json::Value(std::move(x)));
    }
    Table t = with_columns(tables[k], common);
    if (keep_boxes)
      t.set_opaque("shoebox", tables[k].opaque().at("shoebox"));
    if (t.has("id")) {
      Column &id = t.modify_int_column("id", "int", 1);
      for (std::int64_t &i : id.ints)
        if (i >= 0)
          i += first;
    }
    std::map<std::size_t, std::string> identifiers;
    for (const auto &[i, name] : t.identifiers)
      identifiers[i + static_cast<std::size_t>(first)] = name;
    t.identifiers = identifiers;
    renumbered.push_back(std::move(t));
  }
  for (auto &[key, array] : models)
    out[key] = json::Value(std::move(array));
  out["experiment"] = json::Value(std::move(experiments));

  Combined combined;
  combined.document = json::Value(std::move(out));
  combined.reflections = concat_rows_with_shoeboxes(renumbered);
  combined.dropped.assign(dropped.begin(), dropped.end());
  return combined;
}

// The experiments in groups that must be refined together: those on one sweep
// (one image set and scan), or sharing a detector, a beam or a crystal model.
// Every other pair shares nothing refinement fits, so each group is refined on
// its own (Graeme: with no shared models, one experiment at a time, not one
// huge joint refinement -- 16 crystals' data sets refined together had been
// the one joined step left in question).
std::vector<std::vector<std::size_t>> independent_groups(const json::Value &doc,
                                                         bool by_models) {
  const json::Array &all = doc["experiment"].as_array();
  const std::size_t n = all.size();
  std::vector<std::size_t> parent(n);
  for (std::size_t i = 0; i < n; ++i)
    parent[i] = i;
  const auto find = [&](std::size_t i) {
    while (parent[i] != i)
      i = parent[i] = parent[parent[i]];
    return i;
  };
  const auto join = [&](std::size_t a, std::size_t b) {
    parent[find(a)] = find(b);
  };
  const auto index_of = [](const json::Value &e, const char *key) {
    return e[key].is_number() ? static_cast<long long>(e[key].as_number())
                              : -1LL;
  };
  std::map<std::pair<long long, long long>, std::size_t> sweep_first;
  std::map<std::pair<std::string, long long>, std::size_t> model_first;
  for (std::size_t i = 0; i < n; ++i) {
    const json::Value &e = all[i];
    const auto sweep =
        std::make_pair(index_of(e, "imageset"), index_of(e, "scan"));
    if (sweep.first >= 0 || sweep.second >= 0) {
      const auto [at, fresh] = sweep_first.emplace(sweep, i);
      if (!fresh)
        join(i, at->second);
    }
    if (!by_models)
      continue;
    for (const char *key : {"detector", "beam", "crystal"}) {
      const long long m = index_of(e, key);
      if (m < 0)
        continue;
      const auto [at, fresh] = model_first.emplace(std::make_pair(key, m), i);
      if (!fresh)
        join(i, at->second);
    }
  }
  std::map<std::size_t, std::size_t> group_of_root;
  std::vector<std::vector<std::size_t>> groups;
  for (std::size_t i = 0; i < n; ++i) {
    const auto [at, fresh] = group_of_root.emplace(find(i), groups.size());
    if (fresh)
      groups.emplace_back();
    groups[at->second].push_back(i);
  }
  return groups;
}

} // namespace mxi
