// The two views of an experiment list: see views.hh.

#include "views.hh"

#include "shoebox.hh"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>

namespace mxi {

namespace {

template <typename T>
std::vector<T> pick(const std::vector<T> &all,
                    const std::vector<std::size_t> &at, std::size_t whole) {
  if (all.size() != whole)
    return {};
  std::vector<T> out;
  for (std::size_t i : at)
    out.push_back(all[i]);
  return out;
}

} // namespace

Lattices lattices_of(const ExperimentList &whole) {
  Lattices l;
  l.list = whole;
  l.list.experiments.clear();
  for (std::size_t i = 0; i < whole.size(); ++i) {
    if (whole[i].crystal) {
      l.position.push_back(i);
      l.list.experiments.push_back(whole[i]);
    } else {
      l.crystal_less.push_back(i);
    }
  }
  const std::size_t n = whole.size();
  l.list.beam_of = pick(whole.beam_of, l.position, n);
  l.list.detector_of = pick(whole.detector_of, l.position, n);
  l.list.sweep_of = pick(whole.sweep_of, l.position, n);
  l.list.profiles = pick(whole.profiles, l.position, n);
  return l;
}

Table lattice_rows(Table table, const Lattices &lattices) {
  if (!lattices.any_crystal_less())
    return table;
  std::vector<long> renumber;
  for (std::size_t k = 0; k < lattices.position.size(); ++k) {
    const std::size_t at = lattices.position[k];
    if (renumber.size() <= at)
      renumber.resize(at + 1, -1);
    renumber[at] = static_cast<long>(k);
  }
  const Column &id = table.at("id");
  std::vector<std::size_t> rows;
  for (std::size_t r = 0; r < table.nrows; ++r) {
    const std::int64_t i = id.integer(r);
    if (i >= 0 && static_cast<std::size_t>(i) < renumber.size() &&
        renumber[static_cast<std::size_t>(i)] >= 0)
      rows.push_back(r);
  }
  // Every row a lattice's -- the usual case, the sweep as imported carried as
  // a crystal-less experiment with none of the rows -- and the table moved in
  // is the view, its ids renumbered in place: a selection here copied the whole
  // table (Graeme's 36000-image sweep, symmetry 22 GB).
  const auto names = table.identifiers; // before any move
  Table out = rows.size() == table.nrows
                  ? std::move(table)
                  : select_rows_with_shoeboxes(table, rows);
  Column &ids = out.modify_int_column("id", "int", 1);
  for (std::int64_t &i : ids.ints)
    i = renumber[static_cast<std::size_t>(i)];
  out.identifiers.clear();
  for (const auto &[i, name] : names)
    if (i < renumber.size() && renumber[i] >= 0)
      out.identifiers[static_cast<std::size_t>(
          renumber[static_cast<std::size_t>(i)])] = name;
  return out;
}

Table crystal_less_rows(const Table &table, const Lattices &lattices) {
  const Column &id = table.at("id");
  std::vector<std::size_t> rows;
  for (std::size_t r = 0; r < table.nrows; ++r) {
    const std::int64_t i = id.integer(r);
    if (i >= 0 &&
        std::find(lattices.crystal_less.begin(), lattices.crystal_less.end(),
                  static_cast<std::size_t>(i)) != lattices.crystal_less.end())
      rows.push_back(r);
  }
  Table out = select_rows_with_shoeboxes(table, rows);
  out.identifiers.clear();
  for (const auto &[i, name] : table.identifiers)
    if (std::find(lattices.crystal_less.begin(), lattices.crystal_less.end(),
                  i) != lattices.crystal_less.end())
      out.identifiers[i] = name;
  return out;
}

ExperimentList whole_again(const ExperimentList &whole,
                           const Lattices &lattices,
                           const ExperimentList &processed) {
  if (!lattices.any_crystal_less())
    return processed;
  if (processed.size() != lattices.position.size())
    throw std::runtime_error(
        "the lattices came back " + std::to_string(processed.size()) +
        " where " + std::to_string(lattices.position.size()) + " went");
  ExperimentList out = processed;
  out.experiments = whole.experiments;
  out.beam_of = whole.beam_of;
  out.detector_of = whole.detector_of;
  out.sweep_of = whole.sweep_of;
  out.source = whole.source;
  for (std::size_t k = 0; k < lattices.position.size(); ++k)
    out.experiments[lattices.position[k]] = processed[k];
  // Each crystal-less experiment: the models it shares with a lattice of its
  // sweep, as that lattice now has them.
  for (std::size_t s : lattices.crystal_less)
    for (std::size_t k = 0; k < lattices.position.size(); ++k) {
      const std::size_t at = lattices.position[k];
      if (whole.sweep(at) != whole.sweep(s))
        continue;
      Experiment &own = out.experiments[s];
      if (whole.detector_group(at) == whole.detector_group(s))
        own.detector = processed[k].detector;
      if (whole.beam_group(at) == whole.beam_group(s))
        own.beam = processed[k].beam;
      own.goniometer = processed[k].goniometer;
      own.scan = processed[k].scan;
      break;
    }
  // Profile models: the lattices' own where they have them; none for the
  // crystal-less.
  out.profiles.clear();
  if (processed.profiles.size() == processed.size() ||
      processed.profile.present) {
    out.profiles.assign(out.size(), ProfileBlock{});
    for (std::size_t k = 0; k < lattices.position.size(); ++k)
      out.profiles[lattices.position[k]] =
          processed.profiles.size() == processed.size() ? processed.profiles[k]
                                                        : processed.profile;
  }
  return out;
}

Table rows_again(Table lattice_rows, const Table &crystal_less_rows,
                 const Lattices &lattices, const ExperimentList &whole) {
  if (!lattices.any_crystal_less())
    return lattice_rows;
  Table back = std::move(lattice_rows);
  Column &ids = back.modify_int_column("id", "int", 1);
  for (std::int64_t &i : ids.ints)
    i = static_cast<std::int64_t>(
        lattices.position.at(static_cast<std::size_t>(i)));
  // The others with exactly the lattices' columns: what the lattices' rows
  // have and theirs do not, as zeros -- an unindexed spot has no prediction --
  // and what theirs have and the lattices' do not, dropped: the lattices'
  // are what a program wrote. Their shoeboxes, if any, kept as they were.
  // No crystal-less rows to put back -- the usual case -- and the lattices'
  // table is the one written, its ids renumbered in place: concatenated with
  // nothing, it was copied.
  if (crystal_less_rows.nrows == 0) {
    back.identifiers.clear();
    for (std::size_t i = 0; i < whole.size(); ++i)
      if (!whole[i].identifier.empty())
        back.identifiers[i] = whole[i].identifier;
    return back;
  }
  Table rest;
  rest.nrows = crystal_less_rows.nrows;
  rest.version = crystal_less_rows.version;
  rest.identifiers = crystal_less_rows.identifiers;
  for (const auto &[name, value] : crystal_less_rows.opaque())
    rest.set_opaque(name, value);
  for (const std::string &name : back.names()) {
    if (crystal_less_rows.has(name)) {
      rest.set(name, crystal_less_rows.at(name));
      continue;
    }
    Column c = back.at(name);
    if (c.integral) {
      c.ints.assign(rest.nrows * c.width, 0);
      c.reals.clear();
    } else {
      c.reals.assign(rest.nrows * c.width, 0.0);
      c.ints.clear();
    }
    rest.set(name, std::move(c));
  }
  Table out = concat_rows_with_shoeboxes({back, rest});
  out.identifiers.clear();
  for (std::size_t i = 0; i < whole.size(); ++i)
    if (!whole[i].identifier.empty())
      out.identifiers[i] = whole[i].identifier;
  return out;
}

ExperimentList add_sweeps(const ExperimentList &given,
                          const ExperimentList &lattices,
                          const std::vector<std::size_t> &sweep_of_lattice,
                          Table *table) {
  const std::size_t n_sweeps = given.size(), n_lattices = lattices.size();
  if (sweep_of_lattice.size() != n_lattices)
    throw std::runtime_error("add_sweeps: a sweep for every lattice wanted");
  ExperimentList whole;
  whole.image_template = given.image_template;
  // The sweeps: as given, crystal-less, with their first lattice's models.
  std::vector<long> first_lattice(n_sweeps, -1);
  for (std::size_t k = 0; k < n_lattices; ++k)
    if (first_lattice[sweep_of_lattice[k]] < 0)
      first_lattice[sweep_of_lattice[k]] = static_cast<long>(k);
  std::vector<std::string> taken;
  for (std::size_t s = 0; s < n_sweeps; ++s) {
    Experiment e = given[s];
    if (first_lattice[s] >= 0) {
      const Experiment &l =
          lattices[static_cast<std::size_t>(first_lattice[s])];
      e.beam = l.beam;
      e.detector = l.detector;
      e.goniometer = l.goniometer;
      e.scan = l.scan;
    }
    e.crystal.reset();
    if (e.identifier.empty())
      e.identifier = new_identifier();
    taken.push_back(e.identifier);
    whole.experiments.push_back(std::move(e));
  }
  for (std::size_t k = 0; k < n_lattices; ++k) {
    Experiment e = lattices[k];
    if (e.identifier.empty() ||
        std::find(taken.begin(), taken.end(), e.identifier) != taken.end())
      e.identifier = new_identifier();
    taken.push_back(e.identifier);
    whole.experiments.push_back(std::move(e));
  }
  // Sharing: a sweep with its lattices' beam and detector, if it has any;
  // each sweep its own sweep, each lattice its sweep's.
  const std::size_t apart = 1000000;
  for (std::size_t s = 0; s < n_sweeps; ++s) {
    const long k = first_lattice[s];
    whole.beam_of.push_back(
        k >= 0 ? lattices.beam_group(static_cast<std::size_t>(k)) : apart + s);
    whole.detector_of.push_back(
        k >= 0 ? lattices.detector_group(static_cast<std::size_t>(k))
               : apart + s);
    whole.sweep_of.push_back(s);
  }
  for (std::size_t k = 0; k < n_lattices; ++k) {
    whole.beam_of.push_back(lattices.beam_group(k));
    whole.detector_of.push_back(lattices.detector_group(k));
    whole.sweep_of.push_back(sweep_of_lattice[k]);
  }
  // The source's experiments, so that each refers to its sweep's image set:
  // the given ones, then a copy of its sweep's for each lattice.
  whole.source = given.source;
  if (given.source.is_object() && given.source["experiment"].is_array() &&
      given.source["experiment"].as_array().size() == n_sweeps) {
    json::Array entries = given.source["experiment"].as_array();
    for (std::size_t k = 0; k < n_lattices; ++k)
      entries.push_back(
          given.source["experiment"].as_array()[sweep_of_lattice[k]]);
    json::Object document = given.source.as_object();
    document["experiment"] = json::Value(std::move(entries));
    whole.source = json::Value(std::move(document));
  }
  if (table != nullptr && table->has("id")) {
    Column &id = table->modify_int_column("id", "int", 1);
    const bool has_miller = table->has("miller_index");
    for (std::size_t r = 0; r < table->nrows; ++r) {
      const std::int64_t old = id.ints[r];
      const bool indexed =
          has_miller && (table->at("miller_index").integer(r, 0) != 0 ||
                         table->at("miller_index").integer(r, 1) != 0 ||
                         table->at("miller_index").integer(r, 2) != 0);
      if (old < 0 || static_cast<std::size_t>(old) >= n_lattices)
        continue;
      id.ints[r] = indexed
                       ? static_cast<std::int64_t>(n_sweeps) + old
                       : static_cast<std::int64_t>(
                             sweep_of_lattice[static_cast<std::size_t>(old)]);
    }
    table->identifiers.clear();
    for (std::size_t i = 0; i < whole.size(); ++i)
      table->identifiers[i] = whole[i].identifier;
  }
  return whole;
}

json::Value dials_form(const json::Value &document, Table *table) {
  if (!document.is_object() || !document["experiment"].is_array())
    return document;
  const json::Array &all = document["experiment"].as_array();
  std::vector<bool> keep(all.size(), false);
  for (std::size_t i = 0; i < all.size(); ++i)
    keep[i] = all[i]["crystal"].is_number();
  return keep_experiments(document, keep, table);
}

json::Value keep_experiments(const json::Value &document,
                             const std::vector<bool> &keep, Table *table) {
  if (!document.is_object() || !document["experiment"].is_array())
    return document;
  const json::Array &all = document["experiment"].as_array();
  std::vector<long> renumber(all.size(), -1);
  json::Array kept;
  for (std::size_t i = 0; i < all.size(); ++i)
    if (i < keep.size() && keep[i]) {
      renumber[i] = static_cast<long>(kept.size());
      kept.push_back(all[i]);
    }
  // Nothing to leave out; or nothing but crystal-less experiments -- an
  // imported list, before indexing, which DIALS reads as it is.
  if (kept.size() == all.size() || kept.empty())
    return document;
  json::Object out = document.as_object();
  out["experiment"] = json::Value(std::move(kept));
  if (table != nullptr && table->has("id")) {
    Column &id = table->modify_int_column("id", "int", 1);
    for (std::int64_t &i : id.ints)
      i = i >= 0 && static_cast<std::size_t>(i) < renumber.size()
              ? renumber[static_cast<std::size_t>(i)]
              : -1;
    std::map<std::size_t, std::string> identifiers;
    for (const auto &[i, name] : table->identifiers)
      if (i < renumber.size() && renumber[i] >= 0)
        identifiers[static_cast<std::size_t>(renumber[i])] = name;
    table->identifiers = identifiers;
  }
  return json::Value(std::move(out));
}

} // namespace mxi
