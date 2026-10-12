// The two views of an experiment list: see views.hh.

#include "views.hh"

#include "shoebox.hh"

#include <algorithm>
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

Table lattice_rows(const Table &table, const Lattices &lattices) {
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
  Table out = select_rows_with_shoeboxes(table, rows);
  Column &ids = out.modify_int_column("id", "int", 1);
  for (std::int64_t &i : ids.ints)
    i = renumber[static_cast<std::size_t>(i)];
  out.identifiers.clear();
  for (const auto &[i, name] : table.identifiers)
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

Table rows_again(const Table &lattice_rows, const Table &crystal_less_rows,
                 const Lattices &lattices, const ExperimentList &whole) {
  if (!lattices.any_crystal_less())
    return lattice_rows;
  Table back = lattice_rows;
  Column &ids = back.modify_int_column("id", "int", 1);
  for (std::int64_t &i : ids.ints)
    i = static_cast<std::int64_t>(
        lattices.position.at(static_cast<std::size_t>(i)));
  // The others with exactly the lattices' columns: what the lattices' rows
  // have and theirs do not, as zeros -- an unindexed spot has no prediction --
  // and what theirs have and the lattices' do not, dropped: the lattices'
  // are what a program wrote. Their shoeboxes, if any, kept as they were.
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

} // namespace mxi
