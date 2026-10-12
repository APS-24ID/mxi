// The processing's summary in the .rflx (docs/output.md): written, read back,
// and carried from one program's file to the next.

#include <cstdio>
#include <filesystem>

#include "../src/processing.hh"
#include "../src/rflx.hh"
#include "check.hh"

namespace mxi {

TEST(each_program_adds_its_step_to_the_summary_it_was_given) {
  const std::string dir =
      (std::filesystem::temp_directory_path() / "mxi_summary_test").string();
  std::filesystem::create_directories(dir);
  const std::string first = dir + "/indexed.rflx",
                    second = dir + "/refined.rflx";
  const json::Value experiments = json::parse(R"({"experiment": []})");

  // The first program: a figure and a table of its own.
  Summary a("mxi_index", "");
  a.set("indexed_fraction", json::Value(0.982));
  a.table("cycles", {"cycle", "rmsd_x"},
          {{json::Value(1L), json::Value(2L)},
           {json::Value(0.186), json::Value(0.185)}});
  const json::Value da = a.document();
  rflx::write(first, &experiments, nullptr, "mxi_index", false, &da);
  const auto back = rflx::read_summary(first);
  check::is_true(back.has_value(), "the summary read back");
  const auto steps_a = Summary::steps(*back);
  check::equal(static_cast<long long>(steps_a.size()), 1, "one step");
  check::is_true(steps_a[0]["program"].as_string() == "mxi_index",
                 "named for its program");
  const json::Value &cycles = steps_a[0]["cycles"];
  check::is_true(cycles["columns"].as_array().size() == 2 &&
                     cycles["rmsd_x"].as_array()[1].as_number() == 0.185,
                 "the table as it was");

  // The next program, given that file: its step after the first.
  Summary b("mxi_refine", first);
  b.set("reflections", json::Value(11743L));
  const json::Value db = b.document();
  rflx::write(second, &experiments, nullptr, "mxi_refine", false, &db);
  const auto steps_b = Summary::steps(*rflx::read_summary(second));
  check::equal(static_cast<long long>(steps_b.size()), 2, "both steps");
  check::is_true(steps_b[0]["program"].as_string() == "mxi_index" &&
                     steps_b[1]["program"].as_string() == "mxi_refine",
                 "in the order they ran");

  // A file with no summary: a fresh start, not a failure.
  rflx::write(first, &experiments, nullptr, "older");
  check::is_true(!rflx::read_summary(first).has_value(), "none written");
  check::equal(
      static_cast<long long>(
          Summary::steps(Summary("mxi_refine", first).document()).size()),
      1, "this step alone");
  std::filesystem::remove_all(dir);
}

} // namespace mxi
