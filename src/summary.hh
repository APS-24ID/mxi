// The processing's summary (docs/output.md): each program's figures and
// tables -- what a crystallographer reads off its output -- kept in the .rflx
// it writes, under /summary, and carried from file to file, so that a file
// says how it was made and mxi_show can tell the story of a run.
//
// The document is {"steps": [...]}, one entry a program in the order they ran,
// a program run twice there twice: {"program", "time", then its own}. A table
// is {"columns": [names in order], name: [values], ...}, column by column, as
// the file stores it.

#ifndef MXI_SUMMARY_HH
#define MXI_SUMMARY_HH

#include <string>
#include <vector>

#include "json.hh"

namespace mxi {

class Summary {
public:
  //: This program's step, `program` its name; the steps before it from
  //: `input`, a .rflx, if it has any -- none from a pair, which has no place
  //: for them.
  Summary(const std::string &program, const std::string &input);

  //: A figure of this step: a number, a string, a list.
  void set(const std::string &name, json::Value value);
  //: A table of this step, its columns in order, each as long as the others.
  void table(const std::string &name, const std::vector<std::string> &columns,
             const std::vector<std::vector<json::Value>> &values);

  //: Every step, this one last: what is written.
  json::Value document() const;

  //: Every step of a file's summary, `program` and `time` first.
  static std::vector<json::Value> steps(const json::Value &summary);

private:
  std::vector<json::Value> before_;
  json::Object own_;
};

//: Whether the developer's lines are printed (--verbose): to the terminal and
//: the log alike, as everything else is. Without it, the crystallography alone.
void set_verbose(bool on);
bool verbose();
//: printf, under --verbose only.
void detail(const char *format, ...)
#if defined(__GNUC__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;

} // namespace mxi

#endif // MXI_SUMMARY_HH
