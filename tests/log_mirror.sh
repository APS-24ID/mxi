#!/bin/sh
# Every program mirrors its output to mxi_<program>.log where it runs, as
# DIALS writes dials.<program>.log. Checked through the real binaries, since
# the mirror exists only at the level of a process.
#
#   log_mirror.sh MXI_BACKGROUND MXI_INDEX MXI_FIND
set -u
# mxi itself, its programs its subcommands: the links to it by the programs'
# names are not made on Windows (CMakeLists.txt).
mxi=$1
background() { "$mxi" background "$@"; }
index() { "$mxi" index "$@"; }
find_spots() { "$mxi" find "$@"; }
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
cd "$dir" || exit 1
fail() { echo "FAIL: $*"; exit 1; }

# The log is standard output, byte for byte.
echo "3 4 5 3 4 5 4 3 5 4" | background > out 2> err || fail "mxi_background failed"
# Carriage returns aside: on Windows standard output is in text mode, each
# newline arriving in the file as CR LF; the log is written as printed.
tr -d '\r' < out | cmp -s - mxi_background.log || fail "mxi_background.log is not its standard output"

# An error is on standard error and in the log, and not on standard output.
index /nonexistent.expt /nonexistent.refl > out 2> err && fail "mxi_index should fail"
grep -q "cannot open" err || fail "the error is not on standard error"
grep -q "cannot open" mxi_index.log || fail "the error is not in mxi_index.log"
grep -q "cannot open" out && fail "the error is on standard output"

# Asking for help writes no log, and leaves the last run's alone.
cp mxi_index.log saved
index --help > /dev/null 2>&1
cmp -s saved mxi_index.log || fail "--help overwrote mxi_index.log"
find_spots --help > /dev/null 2>&1
[ -f mxi_find.log ] && fail "--help created mxi_find.log"

# The spot finder: its complaint, as every program's, mirrored to its log.
find_spots --no-such-flag > out 2> err
[ -s err ] || fail "mxi_find reported nothing for a bad flag"
grep -q "unknown option" mxi_find.log || fail "mxi_find's complaint is not in mxi_find.log"

echo "log mirror: every check passed"
