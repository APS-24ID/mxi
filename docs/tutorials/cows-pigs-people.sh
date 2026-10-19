#!/bin/sh
# Every step of docs/tutorials/cows-pigs-people.md, its printout kept: the
# output the tutorial shows comes from a run of this.
#
#   docs/tutorials/cows-pigs-people.sh /path/to/data [where]
#
# DATA holds the 36 .nxs files (CIX*, PIX*, X*) from
# https://zenodo.org/records/23423865; WHERE (cows-pigs-people) is made, with
# cows/ and all/ in it, and output/ holding each step's printout as
# <part>-<step>.txt. Stops at the first step that fails. The mxi programs on
# PATH are the ones run; mxi_find --version is printed first, to say which.

set -eu

DATA=$(cd "${1:?usage: $0 DATA [WHERE]}" && pwd)
WHERE=${2:-cows-pigs-people}
mkdir -p "$WHERE/output"
WHERE=$(cd "$WHERE" && pwd)
OUT="$WHERE/output"

step() { # part, name, command...
  part=$1 name=$2
  shift 2
  echo "== $part: $*"
  "$@" > "$OUT/$part-$name.txt" 2>&1 || {
    echo "failed: see $OUT/$part-$name.txt" >&2
    exit 1
  }
}

mxi_find --version | tee "$OUT/version.txt"

# The cows alone.
mkdir -p "$WHERE/cows" && cd "$WHERE/cows"
step cows import mxi_import "$DATA"/CIX*.nxs
step cows find mxi_find imported.rflx
step cows index mxi_index strong.rflx --no-joint
mxi_show indexed.rflx | grep experiment > "$OUT/cows-cells.txt"
step cows refine mxi_refine indexed.rflx
step cows integrate mxi_integrate refined.rflx
step cows symmetry mxi_symmetry integrated.rflx
step cows scale mxi_scale symmetrized.rflx --d-min-auto
step cows correlation mxi_correlation scaled.rflx --html correlation.html

# All thirty-six.
mkdir -p "$WHERE/all" && cd "$WHERE/all"
step all import mxi_import "$DATA"/CIX*.nxs "$DATA"/PIX*.nxs "$DATA"/X*.nxs
step all find mxi_find imported.rflx
step all index mxi_index strong.rflx --no-joint
step all refine mxi_refine indexed.rflx
step all integrate mxi_integrate refined.rflx
step all symmetry mxi_symmetry integrated.rflx
step all scale mxi_scale symmetrized.rflx --d-min-auto
step all correlation mxi_correlation scaled.rflx --html correlation.html \
  --write-clusters
for c in 0 1 2; do
  [ -f "cluster_$c.rflx" ] || continue
  mkdir -p "$c" && cd "$c"
  step all "cluster$c-scale" mxi_scale "../cluster_$c.rflx" --d-min-auto
  cd ..
done
step all export mxi_export scaled.rflx

echo "done: $OUT"
