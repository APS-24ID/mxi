# mxeq -- working notes

Read this before changing anything. It records constraints that were arrived at
rather than chosen, and the reasons are not recoverable from the code.

## What this is for

Checking that a GPU-accelerated MX pipeline produces output **equivalent** to
DIALS', boundary by boundary. DIALS is the acceptance oracle. `mx` is the
development oracle. Neither is imported here.

## Hard constraints

**No cctbx, no DIALS, no dxtbx imports. Ever.** The whole point is that this
runs in a container with no crystallographic software in it, against files
written by a DIALS build somewhere else. Reaching for `dials.array_family` to
read a `.refl` would be more convenient and would destroy the tool.

**No thresholds in version one.** Every check reports distributions. A
tolerance chosen before the distribution has been seen on real data is a
tolerance chosen by guessing, and it will then be tuned until it passes, which
is worse than having none. Thresholds go in a config file against the `--json`
output once `ins10_1` has been through it.

**Equivalence, not identity.** Identity is the right target for the threshold
kernel and is achievable there. It is not achievable at integration and
demanding it would report every difference in summation order as a failure.

**Exit status is not a verdict.** Zero means the comparison ran. If that ever
changes, the change must be explicit and documented here, because scripts will
have been written against the current meaning.

## Format assumptions -- STATUS

The `.refl` layout in `src/mxeq/refl.py` is taken from DIALS' msgpack adapter:
a two-element array tagged `dials::af::reflection_table`, then a map with
`nrows`, the identifiers, and `data` mapping column name to a
`[type_name, payload]` pair, payloads being raw little-endian blobs with
compound components adjacent.

**This round-trips against this package's own writer and has not been validated
against a `.refl` written by a real DIALS.** `mxeq inspect` walks an arbitrary
msgpack document without assuming any key name, and exists precisely for that
first contact. When a real file is available:

1. `mxeq inspect real.refl`
2. If it reads, the assumption held. Record that here and delete this block.
3. If it does not, the structure dump says what the difference is. It will be
   a key name or a payload encoding, and the fix is small.

The identifiers key is read under both `identifiers` and
`experiment_identifiers` because DIALS has spelled it more than one way.

## Things that were got wrong once

**The A matrix is the plain inverse of the real-space matrix, not the inverse
transpose.** The real-space matrix has a, b, c as its *rows*, so `M A = I` by
the definition of the reciprocal basis. The inverse transpose passes every test
on an orthogonal cell, because both are diagonal there. `test_readers.py` uses a
triclinic cell for this and nothing else.

**The candidate operator pool must be closed under transposition explicitly.**
It was originally claimed that the pool was closed already, because it is built
from lattice groups. The cubic group is -- its metric is the identity, so
`M^T = M^-1`. The hexagonal group is not, for the same reason in reverse. A
test caught it. Left unclosed, the reindexing search would have been
convention-dependent for hexagonal and rhombohedral lattices only.

**`fraction_matched` divides by the larger table, not the smaller.** Otherwise a
table containing one reflection that happens to be in the other scores 1.0.

## Design decisions worth not relitigating

**One module per boundary, not one parameterised function.** The join, the
metrics and the failure modes differ at every boundary. Sharing a body between
them would put the differences in flags.

**The reindexing operator is found from the data, never assumed, and never
applied silently.** It is found at the `indexed` boundary and passed explicitly
to later ones with `--operator`. A later boundary given no operator when one is
needed reports a low match fraction and says so; it does not go looking.

**Duplicate keys are resolved by frame order.** Never by intensity: that would
pair whichever two happened to agree.

**CC-half splits randomly with a fixed seed, not on observation parity.**
Observation order is pipeline-dependent, and splitting on it would make CC-half
compare two pipelines' orderings rather than their data.

**Resolution comes from the `.expt` in preference to the `d` column.** `d` is a
pipeline output; binning by it bins the comparison by one of the things being
compared.

## Testing

`tests/fixtures.py` generates synthetic data with *planted* answers -- a known
reindexing operator, a known intensity ratio, a known misorientation. A test
that asserts against a number nobody chose is a test that passes for reasons
nobody knows.

Fixtures are imported by name via `tests/conftest.py`, not shipped in the
package. They generate test data and have no business being installed.

## Style

`black` formatting. Four-space indent, double quotes. No `from x import *`.
Comments explain why, not what; a comment restating the line above it is
deleted, not improved.
