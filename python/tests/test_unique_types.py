"""Every type mxi defines at namespace level has one name of its own.

Two classes of one name in one namespace break C++'s one-definition rule: the
program links without complaint, but the linker keeps one copy of each inline
member -- a destructor, a constructor with default member initialisers -- for
both, and the other class's objects are made or destroyed with it. So it was
(26 November): mxi merge's MergedSet beside mxi correlation's, and
test_correlation's MergedSet destroyed with merge's destructor -- a stack
overflow that crashed the C++ tests on Linux and Windows, not on macOS, whose
linker kept the other. And forward.hh's Moments beside triage.hh's, latent.
Types in an anonymous namespace are their file's own, and exempt.
"""

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]

DEFINITION = re.compile(r"^(?:struct|class)\s+([A-Za-z_]\w*)\s*(?:final\s*)?(?:\{|:(?!:))")
ANONYMOUS_OPEN = re.compile(r"^namespace\s*\{")
ANONYMOUS_CLOSE = re.compile(r"^\}\s*//\s*namespace\s*$")


def definitions():
    found = {}
    for path in sorted(
        list((ROOT / "src").glob("*.hh"))
        + list((ROOT / "src").glob("*.cc"))
        + list((ROOT / "apps").glob("*.cc"))
    ):
        anonymous = 0
        for number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            if ANONYMOUS_OPEN.match(line):
                anonymous += 1
            elif ANONYMOUS_CLOSE.match(line) and anonymous:
                anonymous -= 1
            elif not anonymous:
                m = DEFINITION.match(line)
                if m:
                    found.setdefault(m.group(1), []).append(
                        f"{path.relative_to(ROOT)}:{number}"
                    )
    return found


def test_no_type_is_defined_twice_at_namespace_level():
    found = definitions()
    assert len(found) > 50, "the scan found mxi's types"
    twice = {name: places for name, places in found.items() if len(places) > 1}
    assert not twice, twice
