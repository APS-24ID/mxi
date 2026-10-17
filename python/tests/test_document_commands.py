"""Every mxi command in the documents takes the options it is shown with.

The documents drift as the programs change: docs/multi-sweep.md and the
review ran the chain on .expt/.refl pairs long after each program wrote one
.rflx, and docs/gpu.md showed `mxi_refine --beam --analytic` as "the command
anyone actually runs" after both had become defaults and stopped being
options (Graeme, October 2026). Each command -- in a fenced block, a four-space
block, or inline in prose -- has its options checked against the program's own
--help. Needs MXI_FIND, for the build directory every program is in.
"""

import os
import pathlib
import re
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
FIND = os.environ.get("MXI_FIND")
BUILD = pathlib.Path(FIND).parent if FIND else None

# Records of the past, whose commands were right when they were run.
HISTORY = {"docs/integration_history.md", "CLAUDE.md"}
# (document, option): shown deliberately, with the reason.
ALLOWED = {
    ("docs/backstop.md", "--flag-background"): "a proposed option, not yet written",
}

DOCUMENTS = sorted(
    str(p.relative_to(ROOT))
    for p in list(ROOT.glob("docs/*.md"))
    + [ROOT / "README.md", ROOT / "python/README.md"]
    if str(p.relative_to(ROOT)) not in HISTORY
)

COMMAND = re.compile(r"(?:\./build/|build/)?(mxi_[a-z]+)\b(.*)")


def commands(text):
    """(line number, program, the rest) for every command in a document."""
    in_fence = False
    for n, line in enumerate(text.split("\n"), 1):
        if line.startswith("```"):
            in_fence = not in_fence
            continue
        found = []
        if in_fence or re.match(r"    (?:\$ )?(?:\./build/)?mxi_[a-z]+ ", line):
            found.append(line.strip().removeprefix("$ "))
        found += re.findall(r"`((?:\./build/)?mxi_[a-z]+ [^`]+)`", line)
        for candidate in found:
            m = COMMAND.match(candidate)
            if m:
                yield n, m.group(1), m.group(2).split("#")[0]


_helps = {}


def help_of(program):
    if program not in _helps:
        exe = BUILD / program
        if not exe.exists():
            _helps[program] = None
        else:
            out = subprocess.run([str(exe), "--help"], capture_output=True, text=True)
            _helps[program] = out.stdout + out.stderr
    return _helps[program]


@pytest.mark.skipif(not FIND, reason="set MXI_FIND")
@pytest.mark.parametrize("document", DOCUMENTS)
def test_every_option_shown_is_an_option(document):
    wrong = []
    for n, program, rest in commands((ROOT / document).read_text()):
        text = help_of(program)
        if text is None:
            continue  # programs named but not built are test_documents.py's
        for option in re.findall(r"(?<!\S)(--?[a-zA-Z][\w-]*)", rest):
            if option in ("--help", "-h", "--version") or option in text:
                continue
            if (document, option) in ALLOWED:
                continue
            wrong.append(f"{document}:{n}: {program} {option}")
    assert not wrong, "options the programs do not take:\n" + "\n".join(wrong)
