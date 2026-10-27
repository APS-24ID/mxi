#!/usr/bin/env python3
"""A version's section of RELEASE_NOTES.md, for the GitHub release's text, so
that the two always agree.

    scripts/release_notes.py 1.0.6 > notes.md
    scripts/release_notes.py v1.0.6-rc1 > notes.md   # a candidate: 1.0.6's

Its heading ("## 1.0.6 -- 26 October 2026") left out, the release having its
own title; macOS's quarantine note added, the binaries being unsigned. Fails,
saying so, where the notes have no section for the version.
"""

import pathlib
import re
import sys

MACOS = """
**macOS:** the binaries are not signed, so macOS blocks them once downloaded.
Clear the quarantine flag once, after unpacking:

    xattr -dr com.apple.quarantine mxi-<version>-macos-arm64
"""


def section(text: str, version: str) -> str:
    version = re.sub(r"^v", "", version)
    version = re.sub(r"-(rc|alpha|beta)\d*$", "", version)
    heading = re.compile(r"^## (\S+)", re.M)
    found = [(m.group(1), m.start(), m.end()) for m in heading.finditer(text)]
    for k, (name, start, end) in enumerate(found):
        if name == version:
            stop = found[k + 1][1] if k + 1 < len(found) else len(text)
            return text[end : stop].split("\n", 1)[1].strip() + "\n"
    raise SystemExit(f"RELEASE_NOTES.md has no section for {version}")


def main(argv):
    if len(argv) != 2:
        raise SystemExit("usage: release_notes.py VERSION")
    root = pathlib.Path(__file__).resolve().parent.parent
    notes = section((root / "RELEASE_NOTES.md").read_text(), argv[1])
    version = re.sub(r"^v", "", argv[1])
    sys.stdout.write(notes + MACOS.replace("<version>", version))


if __name__ == "__main__":
    main(sys.argv)
