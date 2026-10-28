#!/usr/bin/env python3
"""A version's section of RELEASE_NOTES.md, for the GitHub release's text, so
that the two always agree.

    scripts/release_notes.py 1.0.6 > notes.md
    scripts/release_notes.py v1.0.6-rc1 > notes.md   # a candidate: 1.0.6's

Its heading ("## 1.0.6 -- 26 October 2026") left out, the release having its
own title; macOS's quarantine note added, the binaries being unsigned.

A release must have its section: without one this fails, saying so, and the
workflow publishes nothing. A candidate (-rc, -alpha, -beta) is made before
its notes are final, so where its version has no section yet it takes the
"Unreleased" one, saying so -- v1.0.7-rc2 had built and tested everything,
then failed here for want of a 1.0.7 section (26 October).
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
    candidate = re.search(r"-(rc|alpha|beta)\d*$", version) is not None
    version = re.sub(r"-(rc|alpha|beta)\d*$", "", version)
    heading = re.compile(r"^## (\S+)", re.M)
    found = [(m.group(1), m.start(), m.end()) for m in heading.finditer(text)]

    def body(k):
        stop = found[k + 1][1] if k + 1 < len(found) else len(text)
        return text[found[k][2] : stop].split("\n", 1)[1].strip() + "\n"

    for k, (name, _, _) in enumerate(found):
        if name == version:
            return body(k)
    if candidate:
        for k, (name, _, _) in enumerate(found):
            if name == "Unreleased":
                return (
                    f"A candidate for {version}: the notes so far, as RELEASE_NOTES.md's "
                    "unreleased section has them.\n\n" + body(k)
                )
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
