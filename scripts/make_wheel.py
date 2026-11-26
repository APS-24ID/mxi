#!/usr/bin/env python3
"""A platform's wheel of mxi: the pure Python wheel (python/, built by
`python -m build --wheel python`) and the platform's mxi binary in one, so that
`pip install mxi` gives mxi, mxi-ui, mxi-view and mxeq at once, of one version
(Graeme, 8 November: everything in the wheel; the release archives mxi alone).

    scripts/make_wheel.py PURE.whl BINARY PLATFORM OUTDIR [TAG] [--licences DIR]

PLATFORM as the release names its archives: linux-x86_64, linux-arm64,
macos-arm64, windows-x86_64, and the -cuda ones -- those a wheel of their own
name, mxi-cuda, mxi's with the CUDA build. The binary goes in the wheel's
scripts (NAME-VERSION.data/scripts/mxi), which pip installs into the
environment's bin (Scripts on Windows) as it is: no Python in front of it.
The tags: manylinux_2_28, the glibc the Linux builds are held to
(scripts/package_release.sh); macOS 13 on Apple silicon; win_amd64.

TAG, the release's git tag, gives the wheel its version as Python writes it
(PEP 440): v1.3.0 is 1.3.0, v1.3.0-rc2 1.3.0rc2, -beta1 b1, -alpha1 a1 -- so
that a candidate's wheel says it is one, as its tag does (Graeme, 8
November). Its release number must be the source's: a tag that is not
refuses, rather than a wheel saying a version it is not.

--licences DIR, the build's notices as scripts/package_release.sh gathers them
for its archive (mxi-VERSION-PLATFORM/licenses): each goes into the wheel too,
NAME-VERSION.dist-info/licenses/notices/, a License-File; and the wheel's
License-Expression is composed from the notices there, each named in NOTICES --
so that FFTW's GPL is declared exactly when the build has FFTW, NVIDIA's terms
for a CUDA build, mimalloc's for Linux (Graeme, 26 November: before PyPI, the
wheels carrying their build's notices as its archive does). A notice not in
NOTICES refuses: a library new to the build needs its licence named here.
Without --licences, the pure wheel's own expression stands (BSD-3-Clause, the
source's), right only for a build without FFTW -- the releases give it.
"""

import base64
import hashlib
import pathlib
import re
import sys
import zipfile

# Each notice package_release.sh may put in licenses/, and its SPDX licence;
# None for the README. HDF5's licence and NVIDIA's EULA have no SPDX
# identifier, so LicenseRef- (PEP 639).
NOTICES = {
    "mxi.txt": "BSD-3-Clause",
    "gemmi.txt": "MPL-2.0",
    "bitshuffle.txt": "MIT",
    "lz4.txt": "BSD-2-Clause",
    "mimalloc.txt": "MIT",
    "hdf5.txt": "LicenseRef-HDF5",
    "libaec.txt": "BSD-2-Clause",
    "zlib.txt": "Zlib",
    "fftw.txt": "GPL-2.0-or-later",
    "nvidia-cuda-runtime.txt": "LicenseRef-NVIDIA-CUDA-EULA",
    "README": None,
}


def notices(directory):
    """The notice files in DIRECTORY, sorted, and the licence expression
    they make: each SPDX identifier once, in NOTICES' order."""
    files = sorted(p for p in pathlib.Path(directory).iterdir() if p.is_file())
    if not files:
        raise SystemExit(f"make_wheel: no notices in {directory}")
    unknown = [f.name for f in files if f.name not in NOTICES]
    if unknown:
        raise SystemExit(f"make_wheel: notices not named in NOTICES: {', '.join(unknown)}")
    present = {f.name for f in files}
    ids = []
    for name, spdx in NOTICES.items():
        if name in present and spdx and spdx not in ids:
            ids.append(spdx)
    return files, " AND ".join(ids)


TAGS = {
    "linux-x86_64": "manylinux_2_28_x86_64",
    "linux-arm64": "manylinux_2_28_aarch64",
    "macos-arm64": "macosx_13_0_arm64",
    "windows-x86_64": "win_amd64",
}


def wheel_version(tag, source):
    """The version a tag gives the wheel, the source's being `source`."""
    m = re.fullmatch(r"v?(\d+\.\d+\.\d+)(?:-(rc|beta|alpha)(\d+))?", tag)
    if not m:
        raise SystemExit(f"make_wheel: a tag not of a release: {tag}")
    release, kind, n = m.groups()
    if release != source:
        raise SystemExit(f"make_wheel: the tag {tag} is not of the source's version {source}")
    return release + ({"rc": "rc", "beta": "b", "alpha": "a"}[kind] + n if kind else "")


def record_line(path, data):
    digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
    return f"{path},sha256={digest},{len(data)}"


def make(pure, binary, platform, outdir, tag=None, licences=None):
    cuda = platform.endswith("-cuda")
    base = platform[: -len("-cuda")] if cuda else platform
    if base not in TAGS:
        raise SystemExit(f"make_wheel: no wheel tag for the platform {platform}")
    platform_tag = TAGS[base]
    m = re.match(r"(?P<name>[^-]+)-(?P<version>[^-]+)-py3-none-any\.whl$", pathlib.Path(pure).name)
    if not m:
        raise SystemExit(f"make_wheel: not a pure wheel: {pure}")
    old, source = m.group("name"), m.group("version")
    version = wheel_version(tag, source) if tag else source
    name = "mxi_cuda" if cuda else "mxi"  # the wheel's normalised name
    exe = "mxi.exe" if base.startswith("windows") else "mxi"
    notice_files, expression = notices(licences) if licences else ([], None)
    entries = []  # (path, bytes, executable)
    with zipfile.ZipFile(pure) as z:
        for info in z.infolist():
            path, data = info.filename, z.read(info)
            if path.startswith(f"{old}-{source}.dist-info/"):
                rest = path[len(f"{old}-{source}.dist-info/"):]
                if rest == "RECORD":
                    continue
                path = f"{name}-{version}.dist-info/{rest}"
                if rest == "WHEEL":
                    text = data.decode()
                    text = re.sub(r"^Tag: .*$", f"Tag: py3-none-{platform_tag}", text, flags=re.M)
                    data = text.encode()
                elif rest == "METADATA":
                    text = data.decode()
                    text = re.sub(r"^Name: .*$", f"Name: {'mxi-cuda' if cuda else 'mxi'}",
                                  text, count=1, flags=re.M)
                    text = re.sub(r"^Version: .*$", f"Version: {version}", text, count=1,
                                  flags=re.M)
                    if expression:
                        text, n = re.subn(r"^License-Expression: .*$",
                                          f"License-Expression: {expression}", text,
                                          count=1, flags=re.M)
                        if n != 1:
                            raise SystemExit("make_wheel: the pure wheel has no License-Expression")
                        # The notices' License-File lines after the last of the
                        # pure wheel's own, in the headers before the README.
                        head, sep, body = text.partition("\n\n")
                        lines = head.split("\n")
                        last = max(i for i, l in enumerate(lines) if l.startswith("License-File: "))
                        lines[last + 1:last + 1] = [f"License-File: notices/{f.name}"
                                                    for f in notice_files]
                        text = "\n".join(lines) + sep + body
                    data = text.encode()
            entries.append((path, data, False))
    for f in notice_files:
        entries.append((f"{name}-{version}.dist-info/licenses/notices/{f.name}",
                        f.read_bytes(), False))
    entries.append((f"{name}-{version}.data/scripts/{exe}",
                    pathlib.Path(binary).read_bytes(), True))
    record = f"{name}-{version}.dist-info/RECORD"
    lines = [record_line(p, d) for p, d, _ in entries] + [f"{record},,"]
    entries.append((record, ("\n".join(lines) + "\n").encode(), False))
    out = pathlib.Path(outdir) / f"{name}-{version}-py3-none-{platform_tag}.whl"
    out.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for path, data, executable in entries:
            info = zipfile.ZipInfo(path, date_time=(2026, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (0o100755 if executable else 0o100644) << 16
            z.writestr(info, data)
    return out


if __name__ == "__main__":
    args = sys.argv[1:]
    licences = None
    if "--licences" in args:
        at = args.index("--licences")
        if at + 1 >= len(args):
            raise SystemExit(__doc__)
        licences = args[at + 1]
        del args[at:at + 2]
    if len(args) not in (4, 5):
        raise SystemExit(__doc__)
    print(make(*args, licences=licences))
