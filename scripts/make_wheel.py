#!/usr/bin/env python3
"""A platform's wheel of mxi: the pure Python wheel (python/, built by
`python -m build --wheel python`) and the platform's mxi binary in one, so that
`pip install mxi` gives mxi, mxi-ui, mxi-view and mxeq at once, of one version
(Graeme, 8 November: everything in the wheel; the release archives mxi alone).

    scripts/make_wheel.py PURE.whl BINARY PLATFORM OUTDIR

PLATFORM as the release names its archives: linux-x86_64, linux-arm64,
macos-arm64, windows-x86_64, and the -cuda ones -- those a wheel of their own
name, mxi-cuda, mxi's with the CUDA build. The binary goes in the wheel's
scripts (NAME-VERSION.data/scripts/mxi), which pip installs into the
environment's bin (Scripts on Windows) as it is: no Python in front of it.
The tags: manylinux_2_28, the glibc the Linux builds are held to
(scripts/package_release.sh); macOS 13 on Apple silicon; win_amd64.
"""

import base64
import hashlib
import pathlib
import re
import sys
import zipfile

TAGS = {
    "linux-x86_64": "manylinux_2_28_x86_64",
    "linux-arm64": "manylinux_2_28_aarch64",
    "macos-arm64": "macosx_13_0_arm64",
    "windows-x86_64": "win_amd64",
}


def record_line(path, data):
    digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
    return f"{path},sha256={digest},{len(data)}"


def make(pure, binary, platform, outdir):
    cuda = platform.endswith("-cuda")
    base = platform[: -len("-cuda")] if cuda else platform
    if base not in TAGS:
        raise SystemExit(f"make_wheel: no wheel tag for the platform {platform}")
    tag = TAGS[base]
    m = re.match(r"(?P<name>[^-]+)-(?P<version>[^-]+)-py3-none-any\.whl$", pathlib.Path(pure).name)
    if not m:
        raise SystemExit(f"make_wheel: not a pure wheel: {pure}")
    old, version = m.group("name"), m.group("version")
    name = "mxi_cuda" if cuda else "mxi"  # the wheel's normalised name
    exe = "mxi.exe" if base.startswith("windows") else "mxi"
    entries = []  # (path, bytes, executable)
    with zipfile.ZipFile(pure) as z:
        for info in z.infolist():
            path, data = info.filename, z.read(info)
            if path.startswith(f"{old}-{version}.dist-info/"):
                rest = path[len(f"{old}-{version}.dist-info/"):]
                if rest == "RECORD":
                    continue
                path = f"{name}-{version}.dist-info/{rest}"
                if rest == "WHEEL":
                    text = data.decode()
                    text = re.sub(r"^Tag: .*$", f"Tag: py3-none-{tag}", text, flags=re.M)
                    data = text.encode()
                elif rest == "METADATA":
                    text = data.decode()
                    text = re.sub(r"^Name: .*$", f"Name: {'mxi-cuda' if cuda else 'mxi'}",
                                  text, count=1, flags=re.M)
                    data = text.encode()
            entries.append((path, data, False))
    entries.append((f"{name}-{version}.data/scripts/{exe}",
                    pathlib.Path(binary).read_bytes(), True))
    record = f"{name}-{version}.dist-info/RECORD"
    lines = [record_line(p, d) for p, d, _ in entries] + [f"{record},,"]
    entries.append((record, ("\n".join(lines) + "\n").encode(), False))
    out = pathlib.Path(outdir) / f"{name}-{version}-py3-none-{tag}.whl"
    out.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for path, data, executable in entries:
            info = zipfile.ZipInfo(path, date_time=(2026, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (0o100755 if executable else 0o100644) << 16
            z.writestr(info, data)
    return out


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    print(make(*sys.argv[1:]))
