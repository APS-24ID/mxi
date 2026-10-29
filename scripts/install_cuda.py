#!/usr/bin/env python3
"""The CUDA toolkit's parts the build needs, from NVIDIA's redistributable
archives, not a package repository.

    scripts/install_cuda.py 12.8.1 linux-x86_64 /usr/local/cuda
    scripts/install_cuda.py 12.8.1 linux-sbsa /usr/local/cuda   # ARM servers

NVIDIA's RHEL 8 repository served an index whose files it did not have --
"all mirrors were already tried without success", the same three files for
hours -- and no retry could get past it (v1.1.2-rc1 and rc2). The
redistributable archives are one per component and platform, the same for
every Linux, listed with their SHA-256 in a manifest for each CUDA release
(redistrib_<version>.json): no index to go stale. mxi's CUDA needs the
compiler and the runtime alone -- cuda_runtime.h, the runtime linked
statically -- with the headers the runtime may draw on: cuda_nvcc,
cuda_cudart, cuda_cccl, and cuda_crt and libnvvm where a release splits them
out of the compiler.

Each archive is checked against the manifest's SHA-256 and unpacked into the
prefix, its top directory dropped, the components merged; lib64 is made a link
to lib, where nvcc looks. Downloads are tried five times. Python's standard
library alone: it runs in the bare manylinux image.
"""

import hashlib
import json
import os
import sys
import tarfile
import tempfile
import time
import urllib.request

BASE = "https://developer.download.nvidia.com/compute/cuda/redist/"
REQUIRED = ["cuda_nvcc", "cuda_cudart", "cuda_cccl"]
IF_PRESENT = ["cuda_crt", "libnvvm"]


def fetch(url: str, attempts: int = 5) -> bytes:
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001 -- any failure is tried again
            if attempt == attempts:
                raise SystemExit(f"could not download {url}: {e}")
            wait = 15 * attempt
            print(f"  {url}: {e}; trying again in {wait} s", file=sys.stderr)
            time.sleep(wait)
    raise AssertionError


def unpack(archive: bytes, prefix: str) -> None:
    """The archive's contents into prefix, its one top directory dropped."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "a.tar.xz")
        with open(path, "wb") as f:
            f.write(archive)
        with tarfile.open(path) as t:
            members = t.getmembers()
            top = {m.name.split("/", 1)[0] for m in members}
            if len(top) != 1:
                raise SystemExit(f"an archive with more than one top directory: {top}")
            for m in members:
                parts = m.name.split("/", 1)
                if len(parts) < 2 or not parts[1]:
                    continue
                if m.name.startswith("/") or ".." in m.name.split("/"):
                    raise SystemExit(f"an unsafe path in the archive: {m.name}")
                m.name = parts[1]
                if m.islnk():
                    m.linkname = m.linkname.split("/", 1)[1]
                t.extract(m, prefix)


def install(version: str, platform: str, prefix: str, base: str = BASE) -> list:
    manifest = json.loads(fetch(f"{base}redistrib_{version}.json"))
    wanted = REQUIRED + [c for c in IF_PRESENT if c in manifest]
    os.makedirs(prefix, exist_ok=True)
    done = []
    for component in wanted:
        entry = manifest.get(component, {}).get(platform)
        if entry is None:
            raise SystemExit(f"CUDA {version} has no {component} for {platform}")
        url = base + entry["relative_path"]
        print(f"{component} {manifest[component].get('version', '')}: {url}")
        data = fetch(url)
        digest = hashlib.sha256(data).hexdigest()
        if digest != entry["sha256"]:
            raise SystemExit(
                f"{component}: SHA-256 {digest}, the manifest says {entry['sha256']}"
            )
        unpack(data, prefix)
        done.append(component)
    lib, lib64 = os.path.join(prefix, "lib"), os.path.join(prefix, "lib64")
    if os.path.isdir(lib) and not os.path.exists(lib64):
        os.symlink("lib", lib64)
    nvcc = os.path.join(prefix, "bin", "nvcc")
    if not os.access(nvcc, os.X_OK):
        raise SystemExit(f"no nvcc at {nvcc} after installing {done}")
    return done


def main(argv):
    if len(argv) != 4:
        raise SystemExit("usage: install_cuda.py VERSION PLATFORM PREFIX")
    version, platform, prefix = argv[1:]
    done = install(version, platform, prefix)
    print(f"CUDA {version} for {platform} in {prefix}: {', '.join(done)}")


if __name__ == "__main__":
    main(sys.argv)
