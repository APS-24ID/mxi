"""scripts/install_cuda.py, against a manifest and archives served here in
NVIDIA's layout: components merged into one prefix, nvcc executable, lib64 a
link to lib; an archive whose SHA-256 is not the manifest's refused; a
platform the manifest lacks refused; optional components where present."""

import functools
import hashlib
import http.server
import importlib.util
import io
import json
import os
import pathlib
import tarfile
import threading

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]


def load():
    spec = importlib.util.spec_from_file_location(
        "install_cuda", ROOT / "scripts/install_cuda.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def archive(top, files):
    """A .tar.xz with one top directory, as NVIDIA's are: {path: (bytes, mode)}."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:xz") as t:
        for name, (data, mode) in files.items():
            info = tarfile.TarInfo(f"{top}/{name}")
            info.size, info.mode = len(data), mode
            t.addfile(info, io.BytesIO(data))
    return buf.getvalue()


@pytest.fixture
def served(tmp_path):
    root = tmp_path / "redist"
    components = {
        "cuda_nvcc": {
            "bin/nvcc": (b"#!/bin/sh\necho nvcc\n", 0o755),
            "nvvm/libnvvm.so": (b"x", 0o644),
        },
        "cuda_cudart": {
            "include/cuda_runtime.h": (b"// runtime\n", 0o644),
            "lib/libcudart_static.a": (b"!<arch>\n", 0o644),
        },
        "cuda_cccl": {"include/cuda/std/version": (b"// cccl\n", 0o644)},
        "cuda_crt": {"include/crt/host_config.h": (b"// crt\n", 0o644)},
    }
    manifest = {"release_label": "12.8.1"}
    for name, files in components.items():
        rel = f"{name}/linux-x86_64/{name}-linux-x86_64-12.8.90-archive.tar.xz"
        data = archive(f"{name}-linux-x86_64-12.8.90-archive", files)
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_bytes(data)
        manifest[name] = {
            "version": "12.8.90",
            "linux-x86_64": {
                "relative_path": rel,
                "sha256": hashlib.sha256(data).hexdigest(),
            },
        }
    (root / "redistrib_12.8.1.json").write_text(json.dumps(manifest))
    handler = functools.partial(
        http.server.SimpleHTTPRequestHandler, directory=str(root)
    )
    handler.log_message = lambda *a: None
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{server.server_port}/", root, manifest
    server.shutdown()


def test_the_components_are_merged_into_one_prefix(served, tmp_path):
    base, root, manifest = served
    prefix = tmp_path / "cuda"
    done = load().install("12.8.1", "linux-x86_64", str(prefix), base=base)
    assert done == ["cuda_nvcc", "cuda_cudart", "cuda_cccl", "cuda_crt"]  # crt: present
    assert os.access(prefix / "bin/nvcc", os.X_OK)
    assert (prefix / "include/cuda_runtime.h").exists() and (
        prefix / "include/cuda/std/version"
    ).exists()
    assert (prefix / "lib64").is_symlink() and (
        prefix / "lib64/libcudart_static.a"
    ).exists()


def test_an_archive_not_the_manifests_is_refused(served, tmp_path):
    base, root, manifest = served
    rel = manifest["cuda_cudart"]["linux-x86_64"]["relative_path"]
    (root / rel).write_bytes((root / rel).read_bytes() + b"tampered")
    with pytest.raises(SystemExit, match="SHA-256"):
        load().install("12.8.1", "linux-x86_64", str(tmp_path / "cuda"), base=base)


def test_a_platform_the_manifest_lacks_is_refused(served, tmp_path):
    base, root, manifest = served
    with pytest.raises(SystemExit, match="no cuda_nvcc for linux-sbsa"):
        load().install("12.8.1", "linux-sbsa", str(tmp_path / "cuda"), base=base)
