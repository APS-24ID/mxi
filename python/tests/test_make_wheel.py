"""scripts/make_wheel.py: a platform's wheel, the Python and mxi in one
(Graeme, 8 November: everything in the wheel, the release archives mxi
alone)."""

import base64
import hashlib
import importlib.util
import os
import pathlib
import subprocess
import sys
import zipfile

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "make_wheel", ROOT / "scripts/make_wheel.py"
)
make_wheel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(make_wheel)


@pytest.fixture(scope="module")
def pure(tmp_path_factory):
    out = tmp_path_factory.mktemp("pure")
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "wheel",
            "--no-deps",
            "--no-build-isolation",
            "-q",
            "-w",
            str(out),
            str(ROOT / "python"),
        ],
        capture_output=True,
        text=True,
    )
    if run.returncode != 0:
        # Built without isolation, with the setuptools installed: the licence
        # as an SPDX expression (PEP 639) needs 77 or later -- older, and
        # every test here skipped quietly (26 November).
        pytest.skip(
            "no wheel built here -- python/pyproject.toml needs setuptools 77 "
            f"or later (pip install -U setuptools): {run.stderr[-300:]}"
        )
    (wheel,) = out.glob("mxi-*-py3-none-any.whl")
    return wheel


@pytest.fixture
def program(tmp_path):
    # A stand-in for mxi: what matters is where it lands and that it runs.
    p = tmp_path / "mxi"
    p.write_text('#!/bin/sh\necho stand-in mxi "$@"\n')
    os.chmod(p, 0o755)
    return p


def entries(wheel):
    with zipfile.ZipFile(wheel) as z:
        return {i.filename: (z.read(i), i.external_attr >> 16) for i in z.infolist()}


@pytest.mark.parametrize(
    "platform, name, tag, exe",
    [
        ("linux-x86_64", "mxi", "manylinux_2_28_x86_64", "mxi"),
        ("linux-arm64-cuda", "mxi_cuda", "manylinux_2_28_aarch64", "mxi"),
        ("macos-arm64", "mxi", "macosx_13_0_arm64", "mxi"),
        ("windows-x86_64", "mxi", "win_amd64", "mxi.exe"),
    ],
)
def test_a_platforms_wheel_holds_mxi_and_says_what_it_is(
    pure, program, tmp_path, platform, name, tag, exe
):
    version = pure.name.split("-")[1]
    wheel = make_wheel.make(pure, program, platform, tmp_path / "out")
    assert wheel.name == f"{name}-{version}-py3-none-{tag}.whl"
    e = entries(wheel)
    data, mode = e[f"{name}-{version}.data/scripts/{exe}"]
    assert data == program.read_bytes() and mode == 0o100755, "mxi, executable"
    info = f"{name}-{version}.dist-info"
    assert f"Tag: py3-none-{tag}" in e[f"{info}/WHEEL"][0].decode()
    metadata = e[f"{info}/METADATA"][0].decode()
    assert f"Name: {name.replace('_', '-')}\n" in metadata
    assert "mxi_ui/app.py" in e and "mxeq/__init__.py" in e, "the Python too"
    # Every file's hash and size as RECORD says; RECORD itself unhashed.
    for line in e[f"{info}/RECORD"][0].decode().strip().split("\n"):
        path, digest, size = line.rsplit(",", 2)
        if path == f"{info}/RECORD":
            assert digest == "" and size == ""
            continue
        body = e[path][0]
        want = base64.urlsafe_b64encode(hashlib.sha256(body).digest()).rstrip(b"=")
        assert digest == "sha256=" + want.decode() and int(size) == len(body), path


@pytest.mark.skipif(os.name == "nt", reason="the stand-in mxi is a shell script")
def test_installed_mxi_runs_from_the_environment_not_activated(pure, program, tmp_path):
    wheel = make_wheel.make(pure, program, "linux-x86_64", tmp_path / "out")
    venv = tmp_path / "venv"
    subprocess.run([sys.executable, "-m", "venv", str(venv)], check=True)
    # As a user's pip: not this suite's PYTHONPATH, whose src holds the
    # egg-info of the wheel's own build -- through which pip took mxi as
    # installed already and, quietly, installed nothing.
    clean = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    pip = subprocess.run(
        [str(venv / "bin/pip"), "install", "-q", "--no-deps", str(wheel)],
        capture_output=True,
        text=True,
        env=clean,
    )
    if pip.returncode != 0 and "is not a supported wheel" in pip.stderr:
        pytest.skip("not a manylinux_2_28 x86_64 system")
    assert pip.returncode == 0, pip.stderr
    run = subprocess.run(
        [str(venv / "bin/mxi"), "--version"],
        capture_output=True,
        text=True,
        env={"PATH": "/usr/bin:/bin"},
    )
    assert run.stdout.strip() == "stand-in mxi --version"
    for command in ("mxi-ui", "mxi-view", "mxeq"):
        assert (venv / "bin" / command).exists(), command


def test_the_tag_gives_the_wheel_its_version():
    # Graeme: a candidate's wheels say rc as its tag does (PEP 440's form).
    v = make_wheel.wheel_version
    assert v("v1.3.0", "1.3.0") == "1.3.0"
    assert v("v1.3.0-rc2", "1.3.0") == "1.3.0rc2"
    assert v("v1.3.0-beta1", "1.3.0") == "1.3.0b1"
    assert v("v1.3.0-alpha3", "1.3.0") == "1.3.0a3"
    for bad in ("v1.2.0-rc1", "v1.3.1", "nightly", "v1.3.0-rc"):
        with pytest.raises(SystemExit):
            v(bad, "1.3.0")


def test_a_candidates_wheel_says_so_throughout(pure, program, tmp_path):
    source = pure.name.split("-")[1]
    wheel = make_wheel.make(
        pure, program, "macos-arm64", tmp_path / "out", f"v{source}-rc2"
    )
    version = f"{source}rc2"
    assert wheel.name == f"mxi-{version}-py3-none-macosx_13_0_arm64.whl"
    e = entries(wheel)
    assert all(
        not p.startswith(f"mxi-{source}.") for p in e
    ), "no trace of the plain version"
    assert f"Version: {version}\n" in e[f"mxi-{version}.dist-info/METADATA"][0].decode()
    assert f"mxi-{version}.data/scripts/mxi" in e
    for line in e[f"mxi-{version}.dist-info/RECORD"][0].decode().strip().split("\n"):
        path, digest, size = line.rsplit(",", 2)
        if digest:
            body = e[path][0]
            want = base64.urlsafe_b64encode(hashlib.sha256(body).digest()).rstrip(b"=")
            assert digest == "sha256=" + want.decode() and int(size) == len(body), path


def _notices(directory, names):
    directory.mkdir()
    for n in names:
        (directory / n).write_text(f"the notice {n}\n")
    return directory


def _metadata_head(wheel):
    with zipfile.ZipFile(wheel) as z:
        meta = z.read(next(n for n in z.namelist() if n.endswith("METADATA"))).decode()
    return meta.split("\n\n")[0].splitlines()


def test_the_wheel_carries_its_builds_notices_and_says_their_licences(
    pure, program, tmp_path
):
    """--licences (Graeme, 26 November): the notices package_release.sh put in
    the build's archive go into its wheel too, and the licence expression is
    composed from them -- FFTW's GPL exactly when the build has FFTW."""
    names = [
        "mxi.txt",
        "gemmi.txt",
        "bitshuffle.txt",
        "lz4.txt",
        "mimalloc.txt",
        "hdf5.txt",
        "libaec.txt",
        "zlib.txt",
        "fftw.txt",
        "README",
    ]
    lic = _notices(tmp_path / "licenses", names)
    wheel = make_wheel.make(
        pure, program, "linux-x86_64", tmp_path / "out", licences=lic
    )
    head = _metadata_head(wheel)
    expression = next(l for l in head if l.startswith("License-Expression: "))
    assert expression == (
        "License-Expression: BSD-3-Clause AND MPL-2.0 AND MIT AND "
        "BSD-2-Clause AND LicenseRef-HDF5 AND Zlib AND GPL-2.0-or-later"
    )
    files = [l.split(": ", 1)[1] for l in head if l.startswith("License-File: ")]
    assert set(f"notices/{n}" for n in names) <= set(files)
    with zipfile.ZipFile(wheel) as z:
        inside = {n.split("licenses/", 1)[1] for n in z.namelist() if "/licenses/" in n}
        record = z.read(next(n for n in z.namelist() if n.endswith("RECORD"))).decode()
    assert set(f"notices/{n}" for n in names) <= inside
    assert all(f"licenses/notices/{n}" in record for n in names)
    licenses = pytest.importorskip("packaging.licenses")
    assert licenses.canonicalize_license_expression(expression.split(": ", 1)[1])


def test_no_gpl_without_fftw_and_an_unknown_notice_refuses(pure, program, tmp_path):
    lic = _notices(tmp_path / "plain", ["mxi.txt", "gemmi.txt", "zlib.txt", "README"])
    wheel = make_wheel.make(
        pure, program, "macos-arm64", tmp_path / "out", licences=lic
    )
    expression = next(
        l for l in _metadata_head(wheel) if l.startswith("License-Expression")
    )
    assert "GPL" not in expression and "MPL-2.0" in expression
    odd = _notices(tmp_path / "odd", ["mxi.txt", "newlib.txt"])
    with pytest.raises(SystemExit, match="newlib.txt"):
        make_wheel.make(pure, program, "linux-x86_64", tmp_path / "out2", licences=odd)


def test_every_notice_package_release_writes_is_named():
    """make_wheel refuses a notice it cannot name; so every notice
    scripts/package_release.sh can write into licenses/ must be in NOTICES --
    checked here, not first in a release's build."""
    import re

    script = (ROOT / "scripts/package_release.sh").read_text()
    written = set(re.findall(r'\$lic/([A-Za-z0-9_.-]+)(?=["\s])', script))
    loop = re.search(r"for package in ([a-z0-9 ]+); do", script)
    assert loop, "the vcpkg notices' loop"
    written |= {f"{p}.txt" for p in loop.group(1).split()}
    assert written, "the notices found"
    assert written <= set(make_wheel.NOTICES), written - set(make_wheel.NOTICES)
