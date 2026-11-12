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
        pytest.skip(f"no wheel built here: {run.stderr[-300:]}")
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
