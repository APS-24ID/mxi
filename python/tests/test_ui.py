"""mxi ui without a window: the dependency check and its install hints, the
directory's state, the runner, the steps' commands, and mxi ui through mxi."""

import os
import subprocess
import sys
import time

import pytest

from mxi_ui import deps, state
from mxi_ui.runner import Runner, as_text
from mxi_ui.steps import STEPS, by_key

BUILD = os.path.dirname(os.environ.get("MXI_FIND", "")) or os.environ.get(
    "MXI_BUILD", ""
)
MXI = os.path.join(BUILD, "mxi")


@pytest.mark.parametrize(
    "system, family, expected",
    [
        ("Linux", "ubuntu", "sudo apt install python3-tk"),
        ("Linux", "debian", "sudo apt install python3-tk"),
        ("Linux", "rhel", "sudo dnf install python3-tkinter"),
        ("Linux", "fedora", "sudo dnf install python3-tkinter"),
        ("Linux", "arch", "sudo pacman -S tk"),
        ("Darwin", "", "brew install python-tk@3.12"),
    ],
)
def test_what_to_install_for_tkinter_is_said_for_the_system(system, family, expected):
    assert expected in deps.tkinter_hint(system=system, family=family, version=(3, 12))


def test_the_family_is_read_from_os_release(tmp_path):
    release = tmp_path / "os-release"
    release.write_text('NAME="Rocky Linux"\nID="rocky"\nID_LIKE="rhel centos fedora"\n')
    assert deps._linux_family(str(release)) == "fedora"


def test_a_missing_package_is_named_with_how_to_install_it(monkeypatch):
    real = deps.importlib.util.find_spec
    monkeypatch.setattr(
        deps.importlib.util,
        "find_spec",
        lambda name: None if name == "h5py" else real(name),
    )
    found = deps.problems(need_display=False)
    assert any("h5py missing" in p and "pip install" in p for p in found)


def test_a_directory_shows_its_steps_done(tmp_path):
    assert state.next_step(str(tmp_path)) == "import"
    for name in ("imported.rflx", "strong.rflx"):
        (tmp_path / name).write_text("")
    done = state.scan(str(tmp_path))
    assert done["import"] and done["find"] and not done["index"]
    assert state.next_step(str(tmp_path)) == "index"


def test_each_step_reads_the_one_before():
    for before, after in zip(STEPS, STEPS[1:]):
        assert after.inputs == (before.output,), after.key
    assert by_key("index").command(extra=["-j", "4"]) == [
        "mxi",
        "index",
        "strong.rflx",
        "-j",
        "4",
    ]
    assert as_text(["mxi", "import", "a b.nxs"]) == "mxi import 'a b.nxs'"


def wait(runner, seconds=20):
    end = time.time() + seconds
    lines = []
    while time.time() < end:
        lines += runner.lines()
        code = runner.finished()
        if code is not None:
            return code, lines + runner.lines()
        time.sleep(0.02)
    raise AssertionError("the command did not finish")


def test_the_runner_streams_a_commands_output_and_its_exit_code(tmp_path):
    r = Runner()
    r.start(
        [
            sys.executable,
            "-c",
            "import sys\nfor i in range(3): print(i, flush=True)\nsys.exit(3)",
        ],
        cwd=str(tmp_path),
    )
    code, lines = wait(r)
    assert code == 3 and lines == ["0", "1", "2"]


def test_a_running_command_can_be_cancelled(tmp_path):
    r = Runner()
    r.start(
        [
            sys.executable,
            "-c",
            "import time\nprint('started', flush=True)\ntime.sleep(60)",
        ]
    )
    end = time.time() + 10
    while time.time() < end and "started" not in r.lines():
        time.sleep(0.02)
    assert r.running()
    r.cancel()
    code, _ = wait(r)
    assert code != 0 and not r.running()
    with pytest.raises(RuntimeError):
        r2 = Runner()
        r2.start([sys.executable, "-c", "import time; time.sleep(5)"])
        try:
            r2.start([sys.executable, "-c", "pass"])
        finally:
            r2.cancel()


@pytest.mark.skipif(not os.path.exists(MXI), reason="mxi is not built")
def test_mxi_ui_runs_through_mxi(tmp_path):
    run = subprocess.run([MXI, "ui", "--version"], capture_output=True, text=True)
    assert run.returncode == 0 and run.stdout.startswith("mxi ui "), run.stderr
    check = subprocess.run([MXI, "ui", "--check"], capture_output=True, text=True)
    assert check.returncode in (0, 1)
    assert "ready" in check.stdout or "mxi ui:" in check.stdout
    no_display = {
        k: v for k, v in os.environ.items() if k not in ("DISPLAY", "WAYLAND_DISPLAY")
    }
    if sys.platform.startswith("linux"):
        # The first problem here: no display, or, where this Python has no
        # tkinter (a release build's container), that -- said either way.
        expected = (deps.problems(need_display=False) or ["no display"])[0][:20]
        run = subprocess.run(
            [MXI, "ui"], capture_output=True, text=True, env=no_display
        )
        assert run.returncode == 1 and expected in run.stderr, run.stderr
