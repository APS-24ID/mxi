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


from mxi_ui.model import build, problems  # noqa: E402
from mxi_ui.steps import TOOLS  # noqa: E402


def test_a_form_builds_its_command():
    scale = by_key("scale")
    argv = build(
        scale,
        ["symmetrized.rflx"],
        {
            "--anomalous": True,
            "--d-min": " 1.8 ",
            "--space-group": "",
            "--absorption-level": "medium",
        },
        output="scaled.rflx",  # the default: not said
        threads="4",
        extra="--shells 10",
    )
    assert argv == [
        "mxi",
        "scale",
        "symmetrized.rflx",
        "--anomalous",
        "--d-min",
        "1.8",
        "--absorption-level",
        "medium",
        "-j",
        "4",
        "--shells",
        "10",
    ]
    assert build(
        scale, ["a.rflx"], {"--anomalous": False}, output="mine.rflx", threads="0"
    ) == [
        "mxi",
        "scale",
        "a.rflx",
        "-o",
        "mine.rflx",
    ], "an unticked box, threads 0 and the default output are all left out"
    verify = by_key("verify")
    assert build(verify, ["indexed.rflx"], {}, output="anything") == [
        "mxi",
        "verify",
        "indexed.rflx",
    ]


def test_a_forms_wrong_values_are_said():
    found = problems(
        by_key("scale"), {"--d-min": "two", "--absorption-level": "huge"}, threads="all"
    )
    assert len(found) == 3
    assert any("Resolution limit" in p for p in found) and any(
        "low, medium, high" in p for p in found
    )
    assert (
        problems(
            by_key("index"), {"--d-min": "1.5", "--max-lattices": "2"}, threads="8"
        )
        == []
    )


def takes(help_text, flag):
    """Whether a command's --help offers the flag: as a whole word on one of
    its option lines, those that begin with a dash."""
    import re

    word = re.compile(rf"(?<![\w-]){re.escape(flag)}(?![\w-])")
    return any(
        line.lstrip().startswith("-") and word.search(line)
        for line in help_text.splitlines()
    )


def test_takes_tells_an_option_from_one_that_is_not():
    text = "  -o, --output PATH  where\n  --d-min D  limit\n  --sigma-b B --sigma-m M   profile\n"
    assert takes(text, "-o") and takes(text, "--d-min") and takes(text, "--sigma-m")
    assert (
        not takes(text, "--d")
        and not takes(text, "--d-max")
        and not takes(text, "--output-expt")
    )


@pytest.mark.skipif(not os.path.exists(MXI), reason="mxi is not built")
@pytest.mark.parametrize("step", STEPS + TOOLS, ids=lambda s: s.key)
def test_every_option_on_a_form_is_one_the_command_takes(step):
    # The forms cannot drift from mxi: each option offered, and its output
    # option, is on one of the command's own --help option lines.
    helped = subprocess.run([MXI, step.key, "--help"], capture_output=True, text=True)
    text = helped.stdout + helped.stderr
    offered = [o.flag for o in step.options] + (
        [step.output_flag] if step.output_flag else []
    )
    missing = [flag for flag in offered if not takes(text, flag)]
    assert not missing, f"mxi {step.key} does not take {missing}"


@pytest.mark.skipif(not os.path.exists(MXI), reason="mxi is not built")
def test_a_summary_with_compressed_arrays_reads(tmp_path):
    # 300 images' spots per image is compressed (bitshuffle-LZ4, as every
    # array of 256 or more), as Graeme's 3600 were: read with the filters.
    import json

    import h5py

    from mxi_ui import results

    subprocess.run(
        [sys.executable, "-m", "mxeq.fixtures.nxmx", str(tmp_path), "300"],
        check=True,
        capture_output=True,
    )
    from mxeq.fixtures import nxmx

    nxmx.write_expt(str(tmp_path / "one.expt"), 300)
    doc = json.load(open(tmp_path / "one.expt"))
    doc["experiment"][0]["imageset"] = 0
    doc["imageset"] = [
        {"__id__": "ImageSequence", "template": str(tmp_path / "series.nxs")}
    ]
    json.dump(doc, open(tmp_path / "l.expt", "w"))
    run = subprocess.run(
        [MXI, "find", str(tmp_path / "l.expt"), "-o", str(tmp_path / "strong.rflx")],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    with h5py.File(tmp_path / "strong.rflx") as f:
        assert (
            "32008" in f["summary/steps/0/spots_per_image"]._filters
        ), "compressed, as it must be"
    # In a fresh process, as mxi ui is: hdf5plugin registers its filters for
    # the whole process, and this one has them already (from the fixtures),
    # so read here a missing import would pass unseen -- it had.
    src = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"
    )
    probe = (
        "from mxi_ui import results\n"
        f"path, rec = results.record({str(tmp_path)!r}, 'find')\n"
        "print(rec['program'], len(rec['spots_per_image'][0]))\n"
    )
    fresh = subprocess.run(
        [sys.executable, "-c", probe],
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": src},
    )
    assert fresh.returncode == 0, fresh.stderr[-500:]
    assert fresh.stdout.split() == ["mxi_find", "300"]


def test_the_correlation_tool_reads_the_scaled_file_and_writes_its_json():
    from mxi_ui.steps import by_key

    step = by_key("correlation")
    assert step.inputs == ("scaled.rflx",) and step.output == "mxi_correlation.json"
    from mxi_ui.model import build

    # --json writes nothing unless given: always given, as -o need not be.
    assert build(step, step.inputs, {}, output=step.output) == [
        "mxi",
        "correlation",
        "scaled.rflx",
        "--json",
        "mxi_correlation.json",
    ]
    find = by_key("find")
    assert "-o" not in build(
        find, find.inputs, {}, output=find.output
    ), "-o's default unsaid"
