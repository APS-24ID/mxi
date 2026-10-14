"""Regression runs over real data sets (docs/regression.md).

Everything lives in one directory, ~/.mxi_regression unless MXI_REGRESSION
says otherwise:

    datasets.toml    the data sets: their images, any options a step needs,
                     a tier
    build.toml       how this machine builds mxi: the source, CMake's options
    builds/<commit>  each commit built once, in a worktree, and kept
    runs/<run>/<ds>  a run's logs, files and report, a directory a data set
    history.sqlite   every run: its identity, each step's result, time and
                     memory, and the figures from its processing summary

    mxeq regress init                 write example datasets.toml, build.toml
    mxeq regress run [--commit REV | --bin DIR] [--datasets A,B] [--tier T]
    mxeq regress list
    mxeq regress compare [RUN] [--against RUN|reference]
    mxeq regress tag RUN NAME

Python only, the standard library and mxeq, so the same on macOS and Linux.
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import platform
import shutil
import socket
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

try:
    import tomllib  # Python 3.11
except ModuleNotFoundError:  # pragma: no cover - 3.10
    try:
        import tomli as tomllib
    except ModuleNotFoundError:
        tomllib = None

# The chain, each step's program, input and output; a data set's options for a
# step are added after these.
STEPS = [
    ("import", "mxi_import", None, "imported.rflx"),
    ("find", "mxi_find", "imported.rflx", "strong.rflx"),
    ("index", "mxi_index", "strong.rflx", "indexed.rflx"),
    ("refine", "mxi_refine", "indexed.rflx", "refined.rflx"),
    ("integrate", "mxi_integrate", "refined.rflx", "integrated.rflx"),
    ("symmetry", "mxi_symmetry", "integrated.rflx", "symmetrized.rflx"),
    ("scale", "mxi_scale", "symmetrized.rflx", "scaled.rflx"),
    ("report", "mxi_report", "scaled.rflx", "mxi_report.html"),
]

EXAMPLE_DATASETS = """\
# The data sets a regression run takes (docs/regression.md). Each a table:
#   images   one master file or several (several sweeps), absolute paths
#   tier     "quick", "full" (the default) or "slow": `mxeq regress run`
#            takes quick and full, --quick the quick ones, --slow all three
#   options  a step's extra arguments, by step: import find index refine
#            integrate symmetry scale
#
# [insulin]
# images = ["/data/insulin/ins10_1.nxs"]
# tier = "quick"
#
# [two-crystals]
# images = ["/data/ins2_27/ins2_27.nxs"]
# tier = "full"
# options.index = ["--max-lattices", "2"]
"""

EXAMPLE_BUILD = """\
# How this machine builds mxi, for `mxeq regress run --commit REV`: the commit
# is checked out in a worktree under builds/ and built with these.
source = "/path/to/mxi"          # the git repository
cmake = ["-DCMAKE_BUILD_TYPE=Release"]
# cmake = ["-DCMAKE_BUILD_TYPE=Release", "-DSPOTFINDER_METAL=ON", "-DMXI_FFTW=ON"]
jobs = 8
# [environment]                  # set for cmake and the build
# CC = "clang"
# CXX = "clang++"
"""


# --- Where things are --------------------------------------------------------


def home() -> Path:
    return Path(os.environ.get("MXI_REGRESSION", Path.home() / ".mxi_regression"))


def read_toml(path: Path) -> dict:
    if tomllib is None:
        raise SystemExit("mxeq regress needs Python 3.11, or the tomli package")
    if not path.exists():
        raise SystemExit(f"{path} is missing: `mxeq regress init` writes an example")
    with open(path, "rb") as f:
        return tomllib.load(f)


def database() -> sqlite3.Connection:
    home().mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(home() / "history.sqlite")
    db.executescript("""
        CREATE TABLE IF NOT EXISTS runs (
            id TEXT PRIMARY KEY, started TEXT, finished TEXT, commit_id TEXT,
            dirty INTEGER, host TEXT, platform TEXT, build TEXT, tag TEXT);
        CREATE TABLE IF NOT EXISTS steps (
            run TEXT, dataset TEXT, step TEXT, status INTEGER,
            wall REAL, cpu REAL, memory_mb REAL);
        CREATE TABLE IF NOT EXISTS metrics (
            run TEXT, dataset TEXT, name TEXT, value REAL);
        """)
    return db


# --- Building a commit -------------------------------------------------------


def git(source: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(source), *args], check=True, capture_output=True, text=True
    ).stdout.strip()


def build_commit(rev: str) -> tuple[Path, str]:
    """The directory with the programs built from `rev`, building it first if
    this machine has not: a worktree of the source at the commit, its
    submodules linked from the source's own, configured and built as
    build.toml says. Each commit is built once."""
    config = read_toml(home() / "build.toml")
    source = Path(config["source"]).expanduser()
    commit = git(source, "rev-parse", rev)
    where = home() / "builds" / commit[:12]
    built = where / "build"
    if (where / "ok").exists():
        return built, commit
    tree = where / "src"
    if not tree.exists():
        where.mkdir(parents=True, exist_ok=True)
        git(source, "worktree", "add", "--detach", str(tree), commit)
        # The submodules from the source's checkout, not the network.
        for sub in (source / "third_party").iterdir():
            mine = tree / "third_party" / sub.name
            if sub.is_dir() and any(sub.iterdir()):
                if mine.exists() and not any(mine.iterdir()):
                    mine.rmdir()
                if not mine.exists():
                    mine.symlink_to(sub.resolve())
    built.mkdir(parents=True, exist_ok=True)
    env = dict(
        os.environ, **{k: str(v) for k, v in config.get("environment", {}).items()}
    )
    log = open(where / "build.log", "w")
    print(f"building {commit[:12]} in {built}", flush=True)
    for words in (
        ["cmake", str(tree), *config.get("cmake", [])],
        ["cmake", "--build", ".", "-j", str(config.get("jobs", 4))],
    ):
        r = subprocess.run(
            words, cwd=built, env=env, stdout=log, stderr=subprocess.STDOUT
        )
        if r.returncode != 0:
            raise SystemExit(
                f"the build of {commit[:12]} failed: {where / 'build.log'}"
            )
    (where / "ok").write_text(commit + "\n")
    return built, commit


# --- Running a data set ------------------------------------------------------


def run_step(words: list[str], cwd: Path, log: Path) -> tuple[int, float, float, float]:
    """One program: its exit status, wall and CPU seconds, and peak memory in
    MB -- its own, from wait4, on Linux (kB) and macOS (bytes) alike."""
    start = time.perf_counter()
    with open(log, "w") as out:
        p = subprocess.Popen(words, cwd=cwd, stdout=out, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(p.pid, 0)
    wall = time.perf_counter() - start
    code = os.waitstatus_to_exitcode(status)
    p.returncode = code
    rss = usage.ru_maxrss / (1024 * 1024 if sys.platform == "darwin" else 1024)
    return code, wall, usage.ru_utime + usage.ru_stime, rss


def warm_up(bin_dir: Path) -> None:
    """Each program run once, asked only for its help, before anything is
    timed: macOS checks a binary the first time it runs, a fraction of a second
    each, some 2 s over the chain that would otherwise be counted against the
    first data set (Graeme). --help, so no program writes its log; in a scratch
    directory besides. Harmless elsewhere."""
    import tempfile

    with tempfile.TemporaryDirectory() as scratch:
        for _, program, _, _ in STEPS:
            exe = bin_dir / program
            if exe.exists():
                try:
                    subprocess.run(
                        [str(exe), "--help"],
                        cwd=scratch,
                        stdin=subprocess.DEVNULL,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                        timeout=60,
                    )
                except (OSError, subprocess.TimeoutExpired):
                    pass


def run_dataset(name: str, spec: dict, bin_dir: Path, where: Path) -> list[tuple]:
    where.mkdir(parents=True, exist_ok=True)
    options = spec.get("options", {})
    rows = []
    for step, program, given, written in STEPS:
        exe = bin_dir / program
        words = [str(exe)]
        words += (
            [str(Path(i).expanduser()) for i in spec["images"]]
            if given is None
            else [given]
        )
        words += ["-o", written, *[str(w) for w in options.get(step, [])]]
        code, wall, cpu, rss = run_step(words, where, where / f"{step}.log")
        rows.append((step, code, wall, cpu, rss))
        print(
            f"  {name:20s} {step:10s} {'ok' if code == 0 else 'FAILED':6s} {wall:7.1f} s",
            flush=True,
        )
        if code != 0:
            break
    return rows


# --- The figures -------------------------------------------------------------


def _row(table: dict, label: str, column: str):
    try:
        i = list(table["shell"]).index(label)
        return float(table[column][i])
    except (KeyError, ValueError, IndexError, TypeError):
        return None


def figures(summary: dict) -> dict[str, float]:
    """What a run is judged by, from a file's processing summary: each
    program's last step's figures."""
    last = {}
    for s in summary.get("steps", []):
        last[s.get("program")] = s
    out: dict[str, float] = {}

    def put(name, value):
        try:
            if value is not None:
                out[name] = float(value)
        except (TypeError, ValueError):
            pass

    f = last.get("mxi_find", {})
    put("find.spots", f.get("spots"))
    x = last.get("mxi_index", {})
    put("index.indexed_fraction", x.get("indexed_fraction"))
    if x.get("indexed") is not None and x.get("reflections"):
        put("index.indexed_fraction", float(x["indexed"]) / float(x["reflections"]))
    r = last.get("mxi_refine", {})
    for k, v in zip(("x_px", "y_px", "z_images"), r.get("rmsd_px_px_images", []) or []):
        put(f"refine.rmsd_{k}", v)
    i = last.get("mxi_integrate", {})
    overall = i.get("overall", {})
    for k in ("isigma_prf", "isigma_sum", "cc_prf"):
        try:
            put(f"integrate.{k}", overall[k][0])
        except (KeyError, IndexError, TypeError):
            pass
    s = last.get("mxi_scale", {})
    m = s.get("merging_summary", {})
    for shell in ("overall", "high"):
        for k in (
            "completeness",
            "multiplicity",
            "i_over_sigma",
            "rmerge",
            "rpim",
            "cc_half",
            "cc_anom",
            "anom_slope",
            "anom_slope_random",
        ):
            put(f"scale.{shell}.{k}", _row(m, shell, k))
    for k in (
        "error_model_a",
        "error_model_a_background",
        "error_model_b",
        "resolution_cc_half",
        "observations",
        "unique",
        "outliers",
    ):
        put(f"scale.{k}", s.get(k))
    for k, v in zip(
        ("a", "b", "c", "alpha", "beta", "gamma"), s.get("two_theta_cell", []) or []
    ):
        put(f"cell.{k}", v)
    return out


# Each figure's sense: +1 higher is better, -1 lower, 0 neither -- a change is
# flagged, not judged. Matched by the name's start, the longest first.
SENSE = {
    "index.indexed_fraction": +1,
    "refine.rmsd": -1,
    "integrate.isigma": +1,
    "integrate.cc_prf": +1,
    "scale.overall.completeness": +1,
    "scale.overall.i_over_sigma": +1,
    "scale.overall.cc_half": +1,
    "scale.overall.cc_anom": +1,
    "scale.overall.rmerge": -1,
    "scale.overall.rpim": -1,
    "scale.high.completeness": +1,
    "scale.high.i_over_sigma": +1,
    "scale.high.cc_half": +1,
    "scale.high.rmerge": -1,
    "scale.high.rpim": -1,
    "time.": -1,
    "cpu.": -1,
    "memory.": -1,
}


def sense_of(name: str) -> int:
    best = ""
    for prefix in SENSE:
        if name.startswith(prefix) and len(prefix) > len(best):
            best = prefix
    return SENSE.get(best, 0)


def tolerance_of(name: str) -> tuple[float, float]:
    """Conservative to start (Graeme): the accuracy figures are the same run to
    run on one build, so any change at all is shown, and one past 1 per cent
    marked; times and memory past 10 per cent and half a second or 10 MB."""
    if name.startswith(("time.", "cpu.")):
        return 0.10, 0.5
    if name.startswith("memory."):
        return 0.10, 10.0
    return 0.01, 0.0


# --- The commands ------------------------------------------------------------


TIERS = ("quick", "full", "slow")


def choose(datasets: dict, names=None, tier=None, upto="full") -> dict:
    """The data sets a run takes. Named ones whatever their tier; else one
    tier exactly (--tier), or every tier up to `upto` -- quick < full < slow:
    --quick the quick ones, the default quick and full, --slow all of them
    (Graeme: the 36000-image set, ten minutes, only when asked). No tier is
    full; a tier not one of the three is an error, not a data set never run."""
    out = {}
    for k, v in datasets.items():
        if not isinstance(v, dict):
            continue
        own = v.get("tier", "full")
        if own not in TIERS:
            raise SystemExit(f"{k}: tier {own!r} is not one of {', '.join(TIERS)}")
        if names is not None:
            if k in names:
                out[k] = v
        elif tier is not None:
            if own == tier:
                out[k] = v
        elif TIERS.index(own) <= TIERS.index(upto):
            out[k] = v
    if names is not None:
        missing = [n for n in names if n not in out]
        if missing:
            raise SystemExit(f"no such data set: {', '.join(missing)}")
    return out


def cmd_init(args) -> int:
    home().mkdir(parents=True, exist_ok=True)
    for name, text in (
        ("datasets.toml", EXAMPLE_DATASETS),
        ("build.toml", EXAMPLE_BUILD),
    ):
        path = home() / name
        if path.exists():  # never overwritten: it may hold the user's edits
            print(f"{path} exists, so left as it is (not overwritten)")
        else:
            path.write_text(text)
            print(f"wrote {path}")
    return 0


def cmd_run(args) -> int:
    datasets = read_toml(home() / "datasets.toml")
    chosen = choose(
        datasets,
        names=args.datasets.split(",") if args.datasets else None,
        tier=args.tier,
        upto="quick" if args.quick else "slow" if args.slow else "full",
    )
    if not chosen:
        raise SystemExit("no data sets chosen")
    if args.bin:
        bin_dir = Path(args.bin).expanduser().resolve()
        commit, dirty, build = "", 0, f"bin {bin_dir}"
        try:  # the binaries' own source, if they sit in its build directory
            source = Path(read_toml(home() / "build.toml")["source"]).expanduser()
            commit = git(source, "rev-parse", "HEAD")
            dirty = (
                1 if git(source, "status", "--porcelain", "--untracked-files=no") else 0
            )
        except (SystemExit, Exception):
            pass
    else:
        bin_dir, commit = build_commit(args.commit)
        dirty, build = 0, "built " + " ".join(
            read_toml(home() / "build.toml").get("cmake", [])
        )
    started = datetime.datetime.now()
    run_id = started.strftime("%Y%m%d-%H%M%S") + (f"-{commit[:8]}" if commit else "")
    db = database()
    db.execute(
        "INSERT INTO runs VALUES (?, ?, NULL, ?, ?, ?, ?, ?, NULL)",
        (
            run_id,
            started.isoformat(timespec="seconds"),
            commit,
            dirty,
            socket.gethostname(),
            f"{platform.system()} {platform.machine()}",
            build,
        ),
    )
    db.commit()
    print(f"run {run_id}: {len(chosen)} data sets, programs from {bin_dir}")
    warm_up(bin_dir)
    from . import rflx

    failed = 0
    for name, spec in chosen.items():
        where = home() / "runs" / run_id / name
        rows = run_dataset(name, spec, bin_dir, where)
        for step, code, wall, cpu, rss in rows:
            db.execute(
                "INSERT INTO steps VALUES (?, ?, ?, ?, ?, ?, ?)",
                (run_id, name, step, code, wall, cpu, rss),
            )
            db.execute(
                "INSERT INTO metrics VALUES (?, ?, ?, ?)",
                (run_id, name, f"time.{step}", wall),
            )
            db.execute(
                "INSERT INTO metrics VALUES (?, ?, ?, ?)",
                (run_id, name, f"memory.{step}", rss),
            )
        # The whole: images to scaled data, the figure quoted for mxi, without
        # the report; and the CPU time, for how well the threads are used.
        chain = [r for r in rows if r[0] != "report"]
        if chain:
            total = sum(r[2] for r in chain)
            cpu = sum(r[3] for r in chain)
            print(
                f"  {name:20s} {'total':10s} {'':6s} {total:7.1f} s"
                f"  (CPU {cpu:.1f} s)",
                flush=True,
            )
            db.execute(
                "INSERT INTO metrics VALUES (?, ?, ?, ?)",
                (run_id, name, "time.total", total),
            )
            db.execute(
                "INSERT INTO metrics VALUES (?, ?, ?, ?)",
                (run_id, name, "cpu.total", cpu),
            )
        if rows and rows[-1][1] != 0:
            failed += 1
        last = None
        for _, _, _, written in reversed(STEPS):
            if written.endswith(".rflx") and (where / written).exists():
                last = where / written
                break
        if last is not None:
            # The figures, or the data set failed at them -- a file not to be
            # read is one data set's failure, not the run's.
            try:
                summary = rflx.summary_dict(str(last)) or {}
            except Exception as e:  # noqa: BLE001
                print(f"  {name:20s} figures    FAILED   {e}", flush=True)
                db.execute(
                    "INSERT INTO steps VALUES (?, ?, ?, ?, ?, ?, ?)",
                    (run_id, name, "figures", 1, 0.0, 0.0, 0.0),
                )
                failed += 1
                summary = {}
            for k, v in figures(summary).items():
                db.execute(
                    "INSERT INTO metrics VALUES (?, ?, ?, ?)", (run_id, name, k, v)
                )
        db.commit()
    db.execute(
        "UPDATE runs SET finished = ? WHERE id = ?",
        (datetime.datetime.now().isoformat(timespec="seconds"), run_id),
    )
    db.commit()
    print(
        f"run {run_id} done: {len(chosen) - failed} of {len(chosen)} data sets ran through"
    )
    return 1 if failed else 0


def _size(path: Path) -> int:
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            try:
                total += os.lstat(os.path.join(root, f)).st_size
            except OSError:
                pass
    return total


def _human(n: float) -> str:
    for unit in ("B", "kB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


# What a kept run can do without: every step's file but the last, scaled.rflx,
# which the report and the figures came from.
INTERMEDIATES = [
    w for _, _, _, w in STEPS if w.endswith(".rflx") and w != "scaled.rflx"
]


def cmd_clean(args) -> int:
    """Remove the files of old runs and old builds, which are large -- and,
    with --intermediates, kept runs' files but scaled.rflx, the report and the
    logs. The figures stay in the history, so list and compare go on working.
    Everything is shown with its size, kept or not (Graeme: 'nothing to
    remove' over 17 GB). Without --yes, nothing is removed."""
    db = database()
    runs = [
        r for r in db.execute("SELECT id, commit_id, tag FROM runs ORDER BY started")
    ]
    kept = {r[0] for r in runs[-args.keep_runs :]} if args.keep_runs > 0 else set()
    kept |= {r[0] for r in runs if r[2]}  # tagged runs keep their files
    tags = {r[0]: r[2] for r in runs}
    plan = []  # (kind, path, remove?, why)
    run_dir = home() / "runs"
    if run_dir.exists():
        for d in sorted(run_dir.iterdir()):
            if not d.is_dir():
                continue
            if d.name not in kept:
                plan.append(("run", d, True, "old"))
                continue
            why = "tagged " + tags[d.name] if tags.get(d.name) else "recent"
            plan.append(("run", d, False, why))
            if args.intermediates:
                for ds in sorted(x for x in d.iterdir() if x.is_dir()):
                    for name in INTERMEDIATES:
                        f = ds / name
                        if f.exists():
                            plan.append(("file", f, True, "intermediate"))
    wanted = {(r[1] or "")[:12] for r in runs if r[0] in kept}
    build_dir = home() / "builds"
    builds = []
    if build_dir.exists():
        builds = sorted(
            (d for d in build_dir.iterdir() if d.is_dir()),
            key=lambda d: d.stat().st_mtime,
        )
    newest = (
        {d.name for d in builds[-args.keep_builds :]} if args.keep_builds > 0 else set()
    )
    for d in builds:
        if d.name in wanted:
            plan.append(("build", d, False, "used by a kept run"))
        elif d.name in newest:
            plan.append(("build", d, False, "newest"))
        else:
            plan.append(("build", d, True, "old"))
    removed_total = kept_total = 0
    for kind, d, remove, why in plan:
        n = _size(d) if d.is_dir() else d.stat().st_size
        if kind == "file":
            label = f"{d.parent.parent.name}/{d.parent.name}/{d.name}"
        else:
            label = d.name
        print(
            f"  {'remove' if remove else 'keep':6s} {kind:5s} {label:48s} "
            f"{_human(n):>10s}  {why}"
        )
        if remove:
            removed_total += n
        elif kind != "file":
            kept_total += n
    # A kept run's size above includes its intermediates; count them once.
    kept_total -= sum(d.stat().st_size for k, d, r, _ in plan if k == "file")
    removal = [(k, d) for k, d, r, _ in plan if r]
    print(
        f"{len(removal)} to remove, {_human(removed_total)}; "
        f"{_human(max(kept_total, 0))} kept"
        + ("" if args.yes or not removal else "; --yes to remove them")
    )
    if not removal:
        if not args.intermediates:
            print(
                "  (--intermediates would also take kept runs' files but "
                "scaled.rflx, the report and the logs; --keep-runs N fewer runs)"
            )
        return 0
    if not args.yes:
        return 0
    source = None
    try:
        source = Path(read_toml(home() / "build.toml")["source"]).expanduser()
    except (SystemExit, Exception):
        pass
    for kind, d in removal:
        if kind == "file":
            d.unlink(missing_ok=True)
            continue
        tree = d / "src"
        if kind == "build" and source is not None and tree.exists():
            # The worktree, through git: the repository's own record of it too.
            subprocess.run(
                ["git", "-C", str(source), "worktree", "remove", "--force", str(tree)],
                capture_output=True,
            )
        shutil.rmtree(d, ignore_errors=True)
    if source is not None:
        subprocess.run(
            ["git", "-C", str(source), "worktree", "prune"], capture_output=True
        )
    print(f"removed {len(removal)}, {_human(removed_total)}")
    return 0


def cmd_list(args) -> int:
    db = database()
    for row in db.execute(
        "SELECT id, commit_id, dirty, host, tag, "
        "(SELECT COUNT(DISTINCT dataset) FROM steps WHERE run = runs.id) "
        "FROM runs ORDER BY started"
    ):
        run_id, commit, dirty, host, tag, n = row
        print(
            f"{run_id:28s} {(commit or '')[:10]:10s}{'+' if dirty else ' '} {host:20s} "
            f"{n:3d} data sets  {tag or ''}"
        )
    return 0


def cmd_tag(args) -> int:
    db = database()
    if args.name == "reference":  # one reference a host
        host = db.execute("SELECT host FROM runs WHERE id = ?", (args.run,)).fetchone()
        if host:
            db.execute(
                "UPDATE runs SET tag = NULL WHERE tag = 'reference' AND host = ?", host
            )
    db.execute("UPDATE runs SET tag = ? WHERE id = ?", (args.name, args.run))
    db.commit()
    return 0


def _metrics(db, run_id) -> dict[tuple[str, str], float]:
    return {
        (d, n): v
        for d, n, v in db.execute(
            "SELECT dataset, name, value FROM metrics WHERE run = ?", (run_id,)
        )
    }


def _status(db, run_id) -> dict[str, str]:
    out = {}
    for d, step, code in db.execute(
        "SELECT dataset, step, status FROM steps WHERE run = ? ORDER BY rowid",
        (run_id,),
    ):
        out[d] = "ok" if code == 0 else f"failed in {step}"
    return out


def cmd_compare(args) -> int:
    db = database()
    runs = [r[0] for r in db.execute("SELECT id FROM runs ORDER BY started")]
    if not runs:
        raise SystemExit("no runs yet")
    current = args.run or runs[-1]
    host = db.execute("SELECT host FROM runs WHERE id = ?", (current,)).fetchone()[0]
    if args.against in (None, "previous"):
        earlier = [
            r
            for r in db.execute(
                "SELECT id FROM runs WHERE host = ? AND started < "
                "(SELECT started FROM runs WHERE id = ?) ORDER BY started",
                (host, current),
            )
        ]
        if not earlier:
            raise SystemExit(f"no earlier run on {host} to compare {current} with")
        baseline = earlier[-1][0]
    elif args.against == "reference":
        row = db.execute(
            "SELECT id FROM runs WHERE tag = 'reference' AND host = ?", (host,)
        ).fetchone()
        if not row:
            raise SystemExit(
                f"no reference on {host}: `mxeq regress tag RUN reference`"
            )
        baseline = row[0]
    else:
        baseline = args.against
    a, b = _metrics(db, baseline), _metrics(db, current)
    sa, sb = _status(db, baseline), _status(db, current)
    print(f"{current} against {baseline}")
    worse = 0
    for d in sorted(set(sa) | set(sb)):
        lines = []
        if sa.get(d) != sb.get(d):
            lines.append(
                f"    status: {sa.get(d, 'not run')} -> {sb.get(d, 'not run')}"
            )
            if sb.get(d, "").startswith("failed"):
                worse += 1
        for name in sorted({n for dd, n in set(a) | set(b) if dd == d}):
            x, y = a.get((d, name)), b.get((d, name))
            if x is None or y is None:
                lines.append(f"    {name}: {x} -> {y}")
                continue
            if x == y:
                continue
            rel, floor = tolerance_of(name)
            change = y - x
            relative = abs(change) / max(abs(x), 1e-12)
            timing = name.startswith(("time.", "cpu.", "memory."))
            past = relative > rel and abs(change) > floor
            if timing and not past:
                continue
            sense = sense_of(name)
            verdict = (
                ("better" if change * sense > 0 else "worse") if sense else "changed"
            )
            mark = (
                "!!" if past and verdict != "better" else ("  " if not past else "+ ")
            )
            if past and verdict == "worse":
                worse += 1
            lines.append(
                f"  {mark}{name}: {x:.6g} -> {y:.6g} ({change / x * 100 if x else 0:+.2f}%) {verdict}"
            )
        print(f"{d}: {'the same' if not lines else ''}")
        for line in lines:
            print(line)
    print(f"{worse} worse past tolerance" if worse else "nothing worse past tolerance")
    return 1 if worse else 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="mxeq regress",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("init", help="write example datasets.toml and build.toml")
    r = sub.add_parser("run", help="run the chain on the data sets, keep the results")
    g = r.add_mutually_exclusive_group()
    g.add_argument("--commit", default="HEAD", help="build and run this commit (HEAD)")
    g.add_argument("--bin", help="run the programs already built in this directory")
    r.add_argument("--datasets", help="these, comma-separated (all)")
    speed = r.add_mutually_exclusive_group()
    speed.add_argument("--quick", action="store_true", help="only the quick data sets")
    speed.add_argument(
        "--slow",
        action="store_true",
        help="every data set, the slow ones too (the default: " "quick and full)",
    )
    speed.add_argument("--tier", choices=TIERS, help="only the data sets of this tier")
    sub.add_parser("list", help="the runs so far")
    cl = sub.add_parser(
        "clean",
        help="remove old runs' files and old builds " "(says what, unless --yes)",
    )
    cl.add_argument(
        "--keep-runs",
        type=int,
        default=5,
        help="the most recent runs whose files stay (5); tagged "
        "runs' always stay; every run's figures stay",
    )
    cl.add_argument(
        "--keep-builds",
        type=int,
        default=2,
        help="the newest builds that stay (2), besides those the " "kept runs used",
    )
    cl.add_argument(
        "--intermediates",
        action="store_true",
        help="also kept runs' files but scaled.rflx, the report "
        "and the logs: every step's .rflx before scaling",
    )
    cl.add_argument("--yes", action="store_true", help="remove them")
    c = sub.add_parser("compare", help="a run against an earlier one")
    c.add_argument("run", nargs="?", help="the run (the latest)")
    c.add_argument("--against", help="a run, 'previous' (the default) or 'reference'")
    t = sub.add_parser(
        "tag", help="name a run, 'reference' the baseline a host compares with"
    )
    t.add_argument("run")
    t.add_argument("name")
    args = p.parse_args(argv)
    return {
        "init": cmd_init,
        "run": cmd_run,
        "list": cmd_list,
        "compare": cmd_compare,
        "clean": cmd_clean,
        "tag": cmd_tag,
    }[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
