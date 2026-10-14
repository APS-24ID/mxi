# Regression runs

STATUS: built, 13 October 2026: run, list, compare, tag; bisecting and the
history pages next. Graeme's decisions below.

Two of October's regressions -- the profile fit's variance floor, which wrecked
strong reflections on a near-empty background, and the cell's e.s.d.s, which
made refinement very slow on four sweeps -- were found only on Graeme's data.
A run of the whole chain over real data sets, every figure kept and each run
compared with the last, finds such things the day they happen.

## Graeme's decisions

* About ten data sets to start, of every kind, more later; two to five
  minutes each, so every run takes them all (tiers are there for later).
* Running a commit, HEAD by default: each machine says how it builds.
* Everything in `~/.mxi_regression` (or `$MXI_REGRESSION`).
* Tolerances conservative to start, to be set from the history.

## The directory

    datasets.toml    the data sets: images, a step's options, a tier
    build.toml       this machine's build: the source, CMake's options, jobs,
                     environment
    builds/<commit>  each commit built once -- a worktree, its submodules
                     linked from the source's, configured and built -- kept
    runs/<run>/<ds>  each step's log, the .rflx files, mxi_report.html
    history.sqlite   runs (commit, uncommitted changes, host, platform, build,
                     tag), steps (status, wall and CPU seconds, peak memory),
                     metrics (every figure)

`mxeq regress init` writes commented examples of the two files, and never
overwrites one that exists.

`clean` removes what is large and old (Graeme): the files of every run but the
most recent five (`--keep-runs`) and the tagged ones -- every run's figures
stay in the history, so `list` and `compare` go on working -- and every build
but those the kept runs used and the newest two (`--keep-builds`), each
worktree removed through git and the repository's records pruned. `--intermediates`
takes kept runs' files too, all but scaled.rflx, the report and the logs --
every step's .rflx before scaling, most of a run's size, nothing a comparison
needs. Every run and build is listed with its size, kept or not -- Graeme's
first clean said only "nothing to remove" over 17 GB. Without `--yes` nothing
is removed.

## The commands

    mxeq regress run [--commit REV | --bin DIR] [--datasets A,B]
                     [--quick | --slow | --tier T]
    mxeq regress list
    mxeq regress compare [RUN] [--against previous|reference|RUN]
    mxeq regress tag RUN reference
    mxeq regress clean [--keep-runs N] [--keep-builds N] [--intermediates] [--yes]

Tiers are in order, quick < full < slow, a data set with none full: `run`
takes quick and full, `--quick` the quick ones, `--slow` all three (the
36000-image set, ten minutes, only when asked), `--tier T` one tier alone, and
`--datasets` the named ones whatever their tier. A tier not one of the three is
an error.

`run` first runs each program once with --help, untimed: macOS checks a
binary the first time it runs, some 2 s over the chain that would otherwise
count against the first data set (Graeme). Then it takes the chain, mxi_import
to mxi_scale and mxi_report, a step at a time, stopping at a failure, which is recorded as such. Each step's peak memory
and CPU time are its own, from wait4 (kB on Linux, bytes on macOS, kept in MB). Each data set's
total, import to scale without the report -- images to scaled data, the
figure quoted for mxi -- is printed and kept as `time.total`, its CPU time as
`cpu.total`.
The figures come from the last .rflx's processing summary: spots; indexed
fraction; refinement's RMSDs; integration's <I/sigma> and profile
correlation; scaling's completeness, multiplicity, <I/sigma>, R merge, R pim,
CC half and the anomalous figures overall and in the outer shell, the error
model, the resolution, observations, unique reflections, outliers; the cell
from two-theta; and every step's time and memory.

`compare` takes a run against the one before it on the same host -- times are
only comparable on one machine -- or the host's reference, or a named run. The
accuracy figures are the same run to run on one build, so any change is shown,
better or worse by the figure's sense or just changed (the cell, the error
model), and one past 1 per cent marked; times and memory only past 10 per cent
and half a second or 10 MB. Its exit status is 1 if anything is worse past
tolerance or a data set no longer runs through, for scripts.

## Tried

Insulin's 300 images, run twice with the same programs: the chain in about a
minute, 52 figures; compared, the accuracy figures identical, only the times
differing -- this machine's noise, the spot finder 32.5 s then 19.4.

## Next

* `mxeq regress bisect --dataset D --figure F GOOD BAD`: each commit built and
  that data set run until the step where the figure changed is found.
* A history page: each figure's history for each data set, plotly, as
  mxi_report draws.
* DIALS's figures for each data set, kept once, a column beside mxi's.
