"""Which steps a directory has done already: those whose output is there --
so that a directory processed before, here or by hand, opens as it stands."""

import os

from mxi_ui.steps import STEPS


def scan(directory):
    """{step key: the output's path, or None} for every step of the sequence."""
    done = {}
    for step in STEPS:
        path = os.path.join(directory, step.output) if step.output else None
        done[step.key] = path if path and os.path.isfile(path) else None
    return done


def next_step(directory):
    """The first step not yet done, or None when all are."""
    done = scan(directory)
    for step in STEPS:
        if not done[step.key]:
            return step.key
    return None
