"""What a step of the workflow is, and the command its form builds: no
tkinter here, so that it is tested anywhere."""

import shlex
from dataclasses import dataclass


@dataclass(frozen=True)
class Option:
    """One of a command's options, as a field of its form.

    flag    -- as the command takes it: "--d-min"
    label   -- as the form shows it
    kind    -- "flag" (a tick box: the flag alone), "number", "text", or
               "choice" (one of `choices`)
    about   -- a line on what it does
    """

    flag: str
    label: str
    kind: str = "text"
    about: str = ""
    choices: tuple = ()


@dataclass(frozen=True)
class Step:
    """One mxi command as the interface offers it.

    key       -- the mxi command: "find" runs mxi find
    title     -- as the sidebar shows it
    about     -- one line on what it does, for the form and for teaching
    inputs    -- the files it reads by default, the previous step's outputs
    output    -- the file it writes by default; done when it exists
    options   -- its main options, as the form's fields (the rest by hand)
    output_flag -- the option naming its output, or "" where it has none
    output_asked -- the output written only when named (--json): the flag
                 always given, where -o's default needs none
    tool      -- not a step of the main sequence but a tool beside it
    """

    key: str
    title: str
    about: str
    inputs: tuple = ()
    output: str = ""
    options: tuple = ()
    output_flag: str = "-o"
    tool: bool = False
    output_asked: bool = False

    def command(self, inputs=None, extra=()):
        """The command line: mxi KEY INPUTS... EXTRA..., as a list."""
        return ["mxi", self.key, *(self.inputs if inputs is None else inputs), *extra]


def _number(text):
    try:
        float(text)
        return True
    except ValueError:
        return False


def problems(step, values, threads=""):
    """What is wrong with a form's values: [] where nothing is."""
    found = []
    for option in step.options:
        value = str(values.get(option.flag, "")).strip()
        if not value or option.kind == "flag":
            continue
        if option.kind == "number" and not _number(value):
            found.append(f"{option.label}: a number, not {value!r}")
        if option.kind == "choice" and value not in option.choices:
            found.append(f"{option.label}: one of {', '.join(option.choices)}")
    if threads.strip() and not threads.strip().isdigit():
        found.append(
            f"Threads: a whole number, 0 or empty for every core, not {threads!r}"
        )
    return found


def build(step, inputs, values, output="", threads="", extra=""):
    """The command a form says: mxi KEY INPUTS, then its options as set --
    a tick box its flag alone, the rest their flag and value, an empty field
    left out -- then the output where it differs from the default, the
    threads, and anything typed as further options."""
    argv = ["mxi", step.key, *inputs]
    for option in step.options:
        value = values.get(option.flag, "")
        if option.kind == "flag":
            if value:
                argv.append(option.flag)
        elif str(value).strip():
            argv += [option.flag, str(value).strip()]
    if (
        step.output_flag
        and output.strip()
        and (output.strip() != step.output or step.output_asked)
    ):
        argv += [step.output_flag, output.strip()]
    if threads.strip() and threads.strip() != "0":
        argv += ["-j", threads.strip()]
    if extra.strip():
        argv += shlex.split(extra)
    return argv
