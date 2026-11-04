"""What a step of the workflow is: no tkinter here."""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Step:
    """One mxi command as the interface offers it.

    key       -- the mxi command: "find" runs mxi find
    title     -- as the sidebar shows it
    about     -- one line on what it does, for the form and for teaching
    inputs    -- the files it reads by default, the previous step's outputs
    output    -- the file it writes by default; done when it exists
    tool      -- not a step of the main sequence but a tool beside it
    """

    key: str
    title: str
    about: str
    inputs: tuple = ()
    output: str = ""
    tool: bool = False
    options: tuple = field(default=())

    def command(self, inputs=None, extra=()):
        """The command line: mxi KEY INPUTS... EXTRA..., as a list."""
        return ["mxi", self.key, *(self.inputs if inputs is None else inputs), *extra]
