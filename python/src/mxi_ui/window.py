"""The main window: the steps down the side; the chosen step's form -- its
inputs, main options, output, threads and any further options, the command
they make shown as it is built -- above tabs for its images, log, summary and
plots; a working directory at the top, a status line at the foot. The viewer
and the plots fill their tabs in the phases after (docs/ui.md).
"""

import os
import shlex
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter import font as tkfont

from mxi_ui import state
from mxi_ui.locate import find_mxi
from mxi_ui.model import build, problems
from mxi_ui.runner import Runner, as_text
from mxi_ui.steps import STEPS, TOOLS, by_key
from mxi_ui.tooltip import Tooltip
from mxi_ui.version import VERSION

POLL_MS = 100
# Run's shortcut: Command-Return on a Mac, Control-Return elsewhere.
RUN_KEY = "Command" if sys.platform == "darwin" else "Control"
RUN_KEY_SHOWN = "Cmd" if sys.platform == "darwin" else "Ctrl"
DONE, PENDING = "done", ""


class MainWindow:
    def __init__(self, root, directory=None):
        self.root = root
        self.runner = Runner()
        self.mxi = find_mxi()
        self.directory = tk.StringVar(value=os.path.abspath(directory or os.getcwd()))
        self.status = tk.StringVar(value="Ready")
        self.selected = None
        self.forms = {}  # each step's entries, kept while the window is open
        self.tooltips = []
        self.command = tk.StringVar()
        self.trouble = tk.StringVar()
        root.title(f"mxi ui {VERSION}")
        root.geometry("1100x720")
        root.minsize(800, 500)
        self._menus()
        self._layout()
        self.refresh()
        root.after(POLL_MS, self._poll)
        root.protocol("WM_DELETE_WINDOW", self.quit)
        root.bind_all(f"<{RUN_KEY}-Return>", lambda _e: self.run_selected())

    # --- building ---------------------------------------------------------

    def _menus(self):
        bar = tk.Menu(self.root)
        file = tk.Menu(bar, tearoff=False)
        file.add_command(label="Open directory...", command=self.choose_directory)
        file.add_separator()
        file.add_command(label="Quit", command=self.quit)
        bar.add_cascade(label="File", menu=file)
        run = tk.Menu(bar, tearoff=False)
        run.add_command(
            label="Run the selected step",
            command=self.run_selected,
            accelerator=f"{RUN_KEY_SHOWN}+Return",
        )
        run.add_command(label="Cancel", command=self.cancel)
        bar.add_cascade(label="Run", menu=run)
        help_ = tk.Menu(bar, tearoff=False)
        help_.add_command(label="About mxi ui", command=self.about)
        bar.add_cascade(label="Help", menu=help_)
        self.root.config(menu=bar)

    def _layout(self):
        top = ttk.Frame(self.root, padding=(8, 6))
        top.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(top, text="Directory").pack(side=tk.LEFT)
        ttk.Entry(top, textvariable=self.directory).pack(
            side=tk.LEFT, fill=tk.X, expand=True, padx=6
        )
        ttk.Button(top, text="Browse...", command=self.choose_directory).pack(
            side=tk.LEFT
        )

        foot = ttk.Frame(self.root, padding=(8, 4))
        foot.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Label(foot, textvariable=self.status).pack(side=tk.LEFT)
        self.cancel_button = ttk.Button(
            foot, text="Cancel", command=self.cancel, state=tk.DISABLED
        )
        self.cancel_button.pack(side=tk.RIGHT)

        panes = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        panes.pack(fill=tk.BOTH, expand=True)
        side = ttk.Frame(panes, padding=(6, 0))
        # Rows as tall as the font needs: fixed, they overlapped where Tk draws
        # text larger (a Mac; a higher scaling).
        linespace = tkfont.nametofont("TkDefaultFont").metrics("linespace")
        ttk.Style(self.root).configure("Treeview", rowheight=linespace + 6)
        self.tree = ttk.Treeview(
            side, columns=("state",), show="tree headings", selectmode="browse"
        )
        self.tree.heading("#0", text="Step")
        self.tree.heading("state", text="")
        # As wide as the longest title in the font drawn, its indent beside.
        widest = max(
            tkfont.nametofont("TkDefaultFont").measure(step.title)
            for step in STEPS + TOOLS
        )
        self.tree.column("#0", width=widest + 3 * linespace, minwidth=widest)
        self.tree.column("state", width=60, anchor=tk.CENTER)
        self.tree.insert("", tk.END, iid="steps", text="Processing", open=True)
        for step in STEPS:
            self.tree.insert("steps", tk.END, iid=step.key, text=step.title)
        self.tree.insert("", tk.END, iid="tools", text="Tools", open=True)
        for step in TOOLS:
            self.tree.insert("tools", tk.END, iid=step.key, text=step.title)
        self.tree.pack(fill=tk.BOTH, expand=True)
        self.tree.bind("<<TreeviewSelect>>", self._on_select)
        panes.add(side, weight=0)

        right = ttk.PanedWindow(panes, orient=tk.VERTICAL)
        self.form = ttk.LabelFrame(right, text="", padding=(10, 6))
        right.add(self.form, weight=0)
        self._form_foot()
        self.tabs = ttk.Notebook(right)
        self.images = ttk.Frame(self.tabs, padding=12)
        ttk.Label(
            self.images,
            text="The viewer is the next phase of mxi ui; until then, mxi view.",
        ).pack(anchor=tk.NW)
        self.log = self._text(self.tabs)
        self.summary = self._text(self.tabs)
        self.plots = ttk.Frame(self.tabs, padding=12)
        ttk.Label(self.plots, text="Plots come with the workflow's results.").pack(
            anchor=tk.NW
        )
        self.tabs.add(self.images, text="Images")
        self.tabs.add(self.log.master, text="Log")
        self.tabs.add(self.summary.master, text="Summary")
        self.tabs.add(self.plots, text="Plots")
        right.add(self.tabs, weight=1)
        panes.add(right, weight=1)

    def _form_foot(self):
        """The command the form makes, its problems, Copy and Run: built once,
        below the fields that each step brings."""
        self.fields = ttk.Frame(self.form)
        self.fields.pack(fill=tk.X)
        foot = ttk.Frame(self.form, padding=(0, 8, 0, 0))
        foot.pack(fill=tk.X)
        # Run and Copy first, from the right, so that a narrow window narrows
        # the command, never drops the buttons (as it had the inputs' Browse).
        self.run_button = ttk.Button(foot, text="Run", command=self.run_selected)
        self.run_button.pack(side=tk.RIGHT, padx=(6, 0))
        ttk.Button(foot, text="Copy", command=self.copy_command).pack(
            side=tk.RIGHT, padx=(6, 0)
        )
        shown = ttk.Entry(
            foot, textvariable=self.command, font="TkFixedFont", state="readonly"
        )
        shown.pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(self.form, textvariable=self.trouble, foreground="#b00020").pack(
            anchor=tk.W
        )

    def _form_of(self, step):
        """The step's entries as Tk variables, made the first time it is shown
        and kept: inputs, each option, output, threads, further options."""
        if step.key not in self.forms:
            form = {
                "inputs": tk.StringVar(value=shlex.join(step.inputs)),
                "output": tk.StringVar(value=step.output),
                "threads": tk.StringVar(value=""),
                "extra": tk.StringVar(value=""),
                "options": {
                    o.flag: (
                        tk.BooleanVar(value=False)
                        if o.kind == "flag"
                        else tk.StringVar(value="")
                    )
                    for o in step.options
                },
            }
            for var in [
                form["inputs"],
                form["output"],
                form["threads"],
                form["extra"],
            ] + list(form["options"].values()):
                var.trace_add("write", lambda *_: self._update_command())
            self.forms[step.key] = form
        return self.forms[step.key]

    def _show_form(self, step):
        for tooltip in self.tooltips:
            tooltip.hide()
        self.tooltips = []
        for child in self.fields.winfo_children():
            child.destroy()
        form = self._form_of(step)
        self.form.configure(text=f"{step.title}: {step.about}")
        self.run_button.configure(text=f"Run {step.key}")
        row = 0

        def label(text):
            ttk.Label(self.fields, text=text).grid(
                row=row, column=0, sticky=tk.W, padx=(0, 8), pady=2
            )

        def tip(widget, text):
            # What each field is for, on hover: the explanation lines beside
            # the fields are gone (Graeme: almost every user would know).
            self.tooltips.append(Tooltip(widget, text))
            return widget

        label("Inputs")
        box = ttk.Frame(self.fields)
        box.grid(row=row, column=1, sticky=tk.EW, pady=2)
        # The button first, from the right: packed after the entry, it was the
        # one dropped where the column was narrower than the row asks.
        ttk.Button(
            box, text="Browse...", command=lambda: self._browse_inputs(step)
        ).pack(side=tk.RIGHT, padx=(4, 0))
        inputs = ttk.Entry(box, textvariable=form["inputs"], width=48)
        inputs.pack(side=tk.LEFT, fill=tk.X, expand=True)
        tip(
            inputs,
            (
                "The images' master files, a sweep each"
                if step.key == "import"
                else f"What mxi {step.key} reads: by default {' '.join(step.inputs)}"
            ),
        )
        row += 1
        for option in step.options:
            label(option.label)
            var = form["options"][option.flag]
            if option.kind == "flag":
                widget = ttk.Checkbutton(self.fields, variable=var)
            elif option.kind == "choice":
                widget = ttk.Combobox(
                    self.fields,
                    textvariable=var,
                    values=("",) + option.choices,
                    state="readonly",
                    width=20,
                )
            else:
                widget = ttk.Entry(self.fields, textvariable=var, width=24)
            widget.grid(row=row, column=1, sticky=tk.W, pady=2)
            tip(
                widget,
                f"{option.flag}: {option.about}" if option.about else option.flag,
            )
            row += 1
        if step.output_flag:
            label("Output")
            output = ttk.Entry(self.fields, textvariable=form["output"], width=24)
            output.grid(row=row, column=1, sticky=tk.W, pady=2)
            tip(
                output,
                f"{step.output_flag}: where to write"
                + (f" (by default {step.output})" if step.output else ""),
            )
            row += 1
        label("Threads")
        threads = ttk.Entry(self.fields, textvariable=form["threads"], width=8)
        threads.grid(row=row, column=1, sticky=tk.W, pady=2)
        tip(threads, "-j: threads to use; empty or 0, every core this job may use")
        row += 1
        label("Further options")
        extra = ttk.Entry(self.fields, textvariable=form["extra"], width=48)
        extra.grid(row=row, column=1, sticky=tk.EW, pady=2)
        tip(extra, f"Any other options, as mxi {step.key} --help lists them")
        self.fields.columnconfigure(1, weight=1)
        self._update_command()

    def _browse_inputs(self, step):
        if step.key == "import":
            chosen = filedialog.askopenfilenames(
                initialdir=self.directory.get(),
                title="The images: each sweep's master file",
                filetypes=[("NXmx master", "*.nxs *.h5"), ("All files", "*")],
            )
        else:
            chosen = filedialog.askopenfilenames(
                initialdir=self.directory.get(),
                title=f"What mxi {step.key} reads",
                filetypes=[("mxi's files", "*.rflx *.expt *.refl"), ("All files", "*")],
            )
        if chosen:
            here = self.directory.get()
            names = [
                os.path.relpath(c, here) if c.startswith(here) else c for c in chosen
            ]
            self._form_of(step)["inputs"].set(shlex.join(names))

    def form_argv(self, step=None):
        """The command the step's form makes, and what is wrong with it."""
        step = step or self.selected
        form = self._form_of(step)
        values = {flag: var.get() for flag, var in form["options"].items()}
        found = problems(step, values, form["threads"].get())
        try:
            inputs = shlex.split(form["inputs"].get())
        except ValueError as e:
            return None, [f"Inputs: {e}"]
        if not inputs:
            found.append(
                "Inputs: "
                + ("choose the images" if step.key == "import" else "none given")
            )
        try:
            argv = build(
                step,
                inputs,
                values,
                form["output"].get(),
                form["threads"].get(),
                form["extra"].get(),
            )
        except ValueError as e:
            return None, found + [f"Further options: {e}"]
        return argv, found

    def _update_command(self):
        if self.selected is None:
            return
        argv, found = self.form_argv()
        self.command.set(as_text(argv) if argv else "")
        self.trouble.set("; ".join(found))
        if not self.runner.running():
            self.run_button.configure(state=tk.DISABLED if found else tk.NORMAL)

    def copy_command(self):
        self.root.clipboard_clear()
        self.root.clipboard_append(self.command.get())

    def _text(self, parent):
        frame = ttk.Frame(parent)
        text = tk.Text(
            frame, wrap=tk.NONE, font="TkFixedFont", state=tk.DISABLED, undo=False
        )
        scroll = ttk.Scrollbar(frame, command=text.yview)
        text.configure(yscrollcommand=scroll.set)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        text.pack(fill=tk.BOTH, expand=True)
        return text

    # --- state ------------------------------------------------------------

    def refresh(self):
        """The steps' state from the directory: done where the output is."""
        done = state.scan(self.directory.get())
        for step in STEPS:
            self.tree.set(step.key, "state", DONE if done[step.key] else PENDING)
        if self.selected is None:
            first = state.next_step(self.directory.get()) or STEPS[-1].key
            self.tree.selection_set(first)
            self.tree.see(first)

    def choose_directory(self):
        chosen = filedialog.askdirectory(
            initialdir=self.directory.get(), mustexist=True
        )
        if chosen:
            self.directory.set(chosen)
            self.refresh()

    def _on_select(self, _event=None):
        chosen = self.tree.selection()
        if not chosen or chosen[0] in ("steps", "tools"):
            self.selected = None
            return
        self.selected = by_key(chosen[0])
        self._show_form(self.selected)
        self._show_summary(self.selected)

    def _show_summary(self, step):
        lines = [step.title, "", step.about, "", "  " + as_text(step.command())]
        if step.output:
            path = os.path.join(self.directory.get(), step.output)
            lines += [
                "",
                f"{step.output}: " + ("written" if os.path.isfile(path) else "not yet"),
            ]
        self._set(self.summary, "\n".join(lines) + "\n")

    # --- running ----------------------------------------------------------

    def run_selected(self):
        if self.selected is None:
            return
        argv, found = self.form_argv()
        if found or not argv:
            messagebox.showinfo("mxi ui", "\n".join(found) or "Nothing to run.")
            return
        self.run(argv, title=self.selected.title)

    def run(self, argv, title=None):
        """Runs a command in the working directory, its output in the log;
        argv[0] "mxi" is the mxi found (locate.py)."""
        if self.runner.running():
            messagebox.showinfo("mxi ui", "A command is running: wait, or cancel it.")
            return False
        if argv and argv[0] == "mxi":
            if not self.mxi:
                messagebox.showerror(
                    "mxi ui", "mxi is not on the PATH, nor beside mxi ui."
                )
                return False
            argv = [self.mxi, *argv[1:]]
        shown = as_text(["mxi", *argv[1:]] if argv[0] == self.mxi else argv)
        self._append(self.log, f"$ {shown}\n")
        self.tabs.select(self.log.master)
        try:
            self.runner.start(argv, cwd=self.directory.get())
        except OSError as e:
            self._append(self.log, f"cannot run it: {e}\n")
            return False
        self.status.set(f"Running {title or shown}...")
        self.cancel_button.configure(state=tk.NORMAL)
        self.run_button.configure(state=tk.DISABLED)
        return True

    def cancel(self):
        self.runner.cancel()

    def _poll(self):
        lines = self.runner.lines()
        if lines:
            self._append(self.log, "\n".join(lines) + "\n")
        code = self.runner.finished()
        if code is not None and self.cancel_button.instate(["!disabled"]):
            self._append(self.log, f"[exit {code}]\n\n")
            self.status.set("Done" if code == 0 else f"Failed, exit {code}")
            self.cancel_button.configure(state=tk.DISABLED)
            self.run_button.configure(state=tk.NORMAL)
            self._update_command()
            self.refresh()
            if self.selected is not None:
                self._show_summary(self.selected)
        self.root.after(POLL_MS, self._poll)

    # --- the rest ---------------------------------------------------------

    def _set(self, text, content):
        text.configure(state=tk.NORMAL)
        text.delete("1.0", tk.END)
        text.insert(tk.END, content)
        text.configure(state=tk.DISABLED)

    def _append(self, text, content):
        text.configure(state=tk.NORMAL)
        text.insert(tk.END, content)
        text.see(tk.END)
        text.configure(state=tk.DISABLED)

    def about(self):
        messagebox.showinfo(
            "About mxi ui",
            f"mxi ui {VERSION}\n\nAn interface to mxi, for processing data and for "
            f"teaching.\n\nmxi: {self.mxi or 'not found'}",
        )

    def quit(self):
        if self.runner.running():
            if not messagebox.askyesno(
                "mxi ui", "A command is running. Stop it and quit?"
            ):
                return
            self.runner.cancel()
        self.root.destroy()
