"""The main window: the steps down the side, their images, log, summary and
plots in tabs, a working directory at the top, a status line at the foot.

Phase 1 of the interface (docs/ui.md): the shell, and running a step with its
default command, streamed into the log. The viewer, each step's options and
the plots fill the tabs in the phases after.
"""

import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from mxi_ui import state
from mxi_ui.locate import find_mxi
from mxi_ui.runner import Runner, as_text
from mxi_ui.steps import STEPS, TOOLS, by_key
from mxi_ui.version import VERSION

POLL_MS = 100
DONE, PENDING = "done", ""


class MainWindow:
    def __init__(self, root, directory=None):
        self.root = root
        self.runner = Runner()
        self.mxi = find_mxi()
        self.directory = tk.StringVar(value=os.path.abspath(directory or os.getcwd()))
        self.status = tk.StringVar(value="Ready")
        self.selected = None
        root.title(f"mxi ui {VERSION}")
        root.geometry("1100x720")
        root.minsize(800, 500)
        self._menus()
        self._layout()
        self.refresh()
        root.after(POLL_MS, self._poll)
        root.protocol("WM_DELETE_WINDOW", self.quit)

    # --- building ---------------------------------------------------------

    def _menus(self):
        bar = tk.Menu(self.root)
        file = tk.Menu(bar, tearoff=False)
        file.add_command(label="Open directory...", command=self.choose_directory)
        file.add_separator()
        file.add_command(label="Quit", command=self.quit)
        bar.add_cascade(label="File", menu=file)
        run = tk.Menu(bar, tearoff=False)
        run.add_command(label="Run the selected step", command=self.run_selected)
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
        self.tree = ttk.Treeview(
            side, columns=("state",), show="tree headings", selectmode="browse"
        )
        self.tree.heading("#0", text="Step")
        self.tree.heading("state", text="")
        self.tree.column("#0", width=170)
        self.tree.column("state", width=60, anchor=tk.CENTER)
        self.tree.insert("", tk.END, iid="steps", text="Processing", open=True)
        for step in STEPS:
            self.tree.insert("steps", tk.END, iid=step.key, text=step.title)
        self.tree.insert("", tk.END, iid="tools", text="Tools", open=True)
        for step in TOOLS:
            self.tree.insert("tools", tk.END, iid=step.key, text=step.title)
        self.tree.pack(fill=tk.BOTH, expand=True)
        self.tree.bind("<<TreeviewSelect>>", self._on_select)
        self.run_button = ttk.Button(side, text="Run", command=self.run_selected)
        self.run_button.pack(fill=tk.X, pady=6)
        panes.add(side, weight=0)

        self.tabs = ttk.Notebook(panes)
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
        panes.add(self.tabs, weight=1)

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
        step = self.selected
        inputs = None
        if step.key == "import":
            chosen = filedialog.askopenfilenames(
                initialdir=self.directory.get(),
                title="The images: each sweep's master file",
                filetypes=[("NXmx master", "*.nxs *.h5"), ("All files", "*")],
            )
            if not chosen:
                return
            inputs = list(chosen)
        self.run(step.command(inputs), title=step.title)

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
