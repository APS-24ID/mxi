"""The main window: the steps down the side; for the chosen step, tabs --
Setup, its form (inputs, main options, output, threads, further options),
then its log, summary, images and plots -- with the command the form makes,
Copy and Run in a bar beneath them, in view whichever tab is shown; a working
directory at the top, a status line at the foot. The viewer and the plots
fill their tabs in the phases after (docs/ui.md).
"""

import os
import shlex
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter import font as tkfont

from mxi_ui import results, state
from mxi_ui.locate import find_mxi
from mxi_ui.model import build, problems
from mxi_ui.runner import Runner, as_text
from mxi_ui.steps import STEPS, TOOLS, by_key
from mxi_ui.tooltip import Tooltip
from mxi_ui.version import VERSION


def draw(axes, plot):
    """One of results.plots on matplotlib's axes: a resolution axis spaced in
    1/d^2 and ticked in d, a log axis for R factors, a second axis on the
    right where the plot has one."""
    if plot.kind == "matrix":
        image = axes.imshow(
            plot.matrix,
            vmin=plot.low,
            vmax=1.0,
            cmap="viridis",
            interpolation="nearest",
        )
        axes.figure.colorbar(image, ax=axes, shrink=0.8)
    for series in plot.series:
        pairs = [
            (x, y)
            for x, y in zip(series.x, series.y)
            if x is not None and y is not None
        ]
        if not pairs:
            continue
        xs, ys = zip(*pairs)
        if plot.kind == "points":
            axes.scatter(xs, ys, s=14, label=series.name)
        elif plot.bars:
            if not plot.resolution:
                width = 0.9  # a bar at each whole position: data sets in order
            else:
                width = (
                    (max(xs) - min(xs)) / max(len(xs), 1) * 0.9 if len(xs) > 1 else 0.01
                )
            axes.bar(xs, ys, width=width, label=series.name)
        else:
            axes.plot(
                xs,
                ys,
                marker="o" if len(xs) <= 30 else None,
                markersize=3,
                label=series.name,
            )
    if plot.log:
        axes.set_yscale("log")
    if plot.integer:
        from matplotlib.ticker import MaxNLocator

        axes.xaxis.set_major_locator(MaxNLocator(integer=True))
    axes.set_title(plot.title, fontsize=10)
    axes.set_xlabel(plot.xlabel)
    axes.set_ylabel(plot.ylabel)
    if plot.resolution:
        xs = [x for s in plot.series + plot.right for x in s.x]
        positions, labels = results.d_ticks(xs)
        axes.set_xticks(positions, labels)
    if plot.right:
        twin = axes.twinx()
        for series in plot.right:
            pairs = [
                (x, y)
                for x, y in zip(series.x, series.y)
                if x is not None and y is not None
            ]
            if pairs:
                xs, ys = zip(*pairs)
                twin.plot(
                    xs,
                    ys,
                    color="tab:orange",
                    marker="o",
                    markersize=3,
                    label=series.name,
                )
        twin.set_ylabel(plot.right_label)
        twin.set_ylim(bottom=0)
    if len(plot.series) + len(plot.right) > 1:
        # One legend for both axes, the right's series among the left's.
        handles, labels = axes.get_legend_handles_labels()
        if plot.right:
            more = twin.get_legend_handles_labels()
            handles, labels = handles + more[0], labels + more[1]
        axes.legend(handles, labels, fontsize=8)
    if plot.kind != "matrix":  # no grid over a matrix's cells
        axes.grid(alpha=0.3)


POLL_MS = 100
ROW_HEIGHT = 330  # pixels a row of plots, at least: more scroll
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
        self.viewers = []
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
        # The images a window of their own, for the step chosen: its file's
        # spots and predictions on them (viewer.py), as the DIALS interface
        # opened its viewer for a step.
        self.view_button = ttk.Button(
            side, text="View images", command=self.view_images
        )
        self.view_button.pack(fill=tk.X, pady=(6, 0))
        Tooltip(
            self.view_button, "The chosen step's images, with its spots and predictions"
        )
        panes.add(side, weight=0)

        # The step's form a tab, Setup, beside its log, summary, images and
        # plots, with the whole height; the command it makes, Copy and Run in a
        # bar beneath, in view whichever tab is shown. Above the tabs, in a
        # divided pane, the form had been squeezed until its last row -- the
        # command and Run -- was a sliver (Graeme's screenshot, on a Mac).
        right = ttk.Frame(panes)
        self._action_bar(right)
        self.tabs = ttk.Notebook(right)
        self.tabs.pack(fill=tk.BOTH, expand=True)
        self.setup = ttk.Frame(self.tabs, padding=(12, 10))
        self.heading = ttk.Label(self.setup, text="", font="TkHeadingFont")
        self.heading.pack(anchor=tk.W, pady=(0, 10))
        self.fields = ttk.Frame(self.setup)
        self.fields.pack(fill=tk.X)
        # The step's notes beneath its form, for learning what it does
        # (notes.py): what it does, what to look at, what goes wrong, the
        # guides -- in the space left, scrolling, so that the form keeps its
        # place however long they are.
        holder = ttk.Frame(self.setup, padding=(0, 12, 0, 0))
        holder.pack(fill=tk.BOTH, expand=True)
        self.notes = tk.Text(
            holder,
            wrap=tk.WORD,
            height=8,
            relief=tk.FLAT,
            borderwidth=0,
            highlightthickness=0,
            state=tk.DISABLED,
            cursor="arrow",
            font="TkDefaultFont",
            background=self.root.cget("background"),
        )
        bar = ttk.Scrollbar(holder, command=self.notes.yview)
        self.notes.configure(yscrollcommand=bar.set)
        bar.pack(side=tk.RIGHT, fill=tk.Y)
        self.notes.pack(fill=tk.BOTH, expand=True)
        self.notes.tag_configure("heading", font="TkHeadingFont", spacing1=6)
        self.notes.tag_configure("link", foreground="#1e88e5", underline=True)
        self.notes.tag_bind(
            "link", "<Enter>", lambda _e: self.notes.configure(cursor="hand2")
        )
        self.notes.tag_bind(
            "link", "<Leave>", lambda _e: self.notes.configure(cursor="arrow")
        )
        self.log = self._text(self.tabs)
        self.summary = self._text(self.tabs)
        self.plots = ttk.Frame(self.tabs, padding=6)
        self.canvas = None
        self.toolbar = None
        # The tab's size governs the plots, not the plots' wish: the canvas
        # is made to fill it whenever it changes, and matplotlib draws the
        # figure to the canvas (_fit_plots).
        self.plots.pack_propagate(False)
        self.plots.bind("<Configure>", lambda _e: self._fit_plots())
        # And when it is shown: built while hidden, a tab's frame reads a
        # pixel wide, and showing it changes nothing to raise <Configure>.
        self.tabs.bind(
            "<<NotebookTabChanged>>",
            lambda _e: (
                self._fit_plots() if self.tabs.select() == str(self.plots) else None
            ),
        )
        self.tabs.add(self.setup, text="Setup")
        self.tabs.add(self.log.master, text="Log")
        self.tabs.add(self.summary.master, text="Summary")
        self.tabs.add(self.plots, text="Plots")
        panes.add(right, weight=1)

    def _action_bar(self, parent):
        """The command the form makes, its problems, Copy and Run, beneath the
        tabs. Run and Copy packed first, from the right, so that a narrow
        window narrows the command, never the buttons."""
        bar = ttk.Frame(parent, padding=(0, 6, 0, 2))
        bar.pack(side=tk.BOTTOM, fill=tk.X)
        row = ttk.Frame(bar)
        row.pack(fill=tk.X)
        self.run_button = ttk.Button(row, text="Run", command=self.run_selected)
        self.run_button.pack(side=tk.RIGHT, padx=(6, 0))
        ttk.Button(row, text="Copy", command=self.copy_command).pack(
            side=tk.RIGHT, padx=(6, 0)
        )
        shown = ttk.Entry(
            row, textvariable=self.command, font="TkFixedFont", state="readonly"
        )
        shown.pack(side=tk.LEFT, fill=tk.X, expand=True)
        # Bright enough to read on a dark window (a Mac's dark mode) and a light.
        ttk.Label(bar, textvariable=self.trouble, foreground="#e53935").pack(
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
        self.heading.configure(text=f"{step.title}: {step.about}")
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

    def images_file(self, step=None):
        """The file to view for a step: its output where written, else what it
        reads -- after Import, the images alone; None where neither is there."""
        step = step or self.selected
        if step is None:
            return None
        here = self.directory.get()
        for name in ([step.output] if step.output.endswith(".rflx") else []) + list(
            step.inputs
        ):
            path = os.path.join(here, name)
            if name.endswith((".rflx", ".expt")) and os.path.isfile(path):
                return path
        return None

    def view_images(self):
        path = self.images_file()
        if path is None:
            messagebox.showinfo(
                "mxi ui",
                (
                    "Nothing to view yet: run Import first."
                    if self.selected
                    else "Choose a step."
                ),
            )
            return None
        from mxi_ui.viewer import Viewer

        try:
            viewer = Viewer(
                self.root,
                path,
                f"mxi ui: {self.selected.title}, {os.path.basename(path)}",
                on_saved=self.masked,
            )
        except Exception as e:  # noqa: BLE001 -- said, whatever it is
            messagebox.showerror(
                "mxi ui", f"Cannot show the images of {os.path.basename(path)}:\n{e}"
            )
            return None
        self.viewers.append(viewer)
        return viewer

    def masked(self, path):
        """The viewer saved untrusted regions to `path` (masked.rflx): spot
        finding reads it from now on, as docs/mask.md has it -- the form's
        input set, the step to run next said."""
        name = os.path.basename(path)
        form = self._form_of(by_key("find"))
        form["inputs"].set(name)
        self.status.set(f"Untrusted regions saved in {name}: Find spots reads it")
        self._update_command()

    def copy_command(self):
        self.root.clipboard_clear()
        self.root.clipboard_append(self.command.get())

    def _show_notes(self, step):
        """The step's notes into the Setup tab's text, its guides as links
        that open in the browser."""
        from mxi_ui import notes

        note = notes.note_for(step.key)
        t = self.notes
        t.configure(state=tk.NORMAL)
        t.delete("1.0", tk.END)
        if note is not None:
            for heading, body in (
                ("What it does", note.does),
                ("What to look at", note.look),
                ("If it goes wrong", note.wrong),
            ):
                t.insert(tk.END, heading + "\n", "heading")
                t.insert(tk.END, body + "\n")
            t.insert(tk.END, "Read more\n", "heading")
            for k, (label, path) in enumerate(note.read):
                tag = f"link{k}"
                t.insert(tk.END, label, ("link", tag))
                t.tag_bind(
                    tag,
                    "<Button-1>",
                    lambda _e, url=notes.link(path): self.open_link(url),
                )
                t.insert(tk.END, "    ")
        t.configure(state=tk.DISABLED)

    def open_link(self, url):
        import webbrowser

        webbrowser.open(url)
        self.status.set(f"Opened {url}")

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
        self._show_notes(self.selected)
        self._show_summary(self.selected)

    def _show_summary(self, step):
        """The step, its command, and what its file records: its facts and
        tables (results.py), read from the file, not the log."""
        lines = [step.title, "", step.about, "", "  " + as_text(step.command())]
        try:
            path, rec = results.record(self.directory.get(), step.key)
        except Exception as e:  # noqa: BLE001 -- said, whatever it is
            path, rec = None, None
            lines += ["", f"{step.output} cannot be read: {e}"]
        if step.output:
            there = os.path.isfile(os.path.join(self.directory.get(), step.output))
            lines += ["", f"{step.output}: " + ("written" if there else "not yet")]
        if rec:
            lines += ["", f"Run {rec.get('time', '')}".rstrip(), ""]
            shown = results.facts(rec)
            width = max((len(label) for label, _ in shown), default=0)
            lines += [f"  {label:<{width}}  {text}" for label, text in shown]
            for table in results.tables(rec):
                lines += ["", table.title, ""] + self._table_lines(table)
        elif path and os.path.isfile(path):
            # Which steps it does record, so that a file written some other
            # way -- converted, say, which carries the steps before it and
            # adds none -- says so.
            try:
                held = [s.get("program", "?") for s in results.steps_of(path)]
            except Exception:  # noqa: BLE001 -- said above, if it is the same
                held = []
            lines += [
                "",
                f"{step.output} records no step of mxi {step.key}: "
                + (
                    f"it records {', '.join(held)}."
                    if held
                    else "it records no steps at all."
                ),
            ]
        self._set(self.summary, "\n".join(lines) + "\n")
        self._show_plots(step, rec)

    def _table_lines(self, table):
        widths = [
            (
                max(len(str(c)), *(len(str(r[i])) for r in table.rows))
                if table.rows
                else len(str(c))
            )
            for i, c in enumerate(table.columns)
        ]

        def line(cells):
            return "  " + "  ".join(
                (str(v).ljust(w) if i == 0 else str(v).rjust(w))
                for i, (v, w) in enumerate(zip(cells, widths))
            )

        return [line(table.columns)] + [line(r) for r in table.rows]

    def _show_plots(self, step, rec):
        """The step's plots, drawn with matplotlib: two to a row, its toolbar
        to zoom, pan and save."""
        for child in self.plots.winfo_children():
            child.destroy()
        self.canvas = None
        self.toolbar = None
        self.scroller = None
        made = results.plots(rec)
        if not made:
            text = (
                f"No plots for {step.title.lower()}."
                if rec
                else f"Plots appear once {step.title.lower()} has run."
            )
            ttk.Label(self.plots, text=text).pack(anchor=tk.NW, padx=6, pady=6)
            return
        from matplotlib.backends.backend_tkagg import (
            FigureCanvasTkAgg,
            NavigationToolbar2Tk,
        )
        from matplotlib.figure import Figure

        columns = 1 if len(made) == 1 else 2
        rows = (len(made) + columns - 1) // columns
        # Any size to begin: _fit_plots draws it to the tab's.
        figure = Figure(figsize=(8, 5), layout="constrained")
        for i, plot in enumerate(made):
            draw(figure.add_subplot(rows, columns, i + 1), plot)
        self.figure = figure
        self.plot_rows = rows
        # The figure in a canvas that scrolls: each row of plots its own
        # height, the tab scrolling where they are more than it holds --
        # five or six plots squeezed into it had been too small to read.
        holder = ttk.Frame(self.plots)
        self.scroller = tk.Canvas(holder, highlightthickness=0, borderwidth=0)
        bar = ttk.Scrollbar(holder, command=self.scroller.yview)
        self.scroller.configure(yscrollcommand=bar.set)
        self.canvas = FigureCanvasTkAgg(figure, master=self.scroller)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.plots, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        holder.pack(fill=tk.BOTH, expand=True)
        bar.pack(side=tk.RIGHT, fill=tk.Y)
        self.scroller.pack(fill=tk.BOTH, expand=True)
        self.scroller.create_window(
            0, 0, window=self.canvas.get_tk_widget(), anchor=tk.NW
        )
        widget = self.canvas.get_tk_widget()
        for target in (widget, self.scroller):
            target.bind(
                "<MouseWheel>", lambda e: self._scroll_plots(-1 if e.delta > 0 else 1)
            )
            target.bind("<Button-4>", lambda _e: self._scroll_plots(-1))
            target.bind("<Button-5>", lambda _e: self._scroll_plots(1))
        self._fit_plots()

    def _fit_plots(self):
        """The figure the tab's width, each row of plots ROW_HEIGHT high or the
        tab's height shared, whichever is more; scrolling where it is taller
        than the tab. Sized from its own wish, the figure had been drawn
        wider than the tab and its right-hand plots cut off."""
        if self.canvas is None:
            return
        self.plots.update_idletasks()
        width = self.scroller.winfo_width()
        shown = self.scroller.winfo_height()
        if width < 50 or shown < 50:
            return
        height = max(shown, self.plot_rows * ROW_HEIGHT)
        self.canvas.get_tk_widget().configure(width=width, height=height)
        # And the figure itself: matplotlib redraws to a resize event, and
        # where the widget's size did not change there is none.
        dpi = self.figure.dpi
        self.figure.set_size_inches(width / dpi, height / dpi, forward=False)
        self.scroller.configure(scrollregion=(0, 0, width, height))
        self.canvas.draw_idle()

    def _scroll_plots(self, units):
        if self.canvas is not None:
            self.scroller.yview_scroll(units, "units")

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
