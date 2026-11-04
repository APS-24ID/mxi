"""A tooltip: a line of text in a borderless window beside a widget, shown
after the pointer rests on it a moment, gone when it leaves or clicks."""

import tkinter as tk

DELAY_MS = 500


class Tooltip:
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self._pending = None
        self._window = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self.hide, add="+")
        widget.bind("<ButtonPress>", self.hide, add="+")

    def _schedule(self, _event=None):
        self._cancel()
        self._pending = self.widget.after(DELAY_MS, self.show)

    def _cancel(self):
        if self._pending is not None:
            self.widget.after_cancel(self._pending)
            self._pending = None

    def show(self):
        self._pending = None
        if self._window is not None or not self.text:
            return
        x = self.widget.winfo_rootx() + 12
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        self._window = window = tk.Toplevel(self.widget)
        window.wm_overrideredirect(True)
        window.wm_geometry(f"+{x}+{y}")
        tk.Label(
            window,
            text=self.text,
            justify=tk.LEFT,
            background="#ffffe0",
            foreground="#000000",
            relief=tk.SOLID,
            borderwidth=1,
            padx=6,
            pady=3,
            wraplength=420,
        ).pack()

    def hide(self, _event=None):
        self._cancel()
        if self._window is not None:
            self._window.destroy()
            self._window = None

    def showing(self):
        return self._window is not None
