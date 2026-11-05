"""Scrolling, however it comes: a wheel's notches, X11's buttons 4 and 5, or a
trackpad's two fingers -- Tk 9's own event, <TouchpadScroll>, which a Mac's
trackpad sends there and no <MouseWheel> (Homebrew's Python has Tk 9). Bound
for each, so that nothing need have a wheel (Graeme has none)."""

import tkinter as tk


def touchpad_delta(widget, delta):
    """A <TouchpadScroll>'s (dx, dy) in pixels: tk::PreciseScrollDeltas where
    Tk has it, else the two 16-bit halves of %D, signed."""
    try:
        dx, dy = (
            int(v)
            for v in widget.tk.splitlist(
                widget.tk.call("tk::PreciseScrollDeltas", delta)
            )
        )
        return dx, dy
    except (tk.TclError, ValueError):
        dx, dy = (delta >> 16) & 0xFFFF, delta & 0xFFFF
        return (dx - 0x10000 if dx >= 0x8000 else dx), (
            dy - 0x10000 if dy >= 0x8000 else dy
        )


def wheel_pixels(delta):
    """A <MouseWheel>'s delta as pixels to scroll the view by, up negative:
    multiples of 120 a notch on Windows, a few units at a time on a Mac."""
    notches = delta / 120.0 if abs(delta) >= 120 else float(delta)
    return -notches * 40.0


def bind_scrolling(widgets, scroll):
    """Each widget scrolls by `scroll(pixels)` -- positive down -- from a
    wheel, X11's buttons, or a trackpad; each binding replacing the widget's
    own (matplotlib's canvas binds the wheel)."""
    for w in widgets:
        w.bind("<MouseWheel>", lambda e: scroll(wheel_pixels(e.delta)))
        w.bind("<Button-4>", lambda _e: scroll(-40))
        w.bind("<Button-5>", lambda _e: scroll(40))
        try:
            w.bind(
                "<TouchpadScroll>",
                lambda e, w=w: scroll(-touchpad_delta(w, e.delta)[1]),
            )
        except tk.TclError:  # Tk 8.6 knows no such event
            pass
