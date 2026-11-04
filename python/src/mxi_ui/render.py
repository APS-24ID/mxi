"""What the viewer draws, worked out without tkinter: a region of a frame
reduced to the window, its contrast and colour applied, as an image Tk takes;
where a detector pixel is on the canvas and back; and the resolution at a
pixel, from the experiment list's beam and detector.

A frame of an Eiger 16M is 4148 by 4362 pixels, a window a fifth of that:
each pixel shown is the maximum over the block it stands for, not their mean
or one of them, so that a spot one pixel wide shows at any zoom.
"""

import math

import numpy as np

COLOUR_MAPS = ("grey", "inverted", "heat")


def snap(scale):
    """A scale (detector pixels per screen pixel) made whole: a whole number
    of detector pixels per screen pixel zoomed out, of screen pixels per
    detector pixel zoomed in -- so that the reduced image is the canvas's own
    size, pixel for pixel, and the markers on it where they belong."""
    if scale >= 1:
        return float(max(1, int(math.ceil(scale - 1e-9))))
    return 1.0 / max(1, int(math.floor(1.0 / scale + 1e-9)))


def next_level(scale, inwards):
    """The whole scale next to this one (snap's): 3 to 2, 1 to 1/2 inwards."""
    if inwards:
        return scale - 1 if scale > 1 else 1.0 / (round(1.0 / scale) + 1)
    return (
        scale + 1
        if scale >= 1
        else (1.0 if round(1.0 / scale) <= 2 else 1.0 / (round(1.0 / scale) - 1))
    )


def reduce(frame, region, scale):
    """The region (x0, y0, x1, y1) of the frame, in detector pixels, at a
    scale snap() made: zoomed out, each pixel shown the maximum of the
    scale-by-scale block it stands for; zoomed in, each detector pixel
    repeated 1/scale times. Returns the array and the detector pixel (x, y)
    its top left stands for -- the region clipped to the detector."""
    height, width = frame.shape
    x0, y0, x1, y1 = region
    x0, y0 = max(0, int(math.floor(x0))), max(0, int(math.floor(y0)))
    x1, y1 = min(width, int(math.ceil(x1))), min(height, int(math.ceil(y1)))
    if x1 <= x0 or y1 <= y0:
        return np.zeros((1, 1), frame.dtype), (x0, y0)
    part = frame[y0:y1, x0:x1]
    if scale < 1:
        repeat = int(round(1.0 / scale))
        return np.repeat(np.repeat(part, repeat, axis=0), repeat, axis=1), (x0, y0)
    block = int(round(scale))
    if block == 1:
        return part, (x0, y0)
    # Rows first, whole rows at a time, then the columns of what is left, an
    # eighth the size at a block of eight; the partial blocks at the edges
    # apart. Padded and reduced in one, as at first, a 16M frame took 160 ms
    # at the fit, this 30: the same maxima.
    h, w = part.shape
    hb, wb = h // block, w // block
    rows = part[: hb * block].reshape(hb, block, w).max(axis=1)
    if h % block:
        rows = np.vstack([rows, part[hb * block :].max(axis=0, keepdims=True)])
    out = rows[:, : wb * block].reshape(rows.shape[0], wb, block).max(axis=2)
    if w % block:
        out = np.hstack([out, rows[:, wb * block :].max(axis=1, keepdims=True)])
    return out, (x0, y0)


def colour(values, top, scheme="grey"):
    """Counts to RGB bytes, 0 to `top` the scale: "grey" dark spots on a light
    ground, as diffraction is drawn; "inverted" the other way; "heat" black to
    red to yellow to white."""
    t = np.clip(values.astype(np.float32) / max(float(top), 1e-9), 0.0, 1.0)
    if scheme == "inverted":
        g = (t * 255).astype(np.uint8)
        return np.dstack([g, g, g])
    if scheme == "heat":
        r = np.clip(3 * t, 0, 1)
        g = np.clip(3 * t - 1, 0, 1)
        b = np.clip(3 * t - 2, 0, 1)
        return (np.dstack([r, g, b]) * 255).astype(np.uint8)
    g = ((1.0 - t) * 255).astype(np.uint8)
    return np.dstack([g, g, g])


def ppm(rgb):
    """RGB bytes as a PPM image, which Tk's PhotoImage reads with no other
    library."""
    h, w = rgb.shape[:2]
    return b"P6\n%d %d\n255\n" % (w, h) + np.ascontiguousarray(rgb, np.uint8).tobytes()


class View:
    """Where the detector is on the canvas: a scale (detector pixels per
    canvas pixel) and the detector point at the canvas's top left. Zoom keeps
    the point under the cursor under it; pan moves by canvas pixels."""

    def __init__(self, detector_size, canvas_size):
        self.width, self.height = detector_size
        self.fit(canvas_size)

    def fit(self, canvas_size):
        cw, ch = max(1, canvas_size[0]), max(1, canvas_size[1])
        self.scale = snap(max(self.width / cw, self.height / ch))
        self.x0 = (self.width - cw * self.scale) / 2
        self.y0 = (self.height - ch * self.scale) / 2

    def to_canvas(self, x, y):
        return (x - self.x0) / self.scale, (y - self.y0) / self.scale

    def to_detector(self, cx, cy):
        return self.x0 + cx * self.scale, self.y0 + cy * self.scale

    def zoom(self, factor, cx, cy, smallest=1 / 16, largest=None):
        """In by `factor` > 1, out below, about canvas point (cx, cy), to the
        next whole scale (snap) at least: a wheel's step always moves it."""
        x, y = self.to_detector(cx, cy)
        wanted = self.scale / factor
        if largest is not None:
            wanted = min(wanted, largest)
        scale = snap(max(wanted, smallest))
        if scale == self.scale:
            # A step too small to change the whole scale: the next level --
            # ..., 1/3, 1/2, 1, 2, 3, ... -- the way the wheel turned.
            scale = next_level(self.scale, inwards=factor > 1)
            scale = max(min(scale, largest or scale), smallest)
        self.scale = scale
        self.x0, self.y0 = x - cx * self.scale, y - cy * self.scale

    def pan(self, dx, dy):
        self.x0 -= dx * self.scale
        self.y0 -= dy * self.scale

    def region(self, canvas_size):
        cw, ch = canvas_size
        return self.x0, self.y0, self.x0 + cw * self.scale, self.y0 + ch * self.scale


def resolution(document, experiment, x, y):
    """The resolution, in Å, at detector pixel (x, y) of the experiment's
    first panel -- from its beam (direction towards the source, as DIALS's)
    and detector -- or None where the geometry is missing or the pixel is on
    the beam."""
    try:
        e = document["experiment"][experiment]
        beam = document["beam"][e["beam"]]
        panel = document["detector"][e["detector"]]["panels"][0]
        fast = np.asarray(panel["fast_axis"], float)
        slow = np.asarray(panel["slow_axis"], float)
        origin = np.asarray(panel["origin"], float)
        px, py = panel["pixel_size"]
        towards_source = np.asarray(beam["direction"], float)
        wavelength = float(beam["wavelength"])
    except (KeyError, IndexError, TypeError, ValueError):
        return None
    point = origin + x * px * fast + y * py * slow
    n = np.linalg.norm(point)
    m = np.linalg.norm(towards_source)
    if n == 0 or m == 0:
        return None
    cos_two_theta = float(np.dot(point / n, -towards_source / m))
    two_theta = math.acos(max(-1.0, min(1.0, cos_two_theta)))
    if two_theta < 1e-6:
        return None
    return wavelength / (2 * math.sin(two_theta / 2))
