"""The background figure for mxeq background: four panels.

1. Every reflection's background against 1/d^2, as a density, with the running
   median and its 5th and 95th percentiles: the smooth curve, if it is one.
2. Each reflection's background over that median, against the azimuth around
   the beam: polarisation shows as a slow wave, the backstop and its flare as
   narrow features.
3. The same against the image: a drift through the scan.
4. The same on the detector, near the beam: the backstop shadow dark, its flare
   bright.

Panels 2 and 3 can be limited to a resolution range, to see the innermost apart
from the rest.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

from .. import background as bg  # noqa: E402


def draw(
    b: bg.Backgrounds,
    path: str,
    d_range: tuple[float, float] | None = None,
    zoom: float = 250.0,
    title: str = "",
) -> None:
    c = bg.curve(b)
    r = bg.ratios(b, c)
    s2 = 1.0 / b.d**2
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    if title:
        fig.suptitle(title)

    ax = axes[0, 0]
    positive = b.background > 0
    hb = ax.hexbin(
        s2[positive],
        b.background[positive],
        gridsize=120,
        yscale="log",
        bins="log",
        cmap="Greys",
        mincnt=1,
    )
    fig.colorbar(hb, ax=ax, label="reflections (log)")
    ax.plot(c.s2, c.median, color="tab:red", lw=1.5, label="median")
    ax.plot(
        c.s2, c.low, color="tab:red", lw=0.8, ls="--", label="5th and 95th percentiles"
    )
    ax.plot(c.s2, c.high, color="tab:red", lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_ylabel("background (counts per pixel)")
    ticks = ax.get_xticks()
    ticks = ticks[(ticks > 0) & (ticks <= s2.max())]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{1 / np.sqrt(t):.2f}" for t in ticks])
    ax.set_xlabel("resolution (A), on a 1/d^2 scale")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("background against resolution")

    sel = np.isfinite(r) & (r > 0)
    label = "all resolutions"
    if d_range is not None:
        hi, lo = max(d_range), min(d_range)
        sel &= (b.d <= hi) & (b.d >= lo)
        label = f"{hi:g} to {lo:g} A"
    rr = np.clip(r, 1e-2, 1e2)
    for ax, x, name in (
        (axes[0, 1], b.azimuth, "azimuth around the beam (degrees)"),
        (axes[1, 0], b.z, "image"),
    ):
        if sel.any():
            hb = ax.hexbin(
                x[sel],
                rr[sel],
                gridsize=(120, 60),
                yscale="log",
                bins="log",
                cmap="Greys",
                mincnt=1,
            )
            fig.colorbar(hb, ax=ax, label="reflections (log)")
        ax.axhline(1.0, color="tab:red", lw=1)
        ax.axhline(0.5, color="tab:blue", lw=0.8, ls="--")
        ax.axhline(2.0, color="tab:orange", lw=0.8, ls="--")
        ax.set_yscale("log")
        ax.set_xlabel(name)
        ax.set_ylabel("background over the median at its resolution")
        ax.set_title(f"against {name.split(' (')[0]}, {label}")

    ax = axes[1, 1]
    near = (b.radius <= zoom) & np.isfinite(r) & (r > 0)
    if near.any():
        sc = ax.scatter(
            b.x[near] - b.centre[0],
            b.y[near] - b.centre[1],
            c=np.clip(r[near], 0.05, 20),
            s=6,
            cmap="coolwarm",
            norm=LogNorm(vmin=0.05, vmax=20),
        )
        fig.colorbar(sc, ax=ax, label="background over the median at its resolution")
    ax.plot([0], [0], "+", color="k", ms=12)
    ax.set_xlim(-zoom, zoom)
    ax.set_ylim(zoom, -zoom)  # slow increasing downwards, as an image viewer draws it
    ax.set_aspect("equal")
    ax.set_xlabel("fast from the beam centre (pixels)")
    ax.set_ylabel("slow from the beam centre (pixels)")
    ax.set_title(f"on the detector, within {zoom:g} pixels of the beam")

    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
