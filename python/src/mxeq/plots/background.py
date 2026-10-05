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


def draw_model(
    fit,
    bg,
    x_px,
    y_px,
    centre,
    path: str,
    z_max: float = 5.0,
    zoom: float = 250.0,
    title: str = "",
) -> None:
    """The background model against the data: R(s) through the backgrounds
    with G, P, Omega and Q divided out; G(phi) through them with R, P, Omega
    and Q divided out; z against resolution; z on the detector by the beam; and
    where the reflections were removed and kept, over the whole detector and by
    the beam."""
    fig, axes = plt.subplots(3, 2, figsize=(14, 15))
    if title:
        fig.suptitle(title)
    ok = fit.used & (bg > 0)

    ax = axes[0, 0]
    y = bg[ok] / (fit.g_row[ok] * fit.known[ok])
    hb = ax.hexbin(
        fit.s[ok], y, gridsize=120, yscale="log", bins="log", cmap="Greys", mincnt=1
    )
    fig.colorbar(hb, ax=ax, label="reflections (log)")
    ax.plot(
        fit.s_curve[0], fit.s_curve[1], color="tab:red", lw=1.5, label="R(s), fitted"
    )
    ax.set_yscale("log")
    ax.set_xlabel("s = 1/d (1/A)")
    ax.set_ylabel("background over G . P . Omega . Q")
    ax.legend(fontsize=8)
    ax.set_title("the radial scatter, R")

    ax = axes[0, 1]
    y = bg[ok] / (fit.r_row[ok] * fit.known[ok])
    hb = ax.hexbin(
        fit.phi[ok],
        y,
        gridsize=(120, 60),
        yscale="log",
        bins="log",
        cmap="Greys",
        mincnt=1,
    )
    fig.colorbar(hb, ax=ax, label="reflections (log)")
    ax.plot(
        fit.phi_curve[0],
        fit.phi_curve[1],
        color="tab:red",
        lw=1.5,
        label="G(phi), fitted",
    )
    ax.set_yscale("log")
    ax.set_ylim(0.3, 3.0)
    ax.set_xlabel("rotation (degrees)")
    ax.set_ylabel("background over R . P . Omega . Q")
    ax.legend(fontsize=8)
    ax.set_title("the rotation, G")

    ax = axes[1, 0]
    finite = np.isfinite(fit.z)
    hb = ax.hexbin(
        fit.s[finite],
        np.clip(fit.z[finite], -50, 50),
        gridsize=(120, 80),
        bins="log",
        cmap="Greys",
        mincnt=1,
    )
    fig.colorbar(hb, ax=ax, label="reflections (log)")
    for v in (-z_max, z_max):
        ax.axhline(v, color="tab:red", lw=0.8, ls="--")
    ax.set_xlabel("s = 1/d (1/A)")
    ax.set_ylabel("z, clipped at 50")
    ax.set_title(f"z against resolution, |z| beyond {z_max:g} flagged")

    ax = axes[1, 1]
    radius = np.hypot(x_px - centre[0], y_px - centre[1])
    near = (radius <= zoom) & (np.isfinite(fit.z) | np.isneginf(fit.z))
    zn = np.clip(np.where(np.isneginf(fit.z[near]), -50, fit.z[near]), -20, 20)
    sc = ax.scatter(
        x_px[near] - centre[0],
        y_px[near] - centre[1],
        c=zn,
        s=6,
        cmap="coolwarm",
        vmin=-20,
        vmax=20,
    )
    fig.colorbar(sc, ax=ax, label="z, clipped at 20")
    ax.plot([0], [0], "+", color="k", ms=12)
    ax.set_xlim(-zoom, zoom)
    ax.set_ylim(zoom, -zoom)
    ax.set_aspect("equal")
    ax.set_xlabel("fast from the beam centre (pixels)")
    ax.set_ylabel("slow from the beam centre (pixels)")
    ax.set_title(f"z on the detector, within {zoom:g} pixels of the beam")

    # Where reflections were removed and kept: the kept as a grey density --
    # too many to draw one by one -- and the removed over them, by reason.
    low = fit.flagged_low
    high = fit.flagged_high
    shell = fit.flagged_shell if len(fit.flagged_shell) else np.zeros(len(x_px), bool)
    removed = low | high | shell
    have = np.isfinite(x_px) & np.isfinite(y_px)
    kept = have & ~removed
    for ax, window, name in (
        (axes[2, 0], None, "over the detector"),
        (axes[2, 1], zoom, f"within {zoom:g} pixels of the beam"),
    ):
        sel = (
            kept
            if window is None
            else kept & (np.hypot(x_px - centre[0], y_px - centre[1]) <= window)
        )
        if sel.any():
            ax.hexbin(
                x_px[sel] - centre[0],
                y_px[sel] - centre[1],
                gridsize=150 if window is None else 80,
                bins="log",
                cmap="Greys",
                mincnt=1,
                alpha=0.6,
            )
        for mask, colour, label in (
            (shell, "tab:orange", "a whole shell left out"),
            (high, "tab:red", "high: the flare's kind"),
            (low, "tab:blue", "low: the shadow's kind"),
        ):
            m = mask & have
            if window is not None:
                m &= np.hypot(x_px - centre[0], y_px - centre[1]) <= window
            ax.scatter(
                x_px[m] - centre[0],
                y_px[m] - centre[1],
                s=4 if window is None else 10,
                c=colour,
                label=f"{label}, {int(m.sum())}",
                linewidths=0,
            )
        ax.plot([0], [0], "+", color="k", ms=12)
        if window is not None:
            ax.set_xlim(-window, window)
            ax.set_ylim(window, -window)
        else:
            ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.set_xlabel("fast from the beam centre (pixels)")
        ax.set_ylabel("slow from the beam centre (pixels)")
        ax.set_title(f"removed and kept, {name}: {int(sel.sum())} kept shown in grey")
        ax.legend(fontsize=8, loc="upper right", markerscale=2)

    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
