#!/usr/bin/env python3
"""mxi_view: the images of an experiment list, with what mxi found on them, in
napari -- and the detector's untrusted regions drawn and saved.

    mxi_view imported.rflx            # the images, a sweep at a time
    mxi_view strong.rflx              # with the strong spots
    mxi_view integrated.rflx          # with the spots and the predictions

The images are a stack, a slider through them, each shown as the sum of the
number set in the panel (one to begin with) from it, so that a shadow shows
against the background. Spots (xyzobs) and predictions (xyzcal) are points,
each on its own images. Untrusted regions -- a backstop, its arm, a shadow --
are a shapes layer to draw on: rectangles, ellipses (a circle where round)
and polygons, napari's own editing and undo. Save writes the experiment list
with every sweep's regions to masked.rflx (-o to choose), the input left as it
was, through mxi_untrusted, and keeps them in mxi_mask.json:
`mxi_untrusted other.rflx --set mxi_mask.json` applies them again with no
window. Spot finding and integration leave the regions out.

Needs napari (pip install "napari[all]"), with numpy, h5py and hdf5plugin.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import mxi_tools  # noqa: E402

VERSION = "@PROJECT_VERSION@"


class View:
    """The window: one sweep at a time, chosen in the panel."""

    def __init__(self, path: str, output: str, window: int, sweep: int = 0):
        import napari

        self.napari = napari
        self.path = os.path.abspath(path)
        self.output = output
        self.window = window
        self.document = mxi_tools.read_experiments(self.path)
        self.sweeps = mxi_tools.sweeps(self.document)
        if not self.sweeps:
            raise SystemExit(f"{path} has no image sets")
        self.table = mxi_tools.read_reflections(
            self.path,
            columns={
                "id",
                "xyzobs.px.value",
                "xyzobs.px.variance",
                "xyzcal.px",
                "bbox",
            },
        )
        self.regions = {s["index"]: list(s["shapes"]) for s in self.sweeps}
        self.images = {}
        self.current = None
        self.previous = None
        # Undo: each sweep's regions as they were before each change to them,
        # latest last; napari's shapes layer has none of its own.
        self.history = {}
        self.state = []
        self.restoring = False
        self.viewer = napari.Viewer(title=f"mxi_view {os.path.basename(path)}")
        self.markers = {}
        self.viewer.dims.events.current_step.connect(lambda _e: self.update_markers())
        self.panel()
        self.show(min(max(sweep, 0), len(self.sweeps) - 1))

    # --- the panel
    def panel(self):
        from magicgui.widgets import ComboBox, Container, Label, PushButton, SpinBox

        names = [f"{s['index']}: {os.path.basename(s['master'])}" for s in self.sweeps]
        self.choice = ComboBox(label="sweep", choices=names, value=names[0])
        self.choice.changed.connect(lambda v: self.show(names.index(v)))
        self.summed = SpinBox(label="images summed", value=self.window, min=1, max=1000)
        self.summed.changed.connect(self.set_window)
        same = PushButton(text="Same as last")
        same.changed.connect(self.same_as_last)
        undo = PushButton(text="Undo")
        undo.changed.connect(self.undo)
        clear = PushButton(text="Clear")
        clear.changed.connect(self.clear)
        save = PushButton(text=f"Save to {self.output}")
        save.changed.connect(self.save)
        self.status = Label(value="")
        self.container = Container(
            widgets=[self.choice, self.summed, undo, clear, same, save, self.status]
        )
        self.viewer.window.add_dock_widget(self.container, name="mxi", area="right")

    # --- a sweep's layers
    def keep_regions(self):
        """The shapes layer's regions, as untrusted regions, for the sweep on
        show."""
        if self.current is None or "untrusted" not in self.viewer.layers:
            return
        self.regions[self.sweeps[self.current]["index"]] = self.current_regions()

    def current_regions(self) -> list:
        """The shapes layer's regions, as untrusted regions."""
        layer = self.viewer.layers["untrusted"]
        regions = []
        for kind, vertices in zip(layer.shape_type, layer.data):
            r = mxi_tools.shape_to_untrusted(kind, vertices)
            if r is not None:
                regions.append(r)
        return regions

    def show(self, k: int, keep: bool = True):
        # The layer's regions kept before another sweep's are drawn -- but not
        # after Same as last has just replaced them, which keeping the layer
        # still on show would undo.
        if keep:
            self.keep_regions()
        if self.current is not None and self.current != k:
            self.previous = self.current
        self.current = k
        s = self.sweeps[k]
        for name in ("untrusted", "predictions", "spots", "images"):
            if name in self.viewer.layers:
                self.viewer.layers.remove(name)
        if s["master"] not in self.images:
            self.images[s["master"]] = mxi_tools.Images(s["master"], self.window)
        images = self.images[s["master"]]
        images.window = self.window
        self.viewer.add_image(
            images,
            name="images",
            colormap="gray_r",
            contrast_limits=(
                0,
                mxi_tools.display_top(images.summed(max(0, s["first"] - 1))),
            ),
        )
        self.add_markers()
        kinds, data = [], []
        for region in self.regions[s["index"]]:
            kind, vertices = mxi_tools.untrusted_to_shape(region)
            kinds.append(kind)
            data.append(vertices)
        self.viewer.add_shapes(
            data,
            shape_type=kinds,
            name="untrusted",
            ndim=2,
            face_color=[1, 0, 0, 0.25],
            edge_color="red",
            edge_width=2,
        )
        layer = self.viewer.layers["untrusted"]
        layer.events.data.connect(self.on_change)
        layer.bind_key("Control-Z", lambda _layer: self.undo(), overwrite=True)
        self.state = self.current_regions()
        self.viewer.dims.set_point(0, max(0, s["first"] - 1))
        self.say(
            f"images {s['first']}-{s['last']}; {self.counts[0]} spots, "
            f"{self.counts[1]} predicted; {len(self.regions[s['index']])} regions"
        )

    def add_markers(self):
        """The spots and predictions of the sweep on show, as markers for the
        image shown only (mxi_tools.Markers), picked out afresh as the slider
        or the window moves: a few hundred, not the whole sweep's millions,
        which had taken some 15 s to start on an integrated 16M sweep. Red on a
        spot's centroid's image, a paler red on the others it is recorded on;
        the predictions cyan, hidden unless already shown, and not picked out
        while hidden."""
        s = self.sweeps[self.current]
        visible = (
            self.viewer.layers["predictions"].visible
            if "predictions" in self.viewer.layers
            else False
        )
        for name in ("predictions", "spots"):
            if name in self.viewer.layers:
                self.viewer.layers.remove(name)
        members = s["experiments"]
        self.markers = {
            "spots": mxi_tools.Markers(self.table, "xyzobs.px.value", members),
            "predictions": mxi_tools.Markers(self.table, "xyzcal.px", members),
        }
        self.counts = [len(self.markers["spots"]), len(self.markers["predictions"])]
        empty = mxi_tools.np.zeros((0, 2))
        # At least ten screen pixels, however far out: a marker 12 detector
        # pixels across is a speck on a 16M detector seen whole (23 October).
        limits = (10, 10000)
        if self.counts[0]:
            self.viewer.add_points(
                empty,
                name="spots",
                size=12,
                face_color="transparent",
                border_color="red",
                border_width=0.15,
                ndim=2,
                canvas_size_limits=limits,
            )
        if self.counts[1]:
            layer = self.viewer.add_points(
                empty,
                name="predictions",
                size=8,
                symbol="cross",
                face_color="cyan",
                ndim=2,
                visible=visible,
                canvas_size_limits=limits,
            )
            layer.events.visible.connect(lambda _e: self.update_markers())
        self.update_markers()
        if "untrusted" in self.viewer.layers:  # kept on top, to draw on
            layers = self.viewer.layers
            layers.move(layers.index("untrusted"), len(layers))

    def update_markers(self):
        """The markers on the image shown -- and the window's."""
        frame = int(self.viewer.dims.current_step[0])
        for name, markers in getattr(self, "markers", {}).items():
            if name not in self.viewer.layers or not len(markers):
                continue
            layer = self.viewer.layers[name]
            if not layer.visible:
                continue
            yx, own = markers.on(frame, self.window)
            layer.data = yx
            if name == "spots" and len(yx):
                layer.border_color = mxi_tools.np.where(
                    own[:, None], [1.0, 0.0, 0.0, 1.0], [1.0, 0.55, 0.55, 0.8]
                )

    def set_window(self, n: int):
        self.window = max(1, int(n))
        for images in self.images.values():
            images.window = self.window
            images._cache.clear()
        if "images" in self.viewer.layers:
            self.viewer.layers["images"].refresh()
            self.update_markers()

    # --- undo
    def on_change(self, event):
        """A change finished -- added, removed, or a shape changed -- the
        regions as they were before it kept to undo; the events while a shape
        is drawn or dragged, and Undo's own, not."""
        if self.restoring:
            return
        action = getattr(event, "action", None)
        action = str(getattr(action, "value", action))
        if action not in ("added", "removed", "changed"):
            return
        now = self.current_regions()
        if now == self.state:
            return
        self.history.setdefault(self.sweeps[self.current]["index"], []).append(
            self.state
        )
        self.state = now

    def set_layer(self, regions: list):
        """The shapes layer made these regions, without recording it."""
        layer = self.viewer.layers["untrusted"]
        self.restoring = True
        try:
            layer.selected_data = set(range(len(layer.data)))
            layer.remove_selected()
            for r in regions:
                kind, vertices = mxi_tools.untrusted_to_shape(r)
                layer.add([vertices], shape_type=[kind])
        finally:
            self.restoring = False
        self.state = [dict(r) for r in regions]
        self.regions[self.sweeps[self.current]["index"]] = [dict(r) for r in regions]

    def undo(self):
        past = self.history.get(self.sweeps[self.current]["index"], [])
        if not past:
            self.say("nothing to undo on this sweep")
            return
        self.set_layer(past.pop())
        self.say(f"undone: {len(self.state)} regions; {len(past)} more to undo")

    def clear(self):
        if not self.current_regions():
            return
        self.history.setdefault(self.sweeps[self.current]["index"], []).append(
            self.current_regions()
        )
        self.set_layer([])
        self.say("cleared; Undo brings them back")

    def same_as_last(self):
        if self.previous is None:
            self.say("no other sweep shown yet")
            return
        self.keep_regions()
        source = self.regions[self.sweeps[self.previous]["index"]]
        self.history.setdefault(self.sweeps[self.current]["index"], []).append(
            self.current_regions()
        )
        self.set_layer(source)

    def save(self):
        self.keep_regions()
        try:
            report = mxi_tools.save_untrusted(self.path, self.regions, self.output)
        except RuntimeError as e:
            self.say(f"not saved: {e}")
            return
        self.say(report.strip().splitlines()[-1] + "; and mxi_mask.json")
        print(report, end="")

    def say(self, text: str):
        self.status.value = text


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="mxi_view",
        description="The images of an experiment list, with the spots and "
        "predictions mxi found on them, in napari; and the detector's untrusted "
        "regions drawn and saved to masked.rflx through mxi_untrusted.",
    )
    parser.add_argument(
        "experiments", help="a .rflx or .expt: imported, strong, integrated ..."
    )
    parser.add_argument(
        "-o", default="masked.rflx", help="where Save writes (masked.rflx)"
    )
    parser.add_argument(
        "--sum", type=int, default=1, help="images summed for each shown (1)"
    )
    parser.add_argument(
        "--sweep", type=int, default=0, help="the sweep shown first (0)"
    )
    parser.add_argument("--version", action="version", version=f"mxi_view {VERSION}")
    args = parser.parse_args(argv)
    try:
        import napari  # noqa: F401
    except ImportError:
        print(
            'mxi_view needs napari: pip install "napari[all]" (with numpy, h5py, '
            "hdf5plugin). mxi_untrusted sets regions with no window.",
            file=sys.stderr,
        )
        return 2
    view = View(args.experiments, args.o, args.sum, args.sweep)
    view.napari.run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
