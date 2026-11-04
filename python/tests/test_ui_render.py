"""The viewer's drawing without a window (mxi_ui/render.py): a frame reduced
by block maxima at a whole scale, zoom by whole levels about the cursor, the
colour maps, the PPM image, the resolution at a pixel."""

import math

import numpy as np
import pytest

from mxi_ui import render


def test_a_one_pixel_spot_survives_any_reduction():
    frame = np.zeros((4362, 4148), np.float32)
    frame[1234, 2345] = 50
    for scale in (2.0, 7.0, 22.0):
        image, at = render.reduce(frame, (0, 0, 4148, 4362), scale)
        assert image.max() == 50 and at == (0, 0)
        assert image.shape == (math.ceil(4362 / scale), math.ceil(4148 / scale))
    zoomed, at = render.reduce(frame, (2340, 1230, 2350, 1240), 0.25)
    assert zoomed.shape == (40, 40) and at == (2340, 1230) and zoomed.max() == 50


def test_a_region_off_the_detector_is_clipped_to_it():
    frame = np.ones((100, 200), np.float32)
    image, at = render.reduce(frame, (-50, -20, 150, 80), 2.0)
    assert at == (0, 0) and image.shape == (40, 75)


@pytest.mark.parametrize(
    "given, whole", [(7.3, 8.0), (1.0, 1.0), (0.6, 1.0), (0.3, 1 / 3)]
)
def test_scales_are_whole(given, whole):
    # Towards the coarser level either way: a little more shown, not less.
    assert math.isclose(render.snap(given), whole)


def test_zoom_moves_by_whole_levels_and_keeps_the_point_under_the_cursor():
    v = render.View((4148, 4362), (800, 700))
    levels = [v.scale]
    for _ in range(10):
        x, y = v.to_detector(300, 200)
        v.zoom(1.25, 300, 200)
        levels.append(v.scale)
        assert all(
            math.isclose(a, b, abs_tol=1e-6)
            for a, b in zip(v.to_canvas(x, y), (300, 200))
        )
    assert levels[:8] == [7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.5]
    assert all(b < a for a, b in zip(levels, levels[1:]))
    v.pan(10, -5)
    assert all(
        math.isclose(a, b)
        for a, b in zip(v.to_canvas(*v.to_detector(40, 50)), (40, 50))
    )


def test_the_colour_maps_at_their_ends():
    values = np.array([[0.0, 10.0, 20.0]])
    grey = render.colour(values, 10)
    assert (
        grey[0, 0].tolist() == [255] * 3 and grey[0, 1].tolist() == [0] * 3
    ), "dark spots"
    assert render.colour(values, 10, "inverted")[0, 1].tolist() == [255] * 3
    heat = render.colour(values, 10, "heat")
    assert heat[0, 0].tolist() == [0, 0, 0] and heat[0, 2].tolist() == [255, 255, 255]
    assert (
        render.ppm(grey).startswith(b"P6\n3 1\n255\n")
        and len(render.ppm(grey)) == 11 + 9
    )


def test_the_resolution_at_a_pixel():
    doc = {
        "experiment": [{"beam": 0, "detector": 0}],
        "beam": [{"direction": [0, 0, 1], "wavelength": 1.0}],
        "detector": [
            {
                "panels": [
                    {
                        "fast_axis": [1, 0, 0],
                        "slow_axis": [0, -1, 0],
                        "origin": [-100, 100, -100],
                        "pixel_size": [0.1, 0.1],
                    }
                ]
            }
        ],
    }
    two_theta = math.atan(100 / 100)
    assert math.isclose(
        render.resolution(doc, 0, 2000, 1000), 1 / (2 * math.sin(two_theta / 2))
    )
    assert render.resolution(doc, 0, 1000, 1000) is None, "on the beam"
    assert render.resolution({}, 0, 1, 1) is None, "no geometry"


def test_the_fast_reduction_is_the_block_maxima():
    # Against the plain reduction it replaced, edges and all.
    rng = np.random.default_rng(7)
    frame = rng.poisson(3, (517, 389)).astype(np.float32)
    for block in (2, 3, 8):
        h, w = frame.shape
        padded = np.pad(frame, ((0, -h % block), (0, -w % block)))
        plain = padded.reshape(padded.shape[0] // block, block, -1, block).max(
            axis=(1, 3)
        )
        fast, _ = render.reduce(frame, (0, 0, w, h), float(block))
        assert np.array_equal(fast, plain), block
