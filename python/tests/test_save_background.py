"""mxi_integrate --save-background-parameters: the background's dispersion
written only when asked, every other column the same either way.

While the dispersion is under investigation (docs/backstop.md) the default
table is what it was. Needs MXI_INTEGRATE, MXI_POSTREFINE_EXPT and
MXI_POSTREFINE_REFL.
"""

import os
import subprocess

import numpy as np
import pytest

from mxeq import refl

BINARY = os.environ.get("MXI_INTEGRATE")
EXPT = os.environ.get("MXI_POSTREFINE_EXPT")
REFL = os.environ.get("MXI_POSTREFINE_REFL")
needs = pytest.mark.skipif(
    not (BINARY and EXPT and REFL),
    reason="set MXI_INTEGRATE, MXI_POSTREFINE_EXPT, MXI_POSTREFINE_REFL",
)
COLUMNS = {
    "background.dispersion",
    "background.dispersion_trimmed",
    "num_pixels.background_trimmed",
}


def integrate(tmp_path, name, *extra):
    run = subprocess.run(
        [
            BINARY,
            EXPT,
            REFL,
            "--last-image",
            "30",
            "--threads",
            "2",
            "-o",
            name + ".refl",
            "--output-expt",
            name + ".expt",
            *extra,
        ],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    return refl.load(str(tmp_path / (name + ".refl"))).columns


@needs
def test_the_dispersion_is_written_only_when_asked(tmp_path):
    plain = integrate(tmp_path, "plain")
    saved = integrate(tmp_path, "saved", "--save-background-parameters")
    assert not COLUMNS & set(plain)
    assert set(saved) - set(plain) == COLUMNS
    for name in plain:
        assert np.array_equal(
            np.asarray(plain[name]), np.asarray(saved[name]), equal_nan=True
        ), name
    d = np.asarray(saved["background.dispersion_trimmed"])
    assert 0.9 < np.nanmedian(d) < 1.1
