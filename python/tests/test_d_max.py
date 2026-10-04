"""mxi_scale --d-max: no observation with d above the limit is scaled.

For leaving out the reflections beside a backstop, whose observations in its
shadow integrate to about nothing: the low resolution limit, as --d-min is the
high. Needs MXI_SCALE, MXI_SCALE_EXPT and MXI_SCALE_REFL.
"""

import os
import subprocess

import numpy as np
import pytest

from mxeq import equivalents as eq
from mxeq import refl

BINARY = os.environ.get("MXI_SCALE")
EXPT = os.environ.get("MXI_SCALE_EXPT")
REFL = os.environ.get("MXI_SCALE_REFL")
needs = pytest.mark.skipif(
    not (BINARY and EXPT and REFL),
    reason="set MXI_SCALE, MXI_SCALE_EXPT, MXI_SCALE_REFL",
)


def scaled_d(tmp_path, name, *extra):
    run = subprocess.run(
        [
            BINARY,
            EXPT,
            REFL,
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
    t = refl.load(str(tmp_path / (name + ".refl")))
    flags = np.asarray(t.columns["flags"]).astype(np.int64)
    return np.asarray(t.columns["d"])[(flags & eq.SCALED) != 0]


@needs
def test_nothing_beyond_the_low_resolution_limit_is_scaled(tmp_path):
    every = scaled_d(tmp_path, "every")
    limit = float(np.sort(every)[-20])  # leave out the twenty lowest
    cut = scaled_d(tmp_path, "cut", "--d-max", f"{limit:.4f}")
    assert every.max() > limit
    assert cut.max() <= limit
    assert len(cut) >= len(every) - 40  # only those, near enough
