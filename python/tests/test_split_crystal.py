"""mxi index --max-lattices on a split crystal (26 November): o1328 of
A71EV2AZ, two domains 1.8 degrees apart. Once the second domain was found,
most spots were claimed by both lattices, --ambiguous exclude (the default)
left 0.7 per cent of them, a third lattice of leftovers (a 5.8 A axis) was
taken, and the refinement was wrecked: cells of 64 A, the distance 175.6 mm.
Now the leftover lattice is refused and the exclusion declined.

Run with MXI_SPLIT_CRYSTAL_DATA naming o1328's strong.rflx.
"""

import os
import re
import subprocess

import pytest

BUILD = os.environ.get("MXI_BUILD")
MXI = os.path.join(BUILD, "mxi") if BUILD else None
DATA = os.environ.get("MXI_SPLIT_CRYSTAL_DATA")

pytestmark = pytest.mark.skipif(
    not (MXI and DATA), reason="set MXI_BUILD and MXI_SPLIT_CRYSTAL_DATA"
)


def test_a_split_crystal_indexed_with_lattices_to_spare_is_not_wrecked(
    tmp_path,
):
    run = subprocess.run(
        [MXI, "index", DATA, "--max-lattices", "4", "-o", "idx.rflx"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr[-2000:]
    out = run.stdout
    assert "unlike the first's" in out, "the leftover lattice refused"
    fraction = float(re.search(r"spots indexed \(([\d.]+)%\)", out).group(1))
    assert fraction > 50, f"most spots indexed, not {fraction} per cent"
    cells = re.findall(r"lattice \d+: \d+ spots; ([\d.]+) ([\d.]+) ([\d.]+) A", out)
    assert len(cells) == 2, out
    for cell in cells:
        for edge in map(float, cell):
            assert 55 < edge < 85, f"a sensible cell, not {cell}"
    distance = float(re.search(r"refined together .* and ([\d.]+) mm", out).group(1))
    assert abs(distance - 169.0) < 2.0, f"the distance near 169 mm, not {distance}"
