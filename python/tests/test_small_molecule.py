"""A small molecule's crystal in a centrosymmetric group, symmetry to scaling.

mxi symmetry --small-molecule considers every space group with the Patterson
group found, centrosymmetric among them; mxi scale, in a centrosymmetric group,
has no anomalous statistics to give -- every reflection is centric, its Friedel
mate itself -- and says so, where it had printed them anyway (completeness 0,
multiplicity the multiplicity), and leaves the anomalous columns out of its
summary. Graeme's R -3 c crystal, 13 November.

Needs an integrated small molecule whose group is centrosymmetric:
MXI_SMALL_MOLECULE, the integrated.rflx.
"""

import os
import subprocess

import pytest

BUILD = os.path.dirname(os.environ.get("MXI_FIND", "")) or os.environ.get(
    "MXI_BUILD", ""
)
MXI = os.path.join(BUILD, "mxi")
DATA = os.environ.get("MXI_SMALL_MOLECULE")

pytestmark = pytest.mark.skipif(
    not (os.path.exists(MXI) and DATA),
    reason="set MXI_SMALL_MOLECULE, an integrated.rflx in a centrosymmetric group",
)


def run(*words, cwd):
    r = subprocess.run([MXI, *words], capture_output=True, text=True, cwd=cwd)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    return r.stdout


def test_a_centrosymmetric_group_chosen_and_scaled_without_anomalous_figures(
    tmp_path,
):
    from mxeq import rflx

    out = run("symmetry", DATA, "--small-molecule", "-o", "sym.rflx", cwd=tmp_path)
    assert "Mean |E^2-1|" in out
    assert "centrosymmetric" in out
    out = run("scale", "sym.rflx", "-o", "scaled.rflx", cwd=tmp_path)
    assert "none: the group is centrosymmetric" in out
    assert "Anomalous completeness" not in out
    assert "Anomalous multiplicity" not in out
    step = [
        s
        for s in rflx.summary_dict(str(tmp_path / "scaled.rflx"))["steps"]
        if s.get("program") == "mxi_scale"
    ][-1]
    assert step["centrosymmetric"] is True
    assert not [c for c in step["merging_summary"] if "anom" in c]
