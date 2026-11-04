"""mxeq, the referee, imports nothing of mxi's own Python, though it shares the
wheel, mxi-python, with it: it judges mxi's files, so it must not share mxi's
reading of them -- a fault in mxi_tools would otherwise pass for agreement.
Every module mxeq has, against every mxi_ module beside it in src/: today
mxi_tools and mxi_view, and whatever is added, the interface among them."""

import ast
import pathlib

import pytest

SRC = pathlib.Path(__file__).resolve().parents[1] / "src"
MXEQ = sorted(str(p.relative_to(SRC)) for p in (SRC / "mxeq").rglob("*.py"))


def test_there_is_something_to_check():
    assert len(MXEQ) > 10
    assert (SRC / "mxi_tools.py").exists(), "mxi's own Python, beside mxeq"


@pytest.mark.parametrize("module", MXEQ)
def test_mxeq_imports_nothing_of_mxis(module):
    tree = ast.parse((SRC / module).read_text())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            imported.add(node.module.split(".")[0])
    mxis = sorted(name for name in imported if name.startswith("mxi_") or name == "mxi")
    assert not mxis, f"{module} imports {mxis}"
