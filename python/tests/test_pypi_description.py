"""PyPI's page is the package's readme, python/PYPI.md -- about mxi, the package
it is now, not mxeq (Graeme, 26 November). PyPI resolves no relative link, so
every link is absolute, and those into the repository name files that exist."""

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]
TEXT = (ROOT / "python" / "PYPI.md").read_text()
LINK = re.compile(r"\]\(([^)]+)\)")
BLOB = "https://github.com/graeme-winter/mxi/blob/main/"


def test_the_package_readme_is_the_pypi_description():
    pyproject = (ROOT / "python" / "pyproject.toml").read_text()
    assert re.search(r'^readme = "PYPI.md"$', pyproject, re.M)
    assert "[project.urls]" in pyproject


def test_every_link_is_absolute_and_every_file_linked_exists():
    links = LINK.findall(TEXT)
    assert links
    for link in links:
        assert link.startswith("https://"), link
        if link.startswith(BLOB):
            path = link[len(BLOB) :].split("#")[0]
            assert (ROOT / path).exists(), path


def test_it_says_how_to_install_both_and_not_both():
    assert "pip install mxi " in TEXT and "pip install mxi-cuda" in TEXT
    assert "not both" in TEXT
