"""Prevent broken notebook structure and a return to hidden exercise solutions."""

from pathlib import Path
import re
import sys
import nbformat
import pytest
from IPython.core.inputtransformer2 import TransformerManager

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "maintenance"))
from catalog_data import catalog


def test_catalog_matches_published_notebooks():
    assert set(catalog()) == {p.stem for p in (ROOT / "python").glob("*.ipynb")}
    assert len([t for t in catalog().values() if "canonical" not in t]) == 58


@pytest.mark.parametrize("slug", sorted(catalog()))
def test_notebook_contract(slug):
    nb = nbformat.read(ROOT / "python" / (slug + ".ipynb"), as_version=4)
    nbformat.validate(nb)
    assert nb.cells[0].cell_type == "markdown"
    assert nb.cells[0].source.startswith("# ")
    refs = next(c.source for c in nb.cells if c.id == "nt-references")
    assert 5 <= len(re.findall(r"(?m)^[*-] ", refs)) <= 10
    assert nb.metadata.kernelspec.name == "python3"
    for i, cell in enumerate(nb.cells):
        if cell.cell_type != "code":
            continue
        assert not re.search(
            r"nt_solutions|put your code|insert your code|not coded", cell.source, re.I
        )
        compile(TransformerManager().transform_cell(cell.source), f"{slug}:{i}", "exec")
