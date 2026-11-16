"""Shared pytest fixtures for primer-designer-lite tests."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Make the project root importable when pytest is run from elsewhere.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from primer_designer.io import read_fasta  # noqa: E402


@pytest.fixture(scope="session")
def template_16s() -> str:
    """E. coli 16S rRNA gene (J01695:1518..3059), 1542 bp."""
    path = Path(__file__).resolve().parent.parent / "examples" / "ecoli_16s_rrnB_J01695.fasta"
    records = read_fasta(path)
    assert len(records) == 1
    assert len(records[0][1]) == 1542
    return records[0][1]


@pytest.fixture(scope="session")
def template_synthetic() -> str:
    """Synthetic 600 bp template (random.seed(20261116))."""
    path = Path(__file__).resolve().parent.parent / "examples" / "synthetic_example.fasta"
    records = read_fasta(path)
    assert len(records) == 1
    assert len(records[0][1]) == 600
    return records[0][1]
