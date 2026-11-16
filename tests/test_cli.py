"""CLI smoke tests: every subcommand runs end to end on real files."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "examples"
FASTA_16S = EXAMPLES / "ecoli_16s_rrnB_J01695.fasta"
SEQ_28 = "CGTTCCAAAGATGTGGGCATGAGCTTAC"
PRIMER_27F = "AGAGTTTGATCATGGCTCAG"
PRIMER_1492R = "GGTTACCTTGTTACGACTT"


def run_cli(*args: str, cwd: Path = ROOT) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "primer_designer", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=120,
    )


def test_cli_tm():
    result = run_cli("tm", SEQ_28)
    assert result.returncode == 0
    assert "Tm Wallace  84.00" in result.stdout
    assert "Tm NN       63.32" in result.stdout
    assert "GC content  50.00" in result.stdout


def test_cli_tm_rejects_bad_sequence():
    result = run_cli("tm", "ACGTN")
    assert result.returncode == 2
    assert "error:" in result.stderr


def test_cli_check_universal_pair():
    result = run_cli(
        "check",
        "--forward", PRIMER_27F,
        "--reverse", PRIMER_1492R,
        "--template", str(FASTA_16S),
    )
    assert result.returncode == 0
    assert "cross-dimer" in result.stdout
    assert "delta Tm" in result.stdout
    assert "evaluation:" in result.stdout


def test_cli_design_full_16s_writes_outputs(tmp_path):
    out_tsv = tmp_path / "full.tsv"
    out_json = tmp_path / "full.json"
    result = run_cli(
        "design",
        "--fasta", str(FASTA_16S),
        "--start", "27", "--end", "1492",
        "--name", "16s_full",
        "--out", str(out_tsv),
        "--json", str(out_json),
    )
    assert result.returncode == 0
    assert "template" in result.stdout
    assert "product" in result.stdout
    with open(out_tsv, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 5
    assert all(1470 <= int(r["product_length"]) <= 1530 for r in rows)
    data = json.loads(out_json.read_text(encoding="utf-8"))
    assert len(data) == 5
    assert data[0]["forward"]["sequence"] == rows[0]["forward_seq"]


def test_cli_design_v4(tmp_path):
    result = run_cli(
        "design",
        "--fasta", str(FASTA_16S),
        "--start", "536", "--end", "787",
        "--name", "16s_v4",
        "--top", "3",
    )
    assert result.returncode == 0
    assert result.stdout.count("rank ") == 3


def test_cli_batch(tmp_path):
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "name\tstart\tend\n16s_full\t27\t1492\n16s_v4\t536\t787\n",
        encoding="utf-8",
    )
    out_tsv = tmp_path / "batch.tsv"
    out_json = tmp_path / "batch.json"
    result = run_cli(
        "batch",
        "--fasta", str(FASTA_16S),
        "--targets", str(targets),
        "--top", "2",
        "--out", str(out_tsv),
        "--json", str(out_json),
    )
    assert result.returncode == 0
    assert "16s_full:" in result.stdout
    assert "16s_v4:" in result.stdout
    with open(out_tsv, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 4  # 2 targets x top 2
    data = json.loads(out_json.read_text(encoding="utf-8"))
    assert set(data) == {"16s_full", "16s_v4"}


def test_cli_design_missing_file_fails():
    result = run_cli(
        "design",
        "--fasta", "does-not-exist.fasta",
        "--start", "27", "--end", "1492",
    )
    assert result.returncode == 2
    assert "error:" in result.stderr


def test_cli_design_impossible_constraints_fails():
    result = run_cli(
        "design",
        "--fasta", str(FASTA_16S),
        "--start", "27", "--end", "1492",
        "--min-product", "5000", "--max-product", "6000",
    )
    assert result.returncode == 2
    assert "no primer pair" in result.stderr
