"""Tests for FASTA / targets parsing and TSV / JSON output."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from primer_designer.design import Target, design_pair
from primer_designer.io import (
    PAIR_COLUMNS,
    parse_fasta_text,
    read_fasta,
    read_targets,
    write_batch_json,
    write_batch_tsv,
    write_pairs_json,
    write_pairs_tsv,
)


def test_parse_fasta_multi_record_wrapped_lowercase():
    text = """>seq1 first record
acgtacgtac
gtacgtacgt
>seq2
TTTTCCCCAAAAGGGG
"""
    records = parse_fasta_text(text)
    assert records == [
        ("seq1 first record", "ACGTACGTACGTACGTACGT"),
        ("seq2", "TTTTCCCCAAAAGGGG"),
    ]


def test_parse_fasta_errors():
    with pytest.raises(ValueError):
        parse_fasta_text("ACGTACGT\n")
    with pytest.raises(ValueError):
        parse_fasta_text("")
    with pytest.raises(ValueError):
        parse_fasta_text(">only-header\n")


def test_read_fasta_16s_example():
    path = Path(__file__).resolve().parent.parent / "examples" / "ecoli_16s_rrnB_J01695.fasta"
    records = read_fasta(path)
    assert len(records) == 1
    assert records[0][0].startswith("ecoli_16s_rrnB_J01695")
    assert len(records[0][1]) == 1542
    assert set(records[0][1]) <= set("ACGT")


def test_read_targets_with_header_and_comments(tmp_path):
    tsv = tmp_path / "targets.tsv"
    tsv.write_text(
        "# example targets\n"
        "name\tstart\tend\n"
        "16s_full\t27\t1492\n"
        "16s_v4\t536\t787\n",
        encoding="utf-8",
    )
    targets = read_targets(tsv)
    assert targets == [
        Target("16s_full", 27, 1492),
        Target("16s_v4", 536, 787),
    ]


def test_read_targets_without_header(tmp_path):
    tsv = tmp_path / "targets.tsv"
    tsv.write_text("mygene\t100\t300\n", encoding="utf-8")
    assert read_targets(tsv) == [Target("mygene", 100, 300)]


def test_read_targets_rejects_bad_coordinates(tmp_path):
    tsv = tmp_path / "targets.tsv"
    tsv.write_text("name\tstart\tend\nbad\t300\t100\n", encoding="utf-8")
    with pytest.raises(ValueError):
        read_targets(tsv)
    tsv.write_text("name\tstart\tend\nbad\tx\t100\n", encoding="utf-8")
    with pytest.raises(ValueError):
        read_targets(tsv)


def test_pairs_tsv_round_trip(tmp_path, template_16s):
    pairs = design_pair(template_16s, 535, 787, "16s_v4")
    out = tmp_path / "pairs.tsv"
    assert write_pairs_tsv(pairs, out) == out
    with open(out, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == len(pairs)
    assert list(rows[0].keys()) == PAIR_COLUMNS
    assert rows[0]["target"] == "16s_v4"
    assert rows[0]["rank"] == "1"
    assert rows[0]["forward_seq"] == pairs[0].forward.seq
    assert int(rows[0]["product_length"]) == pairs[0].product_length
    assert rows[0]["product_seq"] == pairs[0].product_seq


def test_pairs_json_round_trip(tmp_path, template_16s):
    pairs = design_pair(template_16s, 535, 787, "16s_v4")
    out = tmp_path / "pairs.json"
    write_pairs_json(pairs, out)
    data = json.loads(out.read_text(encoding="utf-8"))
    assert len(data) == len(pairs)
    assert data[0]["forward"]["sequence"] == pairs[0].forward.seq
    assert data[0]["product_length"] == pairs[0].product_length


def test_batch_outputs_include_error_row(tmp_path, template_16s):
    results = {
        "ok": design_pair(template_16s, 535, 787, "ok"),
        "bad": "no primer pair satisfies the constraints",
    }
    tsv = write_batch_tsv(results, tmp_path / "batch.tsv")
    with open(tsv, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    error_rows = [r for r in rows if r["forward_seq"].startswith("ERROR:")]
    assert len(error_rows) == 1
    assert error_rows[0]["target"] == "bad"
    data = json.loads(
        write_batch_json(results, tmp_path / "batch.json").read_text(
            encoding="utf-8"
        )
    )
    assert data["bad"] == {"error": "no primer pair satisfies the constraints"}
    assert len(data["ok"]) == len(results["ok"])
