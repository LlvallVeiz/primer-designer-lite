"""File I/O: FASTA and target tables in, TSV/JSON reports out.

Formats
-------
* FASTA: standard multi-record files; headers start with ``>``,
  sequences may wrap over lines and may be lower-case. Records are
  returned as ``(header, sequence)`` tuples with the header text
  (without ``>``) and the sequence upper-cased.
* Targets TSV: a header row plus one target per line with columns
  ``name``, ``start``, ``end`` -- coordinates are **1-based inclusive**
  (the convention used by primer names such as 27F/1492R and by
  GenBank feature tables), converted to 0-based half-open coordinates
  by :class:`primer_designer.design.Target`.
* Pairs TSV: one row per designed pair with a ``rank`` column (1 =
  best) and one column per reported value; the two primer sequences
  and the product sequence are included so the output is self-contained.
* Pairs JSON: the same data as a list of pair objects (see
  :meth:`primer_designer.design.PrimerPair.to_dict`).
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from .design import PrimerPair, Target

__all__ = [
    "read_fasta",
    "parse_fasta_text",
    "read_targets",
    "pairs_to_rows",
    "write_pairs_tsv",
    "write_pairs_json",
    "write_batch_tsv",
    "write_batch_json",
]

PAIR_COLUMNS = [
    "target",
    "rank",
    "forward_seq",
    "forward_start_0based",
    "forward_tm_nn_c",
    "forward_tm_wallace_c",
    "forward_gc_percent",
    "forward_score",
    "reverse_seq",
    "reverse_start_0based",
    "reverse_tm_nn_c",
    "reverse_tm_wallace_c",
    "reverse_gc_percent",
    "reverse_score",
    "delta_tm_c",
    "anneal_temp_c",
    "product_length",
    "product_seq",
    "pair_score",
    "penalties",
]


def parse_fasta_text(text: str) -> list[tuple[str, str]]:
    """Parse FASTA text into ``[(header, sequence), ...]``."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                sequence = "".join(chunks).upper()
                if not sequence:
                    raise ValueError(
                        f"FASTA record '{header}' has an empty sequence"
                    )
                records.append((header, sequence))
            header = line[1:].strip()
            chunks = []
        else:
            if header is None:
                raise ValueError("FASTA sequence found before any header")
            chunks.append(line)
    if header is not None:
        sequence = "".join(chunks).upper()
        if not sequence:
            raise ValueError(f"FASTA record '{header}' has an empty sequence")
        records.append((header, sequence))
    if not records:
        raise ValueError("no FASTA records found")
    return records


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Read a FASTA file into ``[(header, sequence), ...]``."""
    return parse_fasta_text(Path(path).read_text(encoding="utf-8"))


def read_targets(path: str | Path) -> list[Target]:
    """Read a targets TSV (columns ``name``, ``start``, ``end``).

    Coordinates are 1-based inclusive. A ``#`` comment line or a
    header row is accepted; extra columns are ignored.
    """
    targets: list[Target] = []
    with open(path, encoding="utf-8", newline="") as handle:
        rows = [
            row
            for row in csv.reader(handle, delimiter="\t")
            if row and not row[0].lstrip().startswith("#")
        ]
    if not rows:
        raise ValueError("targets file is empty")
    header = [cell.strip().lower() for cell in rows[0]]
    if "name" in header and "start" in header and "end" in header:
        rows = rows[1:]
        name_i, start_i, end_i = (
            header.index("name"),
            header.index("start"),
            header.index("end"),
        )
    else:
        name_i, start_i, end_i = 0, 1, 2
    for lineno, row in enumerate(rows, start=2):
        if len(row) <= max(name_i, start_i, end_i):
            raise ValueError(f"targets line {lineno}: not enough columns")
        try:
            start = int(row[start_i])
            end = int(row[end_i])
        except ValueError as exc:
            raise ValueError(
                f"targets line {lineno}: start/end must be integers"
            ) from exc
        if start < 1 or end < start:
            raise ValueError(
                f"targets line {lineno}: invalid region {start}..{end}"
            )
        targets.append(Target(name=row[name_i].strip(), start=start, end=end))
    if not targets:
        raise ValueError("no targets found in targets file")
    return targets


def pairs_to_rows(pairs: list[PrimerPair]) -> list[dict]:
    """Convert designed pairs to TSV/JSON-ready row dicts."""
    rows: list[dict] = []
    for rank, pair in enumerate(pairs, start=1):
        fwd, rev = pair.forward, pair.reverse
        rows.append(
            {
                "target": pair.target,
                "rank": rank,
                "forward_seq": fwd.seq,
                "forward_start_0based": fwd.start,
                "forward_tm_nn_c": round(fwd.tm_nn, 2),
                "forward_tm_wallace_c": round(fwd.tm_wallace, 2),
                "forward_gc_percent": round(fwd.gc_percent, 2),
                "forward_score": round(fwd.score, 2),
                "reverse_seq": rev.seq,
                "reverse_start_0based": rev.start,
                "reverse_tm_nn_c": round(rev.tm_nn, 2),
                "reverse_tm_wallace_c": round(rev.tm_wallace, 2),
                "reverse_gc_percent": round(rev.gc_percent, 2),
                "reverse_score": round(rev.score, 2),
                "delta_tm_c": round(pair.delta_tm, 2),
                "anneal_temp_c": round(pair.anneal_temp, 2),
                "product_length": pair.product_length,
                "product_seq": pair.product_seq,
                "pair_score": round(pair.score, 2),
                "penalties": "; ".join(pair.penalties),
            }
        )
    return rows


def write_pairs_tsv(pairs: list[PrimerPair], path: str | Path) -> Path:
    """Write designed pairs to a TSV file; return the path written."""
    path = Path(path)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PAIR_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(pairs_to_rows(pairs))
    return path


def write_pairs_json(pairs: list[PrimerPair], path: str | Path) -> Path:
    """Write designed pairs to a JSON file; return the path written."""
    path = Path(path)
    data = [pair.to_dict() for pair in pairs]
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return path


def write_batch_tsv(
    results: dict[str, list[PrimerPair] | str], path: str | Path
) -> Path:
    """Write batch results (all targets) to one TSV; errors become rows.

    A failed target contributes one row whose ``forward_seq`` column
    holds ``ERROR: <message>`` so the failure stays visible in the
    table instead of silently dropping the target.
    """
    path = Path(path)
    rows: list[dict] = []
    for name, value in results.items():
        if isinstance(value, str):
            rows.append(
                {column: "" for column in PAIR_COLUMNS}
                | {"target": name, "rank": 0, "forward_seq": f"ERROR: {value}"}
            )
        else:
            rows.extend(pairs_to_rows(value))
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PAIR_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_batch_json(
    results: dict[str, list[PrimerPair] | str], path: str | Path
) -> Path:
    """Write batch results (all targets) to one JSON file."""
    path = Path(path)
    data = {
        name: (
            [pair.to_dict() for pair in value]
            if isinstance(value, list)
            else {"error": value}
        )
        for name, value in results.items()
    }
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return path
