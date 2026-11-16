"""Command-line interface for primer-designer-lite.

Subcommands
-----------
``tm``:
    Print Wallace, GC-basic, and nearest-neighbour Tm plus GC% for one
    primer sequence.
``check``:
    Evaluate a forward/reverse primer pair: individual scores and
    checks, cross-dimer report, and Tm difference.
``design``:
    Design pairs for one target region on a FASTA template.
``batch``:
    Design pairs for many targets listed in a TSV file.

Coordinates on the command line are 1-based inclusive (``--start`` /
``--end``), matching the targets TSV format and the primer-naming
convention (27F binds at 8-27, 1492R at 1492-1510 of the *E. coli*
16S rRNA gene). The product spans from the forward primer's 5' end
to the reverse primer's 5' end on the template.
"""

from __future__ import annotations

import argparse
import sys

from .checks import dimer_report, self_report
from .design import PrimerPair, design_batch, design_pair
from .io import (
    read_fasta,
    read_targets,
    write_batch_json,
    write_batch_tsv,
    write_pairs_json,
    write_pairs_tsv,
)
from .primer import DesignConfig, evaluate_primer
from .thermo import gc_content, tm_gc_basic, tm_nn, tm_wallace, validate_sequence

__all__ = ["main", "build_parser"]


def _add_design_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--min-len", type=int, default=18)
    parser.add_argument("--max-len", type=int, default=25)
    parser.add_argument("--tm-min", type=float, default=55.0)
    parser.add_argument("--tm-max", type=float, default=65.0)
    parser.add_argument("--gc-min", type=float, default=40.0)
    parser.add_argument("--gc-max", type=float, default=60.0)
    parser.add_argument("--na", type=float, default=50.0,
                        help="Na+ concentration in mM (default 50)")
    parser.add_argument("--primer-conc", type=float, default=250.0,
                        help="primer concentration in nM (default 250)")
    parser.add_argument("--flank", type=int, default=300)
    parser.add_argument("--shift", type=int, default=10)
    parser.add_argument("--min-product", type=int, default=70)
    parser.add_argument("--max-product", type=int, default=1500)
    parser.add_argument("--max-delta-tm", type=float, default=3.0)
    parser.add_argument("--top", type=int, default=5,
                        help="number of pairs to report (default 5)")
    parser.add_argument("--out", help="TSV output path (optional)")
    parser.add_argument("--json", dest="json_out",
                        help="JSON output path (optional)")


def _config_from_args(args: argparse.Namespace) -> DesignConfig:
    return DesignConfig(
        min_len=args.min_len,
        max_len=args.max_len,
        tm_min=args.tm_min,
        tm_max=args.tm_max,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        na_mM=args.na,
        primer_nM=args.primer_conc,
        flank=args.flank,
        shift=args.shift,
        min_product=args.min_product,
        max_product=args.max_product,
        max_delta_tm=args.max_delta_tm,
        top_pairs=args.top,
    )


def _print_pair(pair: PrimerPair, rank: int) -> None:
    fwd, rev = pair.forward, pair.reverse
    print(f"rank {rank}  score {pair.score:.2f}  product {pair.product_length} bp")
    print(f"  forward  {fwd.seq}  ({fwd.length} nt, start {fwd.start} 0-based)")
    print(
        f"           Tm NN {fwd.tm_nn:.2f}C  Wallace {fwd.tm_wallace:.2f}C  "
        f"GC {fwd.gc_percent:.1f}%  score {fwd.score:.2f}"
    )
    if fwd.penalties:
        print(f"           penalties: {'; '.join(fwd.penalties)}")
    print(f"  reverse  {rev.seq}  ({rev.length} nt, start {rev.start} 0-based)")
    print(
        f"           Tm NN {rev.tm_nn:.2f}C  Wallace {rev.tm_wallace:.2f}C  "
        f"GC {rev.gc_percent:.1f}%  score {rev.score:.2f}"
    )
    if rev.penalties:
        print(f"           penalties: {'; '.join(rev.penalties)}")
    print(
        f"  delta Tm {pair.delta_tm:.2f}C  anneal ~{pair.anneal_temp:.2f}C"
    )
    if pair.penalties:
        print(f"  pair penalties: {'; '.join(pair.penalties)}")


def _cmd_tm(args: argparse.Namespace) -> int:
    seq = validate_sequence(args.sequence)
    print(f"sequence  {seq}  ({len(seq)} nt)")
    print(f"Tm Wallace  {tm_wallace(seq):.2f} C")
    print(f"Tm GC basic {tm_gc_basic(seq):.2f} C")
    print(
        f"Tm NN       {tm_nn(seq, na_mM=args.na, primer_nM=args.primer_conc):.2f} C"
        f"  (SantaLucia 1998, Na+ {args.na:g} mM, primer {args.primer_conc:g} nM)"
    )
    print(f"GC content  {gc_content(seq):.2f} %")
    return 0


def _cmd_check(args: argparse.Namespace) -> int:
    fwd_seq = validate_sequence(args.forward)
    rev_seq = validate_sequence(args.reverse)
    template = ""
    if args.template:
        records = read_fasta(args.template)
        template = records[0][1]
    config = DesignConfig(na_mM=args.na, primer_nM=args.primer_conc)
    for label, seq, start in (
        ("forward", fwd_seq, 0),
        ("reverse", rev_seq, 0),
    ):
        primer = evaluate_primer(seq, start, "+" if label == "forward" else "-", template, config) if template else None
        print(f"{label}  {seq}  ({len(seq)} nt)")
        print(f"  Tm NN {tm_nn(seq, na_mM=args.na, primer_nM=args.primer_conc):.2f} C"
              f"  Wallace {tm_wallace(seq):.2f} C  GC {gc_content(seq):.1f}%")
        report = self_report(seq)
        print(f"  homopolymer run {report['homopolymer_run']}, "
              f"GC clamp {'yes' if report['gc_clamp']['has_clamp'] else 'no'} "
              f"({report['gc_clamp']['gc_in_last_5']} GC in last 5), "
              f"self-complementarity {report['self_complementarity']}, "
              f"3' self {report['three_prime_self']}, "
              f"hairpin stem {report['hairpin_stem']} bp")
        if primer is not None:
            status = "OK" if primer.score == 0.0 else f"score {primer.score:.2f}"
            print(f"  evaluation: {status}")
            if primer.penalties:
                print(f"  penalties: {'; '.join(primer.penalties)}")
    dimer = dimer_report(fwd_seq, rev_seq)
    delta_tm = abs(
        tm_nn(fwd_seq, na_mM=args.na, primer_nM=args.primer_conc)
        - tm_nn(rev_seq, na_mM=args.na, primer_nM=args.primer_conc)
    )
    print(f"cross-dimer  max run {dimer['max_complementarity']} bp, "
          f"3' forward {dimer['three_prime_a']} bp, "
          f"3' reverse {dimer['three_prime_b']} bp")
    print(f"delta Tm     {delta_tm:.2f} C")
    return 0


def _cmd_design(args: argparse.Namespace) -> int:
    records = read_fasta(args.fasta)
    if len(records) > 1 and not args.record:
        print(
            f"note: {args.fasta} has {len(records)} records; "
            "using the first (pass --record N to pick another)",
            file=sys.stderr,
        )
    index = args.record or 0
    header, template = records[index]
    config = _config_from_args(args)
    start, end = args.start - 1, args.end  # to 0-based half-open
    pairs = design_pair(
        template, start, end, args.name or "target", config
    )
    print(f"template  {header}  ({len(template)} bp)")
    print(f"target    {args.name or 'target'}  {args.start}..{args.end} (1-based)")
    for rank, pair in enumerate(pairs, start=1):
        _print_pair(pair, rank)
    if args.out:
        write_pairs_tsv(pairs, args.out)
        print(f"wrote {args.out}")
    if args.json_out:
        write_pairs_json(pairs, args.json_out)
        print(f"wrote {args.json_out}")
    return 0


def _cmd_batch(args: argparse.Namespace) -> int:
    records = read_fasta(args.fasta)
    index = args.record or 0
    header, template = records[index]
    targets = read_targets(args.targets)
    config = _config_from_args(args)
    results = design_batch(template, targets, config)
    print(f"template  {header}  ({len(template)} bp)")
    for name, value in results.items():
        if isinstance(value, str):
            print(f"{name}: ERROR: {value}")
            continue
        print(f"{name}: {len(value)} pair(s), best score {value[0].score:.2f}, "
              f"product {value[0].product_length} bp")
        for rank, pair in enumerate(value, start=1):
            _print_pair(pair, rank)
    if args.out:
        write_batch_tsv(results, args.out)
        print(f"wrote {args.out}")
    if args.json_out:
        write_batch_json(results, args.json_out)
        print(f"wrote {args.json_out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="primer-designer",
        description="Design PCR primer pairs from a DNA template (FASTA).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_tm = sub.add_parser("tm", help="Tm and GC% for one primer sequence")
    p_tm.add_argument("sequence")
    p_tm.add_argument("--na", type=float, default=50.0)
    p_tm.add_argument("--primer-conc", type=float, default=250.0)
    p_tm.set_defaults(func=_cmd_tm)

    p_check = sub.add_parser("check", help="evaluate a primer pair")
    p_check.add_argument("--forward", required=True)
    p_check.add_argument("--reverse", required=True)
    p_check.add_argument("--template", help="FASTA template for specificity")
    p_check.add_argument("--na", type=float, default=50.0)
    p_check.add_argument("--primer-conc", type=float, default=250.0)
    p_check.set_defaults(func=_cmd_check)

    p_design = sub.add_parser("design", help="design pairs for one target")
    p_design.add_argument("--fasta", required=True)
    p_design.add_argument("--start", type=int, required=True,
                          help="target start, 1-based inclusive")
    p_design.add_argument("--end", type=int, required=True,
                          help="target end, 1-based inclusive")
    p_design.add_argument("--name", default="target")
    p_design.add_argument("--record", type=int, default=0,
                          help="FASTA record index (default 0)")
    _add_design_options(p_design)
    p_design.set_defaults(func=_cmd_design)

    p_batch = sub.add_parser("batch", help="design pairs for targets in a TSV")
    p_batch.add_argument("--fasta", required=True)
    p_batch.add_argument("--targets", required=True,
                         help="TSV with name/start/end (1-based inclusive)")
    p_batch.add_argument("--record", type=int, default=0)
    _add_design_options(p_batch)
    p_batch.set_defaults(func=_cmd_batch)

    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entry point; returns a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, FileNotFoundError, IndexError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
