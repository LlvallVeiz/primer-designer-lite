# Design decisions — primer-designer-lite

## Coordinate conventions

- **Command line and TSV files use 1-based inclusive coordinates**
  (`--start 27 --end 1492`), the convention of primer names (27F binds
  at 8–27) and of GenBank feature tables. Internally everything is
  0-based half-open `[start, end)`.
- **The target region is the stretch between the primers' 3′ ends.**
  The forward primer's 3′ end lands near `start`, the reverse primer's
  3′ end near `end`. The amplicon (product) spans from the forward
  primer's 5′ end to the reverse primer's 5′ end on the template:
  `product = template[forward.start : reverse.end]`. For the full-length
  16S target (27..1492) that gives a ~1499 bp product, matching the
  published ~1.5 kb full-length 16S amplicon.
- **Primer positions are reported 0-based** (`forward_start_0based`)
  in TSV/JSON output so they can be sliced directly in Python.

## Candidate generation

For each side of the target, candidates slide their 3′ end by up to
`--shift` (default 10) bases either side of the boundary and vary in
length from `--min-len` to `--max-len` (default 18–25). That yields
about 8 lengths × 21 shifts ≈ 168 candidates per side, each evaluated
and scored; the best `--top-candidates` (default 30) per side enter
pairing. Sliding the 3′ end — rather than fixing it exactly at the
boundary — is what lets the designer find a primer with a good Tm and
GC clamp a few bases away from the nominal boundary.

## Scoring (additive penalties; lower is better, 0.0 is ideal)

Individual primer (see `primer.evaluate_primer` for the exact code):

| Rule | Penalty |
|---|---|
| NN Tm outside ideal 58–62 °C | 1.0 per °C outside the ideal window |
| NN Tm outside accepted 55–65 °C | 10.0 + 2.0 per °C beyond the range |
| GC% outside 40–60 | 5.0 + 0.5 per point beyond the range |
| Length ≠ 20 nt | 0.2 per nt away from 20 |
| Homopolymer run of exactly 4 | 6.0 |
| No G/C at the 3′ end | 4.0 |
| 0 or ≥ 4 G/C in the last 5 bases | 3.0 |
| 3′-anchored self-complementarity of 4 | 5.0 |
| Any self-complementary run ≥ 6 | 3.0 |
| Hairpin stem of 5–6 bp | 4.0 |
| Each extra exact binding site | 10.0 |
| Each extra 3′-seed binding site | 5.0 |

Hard rejections (the candidate never appears in output):

- Homopolymer run of 5+.
- 3′-anchored self-complementarity of 5+.
- Hairpin stem of 7+ bp.

Pair scoring (see `design.design_pair`):

| Rule | Penalty |
|---|---|
| Forward/reverse ΔTm | 3.0 per °C |
| ΔTm > 3.0 °C (default `--max-delta-tm`) | pair rejected |
| Cross-dimer 3′ run of 4 bp | 8.0 |
| Cross-dimer 3′ run of 5+ bp | pair rejected |
| Product outside 150–500 bp | 0.01 per bp, capped at 5.0 |
| Product outside `--min-product`–`--max-product` | pair rejected |

The recommended annealing temperature is the lower primer Tm minus
5 °C — the standard starting point for a gradient optimisation, not a
final protocol value.

## Why Tm NN is scored and Wallace is not

Wallace ignores base order, so it cannot distinguish two primers with
the same composition; it is reported for comparison only. The NN model
is order-sensitive (two same-composition primers can differ by several
degrees) and is the model Primer3 and Biopython use, so it is the
scored estimator. The basic GC formula is reported for the same
comparison purpose.

## Why no dependencies

The whole tool is standard-library Python. The NN tables are small
constants, the scans are O(n) to O(n²) at primer length, and the only
O(template × primer) work is the specificity scan — fast enough at
kilobase templates without numpy. That keeps installation to
`pip install -r requirements.txt` (pytest only) and makes the tool
usable on any machine with Python 3.9+.

## Determinism

No randomness anywhere: candidates are generated in a fixed order and
sorted by (score, Tm NN, sequence); ties break identically on every
run. The synthetic example template is generated from a fixed seed
(`random.seed(20261116)`) and stored as a file, so examples reproduce
bit-for-bit.

## What was deliberately left out

- **Mg²⁺ correction** (Owczarzy 2008): would need a second salt model
  and dNTP bookkeeping; see RESEARCH.md limitations.
- **Thermodynamic ΔG for dimers/hairpins** (Primer3-style): the
  longest-complementary-run proxy is conservative and explainable
  without a full folding model.
- **Genome-wide specificity (BLAST)**: out of scope for a lite tool;
  specificity is checked against the given template only.
- **Degenerate primers**: NN parameters do not cover IUPAC codes.
