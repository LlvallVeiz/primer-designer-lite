# primer-designer-lite

Design PCR primer pairs from a DNA template (FASTA) — nearest-neighbour
Tm (SantaLucia 1998), quality checks (homopolymer, GC clamp,
self-complementarity, hairpin, primer-dimer), template specificity, and
batch mode. Standard-library Python only; no dependencies at runtime.

## Background

A PCR primer is a short oligonucleotide (18–25 nt) that anneals to a
template so DNA polymerase can extend it. A good primer pair needs
matched melting temperatures (ideally 58–62 °C, within a few degrees of
each other), 40–60 % GC, a G/C at the 3′ end, no long homopolymer run,
no strong self-complementarity or hairpin, no extendable primer-dimer,
and a unique binding site in the template. This tool scores candidates
against exactly those rules (Primer3 guidelines; Untergasser et al.,
2012) and reports the best pairs with their penalties, so you can see
*why* a primer was chosen — or rejected.

Two Tm estimators are reported for every primer:

- **Wallace rule** — 2 °C per A/T + 4 °C per G/C. Fast, but ignores
  base order.
- **Nearest-neighbour (NN) model** — SantaLucia (1998) unified
  parameters with the monovalent-salt correction. This is the scored
  estimator. Validated against Biopython's `MeltingTemp.Tm_NN`: the
  documentation example `CGTTCCAAAGATGTGGGCATGAGCTTAC` gives Tm_NN =
  **60.32 °C** with Biopython's defaults, reproduced exactly.

See [`docs/RESEARCH.md`](docs/RESEARCH.md) for the full background,
the parameter tables, the validation table, and the limitations
(no Mg²⁺ correction, no degenerate bases, template-local specificity
only). See [`docs/DESIGN.md`](docs/DESIGN.md) for the scoring system
and coordinate conventions.

## Installation

Python 3.9+; no runtime dependencies.

```bash
git clone https://github.com/LlvallVeiz/primer-designer-lite
cd primer-designer-lite
pip install -r requirements.txt   # pytest only, for running the tests
```

## Usage

### `tm` — melting temperature of one primer

```bash
python -m primer_designer tm AGAGTTTGATCATGGCTCAG
```

```text
sequence  AGAGTTTGATCATGGCTCAG  (20 nt)
Tm Wallace  58.00 C
Tm GC basic 49.73 C
Tm NN       53.48 C  (SantaLucia 1998, Na+ 50 mM, primer 250 nM)
GC content  45.00 %
```

(That sequence is the universal bacterial primer 27F.)

### `check` — evaluate a primer pair you already have

```bash
python -m primer_designer check \
  --forward AGAGTTTGATCATGGCTCAG \
  --reverse GGTTACCTTGTTACGACTT \
  --template examples/ecoli_16s_rrnB_J01695.fasta
```

```text
forward  AGAGTTTGATCATGGCTCAG  (20 nt)
  Tm NN 53.48 C  Wallace 58.00 C  GC 45.0%
  homopolymer run 3, GC clamp yes (3 GC in last 5), self-complementarity 6, 3' self 2, hairpin stem 3 bp
  evaluation: score 16.03
  penalties: Tm NN 53.5C below 55.0C; self-complementarity run of 6
reverse  GGTTACCTTGTTACGACTT  (19 nt)
  Tm NN 50.80 C  Wallace 54.00 C  GC 42.1%
  homopolymer run 2, GC clamp no (2 GC in last 5), self-complementarity 3, 3' self 1, hairpin stem 2 bp
  evaluation: score 22.59
  penalties: Tm NN 50.8C below 55.0C; no GC clamp (3' base T)
cross-dimer  max run 3 bp, 3' forward 2 bp, 3' reverse 1 bp
delta Tm     2.68 C
```

Note what this shows: the classic 27F/1492R pair — designed decades
before NN models — scores poorly by modern criteria (both Tm NN values
sit below the 55–65 °C window at 50 mM Na+ / 250 nM primer), yet it is
dimer-clean and works in practice. The tool reports the facts; the
judgement about whether to use a published primer anyway is yours.

### `design` — design pairs for one target region

Coordinates are **1-based inclusive**, as in primer names (27F binds at
8–27 of the *E. coli* 16S gene) and GenBank features. The target is the
stretch between the primers' 3′ ends; the product spans from the
forward primer's 5′ end to the reverse primer's 5′ end.

```bash
python -m primer_designer design \
  --fasta examples/ecoli_16s_rrnB_J01695.fasta \
  --start 27 --end 1492 --name 16s_full \
  --out examples/output_16s_full.tsv --json examples/output_16s_full.json
```

```text
template  ecoli_16s_rrnB_J01695 Escherichia coli 16S rRNA gene (rrnB), GenBank J01695.1:1518..3059, 1542 bp  (1542 bp)
target    16s_full  27..1492 (1-based)
rank 1  score 7.10  product 1499 bp
  forward  GATCATGGCTCAGATTGAACGC  (22 nt, start 14 0-based)
           Tm NN 58.00C  Wallace 66.00C  GC 50.0%  score 0.40
           penalties: Tm NN 58.0C below ideal 58C
  reverse  TACGGTTACCTTGTTACGACTTCAC  (25 nt, start 1488 0-based)
           Tm NN 58.23C  Wallace 72.00C  GC 44.0%  score 1.00
  delta Tm 0.23C  anneal ~53.00C
  pair penalties: product 1499 bp outside ideal 150-500 bp
```

The top pair lands on the published universal sites: the forward
primer overlaps 27F (its 3′ end is at the same boundary) and the
reverse primer contains the full 1492R sequence
(`GGTTACCTTGTTACGACTT`) extended by 6 bases. The ~1499 bp product
matches the published ~1.5 kb full-length 16S amplicon.

The same command on the V4 region (`--start 536 --end 787`) gives a
300 bp product with ΔTm 0.14 °C:

```text
rank 1  score 8.81  product 300 bp
  forward  AGCAGCCGCGGTAATACGGA  (20 nt, start 519 0-based)
           Tm NN 62.47C  Wallace 64.00C  GC 60.0%  score 7.47
           penalties: Tm NN 62.5C above ideal 62C; no GC clamp (3' base A); self-complementarity run of 6
  reverse  TCGTTTACGGCGTGGACTACCAG  (23 nt, start 796 0-based)
           Tm NN 62.33C  Wallace 72.00C  GC 56.5%  score 0.93
           penalties: Tm NN 62.3C above ideal 62C
  delta Tm 0.14C  anneal ~57.33C
```

(The forward primer sits on the 515F site; its "no GC clamp" penalty
is real — the best V4 forward primer the tool found ends in A — and
is reported rather than hidden.)

### `batch` — design pairs for many targets at once

```bash
python -m primer_designer batch \
  --fasta examples/ecoli_16s_rrnB_J01695.fasta \
  --targets examples/targets_batch.tsv \
  --top 2 --out examples/output_batch.tsv --json examples/output_batch.json
```

```text
template  ecoli_16s_rrnB_J01695 Escherichia coli 16S rRNA gene (rrnB), GenBank J01695.1:1518..3059, 1542 bp  (1542 bp)
16s_full: 2 pair(s), best score 7.10, product 1499 bp
16s_v4: 2 pair(s), best score 8.81, product 300 bp
```

Targets TSV format (1-based inclusive, header row, `#` comments
allowed):

```tsv
name	start	end
16s_full	27	1492
16s_v4	536	787
```

A target that cannot be designed does not abort the batch: it is
reported as `ERROR: ...` in the TSV/JSON output.

### Options

```
--min-len / --max-len     primer length range (default 18-25)
--tm-min / --tm-max       accepted NN Tm window (default 55-65 C)
--gc-min / --gc-max       accepted GC% window (default 40-60)
--na                      Na+ in mM (default 50)
--primer-conc             primer concentration in nM (default 250)
--flank / --shift         search window around the target boundary
--min-product / --max-product   amplicon size window (default 70-1500)
--max-delta-tm            max forward/reverse Tm difference (default 3.0)
--top                     pairs to report (default 5)
--out / --json            TSV / JSON output paths
```

## Examples

| File | Contents |
|---|---|
| `examples/ecoli_16s_rrnB_J01695.fasta` | *E. coli* 16S rRNA gene (rrnB), GenBank J01695.1:1518..3059, 1542 bp, GC 54.41 % |
| `examples/synthetic_example.fasta` | Synthetic 600 bp template (`random.seed(20261116)`), clearly labelled synthetic — not a real gene |
| `examples/targets_16s.tsv` | 16S targets: full-length (27..1492) and V4 (536..787) |
| `examples/targets_synthetic.tsv` | Synthetic target (201..400) |
| `examples/targets_batch.tsv` | Both 16S targets for batch mode |
| `examples/output_16s_full.tsv` / `.json` | Real `design` output, full-length 16S (best pair: Tm 58.00/58.23 °C, product 1499 bp) |
| `examples/output_16s_v4.tsv` | Real `design` output, V4 (best pair: Tm 62.47/62.33 °C, product 300 bp) |
| `examples/output_synthetic.tsv` | Real `design` output, synthetic target (best pair: Tm 58.82/58.75 °C, product 245 bp, ΔTm 0.07 °C) |
| `examples/output_batch.tsv` / `.json` | Real `batch` output for both 16S targets |

## Tests

```bash
python -m pytest -q
```

**78 passed.** The suite pins the implementation against published
values (Biopython documentation Tm_NN = 60.32 °C; cross-checks against
Biopython `MeltingTemp` for 27F, 1492R, and three further sequences),
hand-computed Wallace/GC values, hand-checked hairpin/dimer cases, the
recovery of the published 27F/1492R sites from the real 16S template,
and CLI smoke tests on the real example files.

## Project structure

```text
primer-designer-lite/
├── primer_designer/
│   ├── thermo.py        Tm Wallace / GC-basic / nearest-neighbour, GC%, reverse complement
│   ├── checks.py        homopolymer, GC clamp, self-complementarity, hairpin, dimer
│   ├── specificity.py   binding-site search (exact / near / 3' seed) in the template
│   ├── primer.py        Primer dataclass, penalty scoring, candidate generation
│   ├── design.py        pair design, batch design, product / annealing temperature
│   ├── io.py            FASTA + targets TSV in, pairs TSV / JSON out
│   ├── cli.py           argparse CLI: tm / check / design / batch
│   └── __main__.py      python -m primer_designer
├── tests/               78 pytest tests (thermo, checks, specificity, primer, design, io, cli)
├── examples/            real 16S rRNA template, synthetic template, targets, real outputs
├── docs/
│   ├── RESEARCH.md      background, parameter tables, validation, references
│   └── DESIGN.md        scoring system, coordinate conventions, design decisions
├── pyproject.toml
├── requirements.txt     pytest only (runtime is stdlib-only)
└── LICENSE              MIT
```

## Limitations

- **Monovalent salt only.** Real PCR buffers contain Mg²⁺, which
  raises Tm by a few degrees (Owczarzy et al., 2008 correction not
  implemented). Use the reported annealing temperature as a gradient
  starting point.
- **No degenerate bases** (N, R, Y, …): use the concrete variant of a
  degenerate primer instead.
- **Specificity is template-local**: off-targets elsewhere in a genome
  need a BLAST-style search.
- **Ungapped dimer/hairpin scans**: longest complementary run, not a
  thermodynamic ΔG.
- Designed pairs are computational predictions — validate with a real
  PCR before ordering.

## Future work

- Mg²⁺ / dNTP salt correction (Owczarzy 2008).
- Thermodynamic ΔG for dimers and hairpins (Primer3-style).
- Degenerate-base support for universal primers.
- Primer efficiency prediction from the 3′ end sequence.
- Optional Biopython-backed cross-check mode.

## References

- SantaLucia J. (1998). *PNAS* 95(4):1460–1465. https://www.pnas.org/doi/10.1073/pnas.95.4.1460
- Wallace RB et al. (1979). *Nucleic Acids Research* 6(11):3543–3557. https://doi.org/10.1093/nar/6.11.3543
- Biopython `Bio.SeqUtils.MeltingTemp` documentation. https://biopython.org/docs/1.79/api/Bio.SeqUtils.MeltingTemp.html
- Untergasser A et al. (2012). Primer3. *Nucleic Acids Research* 40(15):e115. https://doi.org/10.1093/nar/gks596
- Weisburg WG et al. (1991). *Journal of Bacteriology* 173(2):697–703. https://doi.org/10.1128/jb.173.2.697-703.1991
- Caporaso JG et al. (2011). *PNAS* 108(Suppl 1):4516–4522. https://doi.org/10.1073/pnas.1000080107
- NCBI Nucleotide, GenBank J01695.1. https://www.ncbi.nlm.nih.gov/nuccore/J01695.1

Full reference list and research notes: [`docs/RESEARCH.md`](docs/RESEARCH.md).

## License

MIT — see [LICENSE](LICENSE).
