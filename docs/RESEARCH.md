# Research notes — primer-designer-lite

How PCR primer design works, what this tool implements, and where every
number comes from. All claims below are traceable to the references at
the end; anything this tool does *not* model is listed under
[Limitations](#limitations).

## 1. What a good PCR primer looks like

A PCR primer is a short synthetic oligonucleotide (typically 18–25 nt)
that anneals to one strand of a DNA template so that DNA polymerase can
extend it. Primer3 — the reference implementation most labs use —
summarises the standard guidelines (Untergasser et al., 2012):

- **Length 18–25 nt**: long enough to be specific, short enough to
  synthesise cheaply and anneal efficiently.
- **Melting temperature (Tm) 55–65 °C**, ideally 58–62 °C, with the two
  primers of a pair within a few degrees of each other so they anneal
  at the same temperature.
- **GC content 40–60 %**: GC pairs (three hydrogen bonds) bind more
  strongly than AT pairs (two), so GC content drives Tm.
- **A G or C at the 3′ end ("GC clamp")**: the 3′ end is where the
  polymerase starts extension, so a strongly bound 3′ end improves
  efficiency. 1–3 G/C in the last five bases; more over-stabilises the
  3′ end and invites mispriming.
- **No long homopolymer run** (4+ identical bases in a row): runs risk
  polymerase slippage and synthesis errors.
- **No strong self-complementarity, hairpin, or primer-dimer**: a
  primer that folds on itself or dimerises with its partner is
  unavailable for the template — and a dimer with a paired 3′ end can
  be *extended* by the polymerase, consuming reagents.
- **Unique binding site**: the primer (especially its 3′ end) must
  occur once in the template.

## 2. Melting temperature

### 2.1 Wallace rule (quick estimate)

For short oligonucleotides, Wallace et al. (1979) proposed counting
2 °C per A/T pair and 4 °C per G/C pair:

    Tm = 2·(A+T) + 4·(G+C)

It ignores base order entirely — two primers with the same composition
get the same Tm — so it is only a rule of thumb. Example, hand-checked:
`CGTTCCAAAGATGTGGGCATGAGCTTAC` has 14 A/T and 14 G/C, giving
2·14 + 4·14 = **84.00 °C**, matching the Biopython documentation example.

### 2.2 Nearest-neighbour model (SantaLucia 1998) — the scored estimator

The duplex stability of DNA depends mostly on *stacking* between
adjacent base pairs, not just on base counts. SantaLucia (1998) unified
the thermodynamic parameters of ten nearest-neighbour (NN) doublets,
measured on 108 oligonucleotide duplexes, into a single table
(ΔH° in kcal/mol, ΔS° in cal/(mol·K), at 1 M NaCl):

| Doublet | ΔH° | ΔS° | Doublet | ΔH° | ΔS° |
|---------|-----|-----|---------|-----|-----|
| AA/TT | −7.9 | −22.2 | GA/CT | −8.2 | −22.2 |
| AT/AT | −7.2 | −20.4 | CG/CG | −10.6 | −27.2 |
| TA/TA | −7.2 | −21.3 | GC/GC | −9.8 | −24.4 |
| CA/GT | −8.5 | −22.7 | GG/CC | −8.0 | −19.9 |
| GT/AC | −8.4 | −22.4 | CT/AG | −7.8 | −21.0 |

Initiation terms: a terminal A·T pair adds (2.3, 4.1); a terminal G·C
pair adds (0.1, −2.8). A self-complementary sequence takes a symmetry
correction of ΔS = −1.4 cal/(mol·K). The salt correction (monovalent
cations only) is applied to entropy:

    ΔS[Na+] = ΔS[1 M] + 0.368 · (N − 1) · ln[Na+]

and the melting temperature follows from the two-state model:

    Tm = 1000 · ΔH° / (ΔS[Na+] + R · ln(C_eff)) − 273.15

with R = 1.987 cal/(mol·K) and C_eff the effective strand
concentration in molar.

**Concentration convention.** This tool uses C_eff = primer − template/2
(the convention of Biopython's `Tm_NN`), with defaults 250 nM primer,
0 nM template (primer in excess, as in PCR), 50 mM Na+. The choice of
convention moves Tm by a fraction of a degree; it is documented here
so results are reproducible.

### 2.3 Validation of the implementation

| Check | Expected | This tool |
|---|---|---|
| Biopython docs example `CGTTCCAAAGATGTGGGCATGAGCTTAC`, Tm_Wallace | 84.00 °C | 84.00 °C |
| Same sequence, Tm_NN, Biopython defaults (25 nM / 25 nM, 50 mM Na+) | 60.32 °C | 60.32 °C |
| Same sequence, Tm_NN, 250 nM primer, 50 mM Na+ (Biopython cross-check) | 63.32 °C | 63.32 °C |
| 27F `AGAGTTTGATCATGGCTCAG`, Tm_NN, 250 nM, 50 mM Na+ | 53.48 °C | 53.48 °C |
| 1492R `GGTTACCTTGTTACGACTT`, Tm_NN, 250 nM, 50 mM Na+ | 50.80 °C | 50.80 °C |
| `GCGCTCGTCGTCT`, Tm_NN, 250 nM, 50 mM Na+ | 52.07 °C | 52.07 °C |
| `ATGCATGCATGCATGCATGC`, Tm_NN, 250 nM, 50 mM Na+ | 59.59 °C | 59.59 °C |

The Biopython cross-check values were computed with Biopython 1.8x
`Bio.SeqUtils.MeltingTemp.Tm_NN` (an independent implementation of the
same model) on 2026-11-16. The "expected" column for the first two rows
is published in the Biopython documentation itself.

**Known convention difference.** Owczarzy et al. (2008) give a worked
example: `AAGGCGAGTCAGGCTCAGTG` has Tm **76.3 °C at 1 M Na+**. Under
this tool's (Biopython's) concentration convention the same sequence
gives **75.98 °C** — the ~0.3 °C gap is the concentration-convention
difference, not a parameter error. The value 75.98 °C is pinned in
`tests/test_thermo.py` with this explanation.

### 2.4 Basic GC formula (comparison only)

The textbook estimate Tm = 64.9 + 41·(n_GC − 16.4)/n (the form
implemented as `Tm_GC` in Biopython) is included as a third estimator
for comparison. Like Wallace, it ignores base order. It is *not* used
for scoring.

## 3. Primer-dimer, hairpin, and specificity checks

- **Self-complementarity / dimer**: the longest run of paired bases
  between a primer and the reverse complement of another primer (or
  itself), over all ungapped alignments. Two strands pair exactly where
  the first equals the reverse complement of the second. The
  **3′-anchored** run — a run that includes the primer's 3′ end — is
  scored separately, because only a paired 3′ end can be extended.
- **Hairpin**: a stem of *k* bp pairs the 5′ arm `seq[i:i+k]`
  antiparallel with a downstream arm, leaving a loop of ≥ 3 nt.
  The scan is exhaustive at primer length.
- **Specificity**: exact binding sites of the primer (both strands)
  in the template, near sites (≤ 2 mismatches), and sites where the
  3′-most 12 nt ("seed") match exactly — the seed matters because
  extension only needs the 3′ end to be paired.

## 4. The example template: E. coli 16S rRNA

The bundled example template is the *Escherichia coli* 16S rRNA gene
(rrnB), extracted from GenBank record **J01695.1**, feature
`rRNA 1518..3059` — **1542 bp, GC 54.41 %** — on 2026-11-16
(`examples/ecoli_16s_rrnB_J01695.fasta`). The 16S rRNA gene is the
standard bacterial identification marker, and its universal primer
sites are published and position-numbered, which makes it an ideal
reality check for a designer:

| Primer | Sequence (5′→3′) | E. coli position | Source |
|---|---|---|---|
| 27F | `AGAGTTTGATCMTGGCTCAG` (M = A/C) | 8–27 | Lane (1991); Weisburg et al. (1991) |
| 1492R | `GGTTACCTTGTTACGACTT` | 1492–1510 | Weisburg et al. (1991); Turner et al. (1999) |
| 515F | `GTGCCAGCMGCCGCGGTAA` | 515–536 | Caporaso et al. (2011) |
| 806R | `GGACTACHVGGGTWTCTAAT` | 787–806 | Caporaso et al. (2011) |

Verified in the bundled sequence on 2026-11-16: the E. coli
(non-degenerate) 27F `AGAGTTTGATCATGGCTCAG` occurs at 8–27, 1492R's
reverse complement at 1492–1510, 515F's site `GTGCCAGCAGCCGCGGTAATAC`
at 515–536, and 806R's reverse complement `ATTAGATACCCTGGTAGTCC` at
787–806. The designer recovers these sites in its candidate pool
(`tests/test_design.py::test_design_recovers_universal_sites`), and
its top full-length pair lands on the same sites (see README).

The second example template (`examples/synthetic_example.fasta`) is a
**synthetic** 600 bp sequence generated with `random.seed(20261116)`,
clearly labelled synthetic in its header — it is not a real gene.

## 5. Limitations

- **Monovalent salt only.** The NN model here corrects for Na+ (or
  other monovalent ions treated as Na+). Real PCR buffers also contain
  Mg²⁺, which stabilises duplexes far more strongly; the Owczarzy et
  al. (2008) mixed-salt correction for Mg²⁺ is *not* implemented, so
  reported Tm values are systematically a few degrees below the Tm in
  a Mg²⁺-containing buffer. Treat the annealing temperature as a
  starting point for a gradient, not a final protocol value.
- **No mismatch or dangling-end thermodynamics.** Tm is computed for a
  perfectly matched duplex only.
- **Degenerate bases (N, R, Y, …) are rejected.** Universal primers
  with IUPAC degeneracy (e.g. 27F's M) cannot be scored directly; use
  the concrete E. coli variant instead.
- **Ungapped dimer/hairpin scans only.** No thermodynamic ΔG for
  dimers/hairpins (as Primer3 computes); the checks report the longest
  complementary run, a conservative proxy.
- **Specificity is template-local.** Off-targets elsewhere in a genome
  need a BLAST-style search; this tool only checks the given template.
- **No experimental validation.** All "good primer" judgements here
  are computational; a designed pair still needs a real PCR test.

## References

1. SantaLucia J. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *Proceedings of the National Academy of Sciences* 95(4):1460–1465. https://www.pnas.org/doi/10.1073/pnas.95.4.1460
2. Wallace RB, Shaffer J, Murphy RF, Bonner J, Hirose T, Itakura K. (1979). Hybridization of synthetic oligodeoxyribonucleotides to ΦX174 DNA: the effect of single base pair mismatch. *Nucleic Acids Research* 6(11):3543–3557. https://doi.org/10.1093/nar/6.11.3543
3. Biopython — `Bio.SeqUtils.MeltingTemp` module documentation (Tm_Wallace / Tm_GC / Tm_NN, worked example 60.32 °C). https://biopython.org/docs/1.79/api/Bio.SeqUtils.MeltingTemp.html
4. Owczarzy R, You Y, Moreira BG, Manthey JA, Huang L, Behlke MA, Walder JA. (2008). Effects of sodium ions on DNA duplex stability and kinetics. *Biochemistry* 47(19):5336–5353. https://doi.org/10.1021/bi702363u
5. Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG. (2012). Primer3 — new capabilities and interfaces. *Nucleic Acids Research* 40(15):e115. https://doi.org/10.1093/nar/gks596
6. Weisburg WG, Barns SM, Pelletier DA, Lane DJ. (1991). 16S ribosomal DNA amplification for phylogenetic study. *Journal of Bacteriology* 173(2):697–703. https://doi.org/10.1128/jb.173.2.697-703.1991
7. Lane DJ. (1991). 16S/23S rRNA sequencing. In *Nucleic Acid Techniques in Bacterial Systematics* (Stackebrandt E, Goodfellow M, eds.), Wiley, New York, pp. 115–175.
8. Caporaso JG, Lauber CL, Walters WA, Berg-Lyons D, Lozupone CA, Turnbaugh PJ, Fierer N, Knight R. (2011). Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample. *Proceedings of the National Academy of Sciences* 108(Suppl 1):4516–4522. https://doi.org/10.1073/pnas.1000080107
9. NCBI Nucleotide — *Escherichia coli* rRNA operon (rrnB), GenBank J01695.1 (16S rRNA feature 1518..3059). https://www.ncbi.nlm.nih.gov/nuccore/J01695.1
10. von Ahsen N, Wittwer CT, Schütz E. (2001). Oligonucleotide melting temperatures under PCR conditions: nearest-neighbor corrections for Mg²⁺, deoxynucleotide triphosphate, and dimethyl sulfoxide concentrations with comparison to alternative empirical formulas. *Clinical Chemistry* 47(11):1956–1961.
