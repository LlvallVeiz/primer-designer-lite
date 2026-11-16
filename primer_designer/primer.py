"""Primer candidates: evaluation, penalty scoring, and generation.

A :class:`Primer` is one candidate oligonucleotide anchored at a
position on a template. It is scored by a transparent penalty system
(documented in :func:`evaluate_primer`); lower is better, 0.0 is a
primer that meets every default guideline. Candidates that violate a
hard rule -- a homopolymer run of 5+, a hairpin stem of 7+ bp, or more
than one exact binding site -- are rejected outright instead of merely
penalised, so they never appear in design output.

Design guidelines (defaults in :class:`DesignConfig`,
:mod:`primer_designer.design`) follow the Primer3 manual and
Untergasser et al. (2012): length 18-25 nt, NN Tm 55-65 degC (ideal
58-62), GC 40-60 %, a G/C 3' clamp, no long homopolymer run, and a
small forward/reverse Tm difference.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .checks import (
    dimer_report,
    gc_clamp,
    hairpin_stem,
    longest_homopolymer_run,
    three_prime_complementarity,
)
from .specificity import specificity_report
from .thermo import gc_content, reverse_complement, tm_nn, tm_wallace, validate_sequence

__all__ = ["Primer", "DesignConfig", "evaluate_primer", "generate_candidates"]

# Hard-rejection limits. A candidate breaking one of these is unusable
# regardless of its other properties.
MAX_HOMOPOLYMER_RUN = 4  # runs of 5+ are rejected
MAX_HAIRPIN_STEM = 6  # stems of 7+ bp are rejected
MAX_SELF_3PRIME = 4  # 3'-anchored self-complementarity of 5+ is rejected


@dataclass
class DesignConfig:
    """Tunable design parameters (all defaults documented in the README)."""

    min_len: int = 18
    max_len: int = 25
    tm_min: float = 55.0
    tm_max: float = 65.0
    tm_ideal_low: float = 58.0
    tm_ideal_high: float = 62.0
    gc_min: float = 40.0
    gc_max: float = 60.0
    na_mM: float = 50.0
    primer_nM: float = 250.0
    flank: int = 300
    shift: int = 10
    min_product: int = 70
    max_product: int = 1500
    max_delta_tm: float = 3.0
    top_pairs: int = 5
    top_candidates: int = 30


@dataclass
class Primer:
    """One evaluated primer candidate.

    ``start`` is the 0-based position of the primer's 5' end on the
    template for a forward primer; for a reverse primer it is the
    0-based start of the template segment whose reverse complement is
    the primer (so the primer's 5' end sits at ``start + len(seq)`` on
    the template). ``strand`` is ``"+"`` (forward) or ``"-"`` (reverse).
    """

    seq: str
    start: int
    strand: str
    tm_wallace: float = 0.0
    tm_nn: float = 0.0
    gc_percent: float = 0.0
    score: float = 0.0
    penalties: list[str] = field(default_factory=list)
    checks: dict = field(default_factory=dict)

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def end(self) -> int:
        """0-based exclusive end of the primer footprint on the template."""
        return self.start + self.length

    def to_dict(self) -> dict:
        """JSON-friendly representation (rounded floats)."""
        return {
            "sequence": self.seq,
            "start_0based": self.start,
            "strand": self.strand,
            "length": self.length,
            "tm_wallace_c": round(self.tm_wallace, 2),
            "tm_nn_c": round(self.tm_nn, 2),
            "gc_percent": round(self.gc_percent, 2),
            "score": round(self.score, 2),
            "penalties": list(self.penalties),
            "checks": self.checks,
        }


def evaluate_primer(
    seq: str,
    start: int,
    strand: str,
    template: str,
    config: DesignConfig | None = None,
) -> Primer | None:
    """Evaluate one candidate; return ``None`` when it breaks a hard rule.

    Penalty scheme (additive; lower is better):

    * NN Tm outside the ideal 58-62 degC window: 1.0 per degC outside
      the ideal window; outside the accepted 55-65 degC range: 10.0
      plus 2.0 per degC beyond the range.
    * GC% outside 40-60: 5.0 plus 0.5 per point beyond the range.
    * Length: 0.2 per nt away from 20 (mild preference for 20-mers).
    * Homopolymer run of exactly 4: 6.0. Runs of 5+ are rejected.
    * No G/C at the 3' end: 4.0. Zero or 4+ G/C in the last five
      bases: 3.0.
    * 3'-anchored self-complementarity of 4: 5.0; of 5+: rejected.
      Any self-complementary run of 6+: 3.0.
    * Hairpin stem of 5-6 bp: 4.0; of 7+ bp: rejected.
    * Extra exact binding sites in the template: 10.0 each; extra
      3'-seed sites: 5.0 each.
    """
    config = config or DesignConfig()
    seq = validate_sequence(seq)
    template = validate_sequence(template)
    primer = Primer(seq=seq, start=start, strand=strand)
    primer.tm_wallace = tm_wallace(seq)
    primer.tm_nn = tm_nn(seq, na_mM=config.na_mM, primer_nM=config.primer_nM)
    primer.gc_percent = gc_content(seq)

    penalties: list[str] = []
    score = 0.0

    # --- Tm (nearest neighbour is the scored estimator) ---
    if primer.tm_nn < config.tm_min:
        over = config.tm_min - primer.tm_nn
        score += 10.0 + 2.0 * over
        penalties.append(f"Tm NN {primer.tm_nn:.1f}C below {config.tm_min}C")
    elif primer.tm_nn > config.tm_max:
        over = primer.tm_nn - config.tm_max
        score += 10.0 + 2.0 * over
        penalties.append(f"Tm NN {primer.tm_nn:.1f}C above {config.tm_max}C")
    elif primer.tm_nn < config.tm_ideal_low:
        score += config.tm_ideal_low - primer.tm_nn
        penalties.append(f"Tm NN {primer.tm_nn:.1f}C below ideal 58C")
    elif primer.tm_nn > config.tm_ideal_high:
        score += primer.tm_nn - config.tm_ideal_high
        penalties.append(f"Tm NN {primer.tm_nn:.1f}C above ideal 62C")

    # --- GC content ---
    if primer.gc_percent < config.gc_min:
        score += 5.0 + 0.5 * (config.gc_min - primer.gc_percent)
        penalties.append(f"GC {primer.gc_percent:.1f}% below 40%")
    elif primer.gc_percent > config.gc_max:
        score += 5.0 + 0.5 * (primer.gc_percent - config.gc_max)
        penalties.append(f"GC {primer.gc_percent:.1f}% above 60%")

    # --- Length preference ---
    score += 0.2 * abs(len(seq) - 20)

    # --- Homopolymer run ---
    run = longest_homopolymer_run(seq)
    if run > MAX_HOMOPOLYMER_RUN:
        return None  # hard reject: run of 5+
    if run == MAX_HOMOPOLYMER_RUN:
        score += 6.0
        penalties.append(f"homopolymer run of {run}")

    # --- GC clamp at the 3' end ---
    clamp = gc_clamp(seq)
    if not clamp["has_clamp"]:
        score += 4.0
        penalties.append(f"no GC clamp (3' base {clamp['three_prime_base']})")
    if clamp["gc_in_last_5"] == 0 or clamp["gc_in_last_5"] >= 4:
        score += 3.0
        penalties.append(f"{clamp['gc_in_last_5']} GC in last 5 bases")

    # --- Self-complementarity and hairpin ---
    self_dimer = dimer_report(seq, seq)
    if self_dimer["three_prime_a"] > MAX_SELF_3PRIME:
        return None  # hard reject: 3' self-dimer of 5+
    if self_dimer["three_prime_a"] == MAX_SELF_3PRIME:
        score += 5.0
        penalties.append(
            f"3' self-complementarity of {self_dimer['three_prime_a']}"
        )
    if self_dimer["max_complementarity"] >= 6:
        score += 3.0
        penalties.append(
            f"self-complementarity run of {self_dimer['max_complementarity']}"
        )
    stem = hairpin_stem(seq)
    if stem > MAX_HAIRPIN_STEM:
        return None  # hard reject: hairpin stem of 7+ bp
    if stem >= 5:
        score += 4.0
        penalties.append(f"hairpin stem of {stem} bp")

    # --- Specificity against the template ---
    spec = specificity_report(seq, template)
    extra_exact = spec["exact_sites"] - 1
    if extra_exact < 0:
        # The primer does not occur in the template at all (e.g. a
        # hand-typed primer for `check`); not an error, just no site.
        extra_exact = 0
    if extra_exact > 0:
        score += 10.0 * extra_exact
        penalties.append(f"{spec['exact_sites']} exact binding sites")
    extra_seed = spec["seed_sites"] - 1
    if extra_seed > 0:
        score += 5.0 * extra_seed
        penalties.append(f"{spec['seed_sites']} 3'-seed binding sites")

    primer.score = round(score, 2)
    primer.penalties = penalties
    primer.checks = {
        "homopolymer_run": run,
        "gc_clamp": clamp,
        "self_complementarity": self_dimer["max_complementarity"],
        "three_prime_self": self_dimer["three_prime_a"],
        "hairpin_stem": stem,
        "exact_sites": spec["exact_sites"],
        "near_sites": spec["near_sites"],
        "seed_sites": spec["seed_sites"],
    }
    return primer


def generate_candidates(
    template: str,
    target_start: int,
    target_end: int,
    strand: str,
    config: DesignConfig | None = None,
) -> list[Primer]:
    """Generate scored candidates for one side of a target region.

    ``target_start`` / ``target_end`` are 0-based, half-open: the
    forward primer's 3' end lands near ``target_start`` and the reverse
    primer's 3' end near ``target_end``. Candidates slide their 3' end
    by up to ``config.shift`` bases either side of the boundary, and
    vary in length from ``config.min_len`` to ``config.max_len``.
    Rejected candidates are dropped; the rest are returned sorted by
    (score, Tm NN, sequence) for deterministic output.
    """
    config = config or DesignConfig()
    template = validate_sequence(template)
    if strand not in ("+", "-"):
        raise ValueError("strand must be '+' or '-'")
    if not 0 <= target_start < target_end <= len(template):
        raise ValueError(
            f"invalid target region [{target_start}, {target_end}) "
            f"for template of length {len(template)}"
        )

    candidates: list[Primer] = []
    boundary = target_start if strand == "+" else target_end
    for length in range(config.min_len, config.max_len + 1):
        for shift in range(-config.shift, config.shift + 1):
            if strand == "+":
                # 3' end at boundary + shift; 5' end length bases upstream.
                end = boundary + shift
                start = end - length
            else:
                # Reverse primer: its 3' end is at boundary + shift on
                # the template; the primer is the reverse complement of
                # template[start:start+length] with start = boundary + shift.
                start = boundary + shift
                end = start + length
            if start < 0 or end > len(template):
                continue
            segment = template[start:end]
            seq = segment if strand == "+" else reverse_complement(segment)
            primer = evaluate_primer(seq, start, strand, template, config)
            if primer is not None:
                candidates.append(primer)

    candidates.sort(key=lambda p: (p.score, p.tm_nn, p.seq))
    return candidates
