"""Pair design: pick the best forward/reverse primer pair for a target.

The designer works on a template sequence and a target region given as
0-based half-open coordinates ``[start, end)``. The forward primer's
3' end lands near ``start`` and the reverse primer's 3' end near
``end`` (see :func:`primer_designer.primer.generate_candidates`), so
the amplicon spans from the forward primer's 5' end to the reverse
primer's 5' end on the template:

    product = template[forward.start : reverse.end]
    product_length = reverse.end - forward.start

Pair scoring adds, on top of the two primers' individual scores:

* 3.0 per degC of forward/reverse NN Tm difference (a matched pair
  anneals evenly; the pair is rejected when |dTm| > ``max_delta_tm``);
* a cross-dimer penalty: 3'-anchored complementarity of 4 between the
  two primers costs 8.0, of 5+ rejects the pair (the 3' end of one
  primer would prime extension on the other);
* a product-length preference: products of 150-500 bp are ideal
  (0.01 per bp outside that window, capped at 5.0).

The recommended annealing temperature is the lower primer Tm minus
5 degC, the standard starting point for a gradient optimisation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .checks import three_prime_complementarity
from .primer import DesignConfig, Primer, generate_candidates
from .thermo import validate_sequence

__all__ = ["PrimerPair", "design_pair", "design_batch", "Target"]

MAX_CROSS_3PRIME = 4  # cross-dimer 3' run of 5+ rejects the pair


@dataclass
class Target:
    """One design target: a named region on a template (1-based inclusive)."""

    name: str
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive

    def to_zero_based(self) -> tuple[int, int]:
        """Return (start, end) as 0-based half-open coordinates."""
        return self.start - 1, self.end


@dataclass
class PrimerPair:
    """A designed forward/reverse pair plus its amplicon."""

    target: str
    forward: Primer
    reverse: Primer
    product_length: int
    product_seq: str
    delta_tm: float
    anneal_temp: float
    score: float
    penalties: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """JSON-friendly representation (rounded floats)."""
        return {
            "target": self.target,
            "forward": self.forward.to_dict(),
            "reverse": self.reverse.to_dict(),
            "product_length": self.product_length,
            "product_seq": self.product_seq,
            "delta_tm_c": round(self.delta_tm, 2),
            "anneal_temp_c": round(self.anneal_temp, 2),
            "score": round(self.score, 2),
            "penalties": list(self.penalties),
        }


def _product_length_penalty(length: int) -> float:
    """Mild preference for 150-500 bp amplicons."""
    if 150 <= length <= 500:
        return 0.0
    distance = 150 - length if length < 150 else length - 500
    return min(5.0, 0.01 * distance)


def design_pair(
    template: str,
    start: int,
    end: int,
    target_name: str = "target",
    config: DesignConfig | None = None,
) -> list[PrimerPair]:
    """Design primer pairs for target region ``[start, end)`` (0-based).

    Returns up to ``config.top_pairs`` pairs sorted by score (best
    first). Raises ``ValueError`` when no pair satisfies the constraints
    (for example, the product-length window or the Tm-difference limit).
    """
    config = config or DesignConfig()
    template = validate_sequence(template)
    forwards = generate_candidates(template, start, end, "+", config)
    reverses = generate_candidates(template, start, end, "-", config)
    if not forwards:
        raise ValueError("no forward primer candidate passed the filters")
    if not reverses:
        raise ValueError("no reverse primer candidate passed the filters")

    pairs: list[PrimerPair] = []
    for fwd in forwards[: config.top_candidates]:
        for rev in reverses[: config.top_candidates]:
            product_length = rev.end - fwd.start
            if not config.min_product <= product_length <= config.max_product:
                continue
            delta_tm = abs(fwd.tm_nn - rev.tm_nn)
            if delta_tm > config.max_delta_tm:
                continue
            cross_3prime = max(
                three_prime_complementarity(fwd.seq, rev.seq),
                three_prime_complementarity(rev.seq, fwd.seq),
            )
            if cross_3prime > MAX_CROSS_3PRIME:
                continue  # pair rejected: extendable cross-dimer
            penalties: list[str] = []
            score = fwd.score + rev.score + 3.0 * delta_tm
            if delta_tm > 2.0:
                penalties.append(f"delta Tm {delta_tm:.2f}C > 2C")
            if cross_3prime == MAX_CROSS_3PRIME:
                score += 8.0
                penalties.append(
                    f"cross-dimer 3' run of {cross_3prime} bp"
                )
            length_penalty = _product_length_penalty(product_length)
            score += length_penalty
            if length_penalty:
                penalties.append(
                    f"product {product_length} bp outside ideal 150-500 bp"
                )
            pairs.append(
                PrimerPair(
                    target=target_name,
                    forward=fwd,
                    reverse=rev,
                    product_length=product_length,
                    product_seq=template[fwd.start : rev.end],
                    delta_tm=round(delta_tm, 2),
                    anneal_temp=round(
                        min(fwd.tm_nn, rev.tm_nn) - 5.0, 2
                    ),
                    score=round(score, 2),
                    penalties=penalties,
                )
            )

    if not pairs:
        raise ValueError(
            "no primer pair satisfies the constraints "
            f"(product {config.min_product}-{config.max_product} bp, "
            f"delta Tm <= {config.max_delta_tm}C, cross-dimer 3' <= 4)"
        )
    pairs.sort(key=lambda p: (p.score, p.product_length, p.forward.seq))
    return pairs[: config.top_pairs]


def design_batch(
    template: str,
    targets: list[Target],
    config: DesignConfig | None = None,
) -> dict[str, list[PrimerPair] | str]:
    """Design pairs for several targets; failures are reported per target.

    Returns ``{target_name: [pairs...]}``; a target that cannot be
    designed maps to its error message string instead of a list, so one
    bad target never aborts the whole batch.
    """
    config = config or DesignConfig()
    template = validate_sequence(template)
    results: dict[str, list[PrimerPair] | str] = {}
    for target in targets:
        start, end = target.to_zero_based()
        try:
            results[target.name] = design_pair(
                template, start, end, target.name, config
            )
        except ValueError as exc:
            results[target.name] = str(exc)
    return results
