"""Binding-site search: how specifically does a primer hit a template?

A primer is specific when it has exactly one binding site in the
template. Two levels are checked:

* **Exact / near sites** -- sliding-window scan of the template for the
  primer sequence (forward strand) and its reverse complement (reverse
  strand), allowing up to ``max_mismatches`` mismatches. A near site
  with few mismatches is still a mispriming risk, especially when the
  mismatches sit away from the 3' end.
* **3' seed sites** -- sites where the 3'-most 12 nt of the primer match
  exactly. Polymerase extension only needs the 3' end to be paired, so
  a perfect 3' seed is a strong off-target signal even when the 5' end
  differs (this mirrors the seed logic used by Primer3's mispriming
  checks and by BLAST-style primer specificity screens).
"""

from __future__ import annotations

from .thermo import reverse_complement, validate_sequence

__all__ = ["find_sites", "specificity_report"]

DEFAULT_SEED_LENGTH = 12


def _mismatches(a: str, b: str, limit: int | None = None) -> int:
    """Hamming distance between equal-length strings, with early exit.

    Stops counting (and returns a value above ``limit``) as soon as
    more than ``limit`` mismatches are seen, so scans of a long
    template stay fast.
    """
    count = 0
    for x, y in zip(a, b):
        if x != y:
            count += 1
            if limit is not None and count > limit:
                return count
    return count


def _seed_site_count(primer: str, template: str) -> int:
    """Count template sites where the primer's 3' seed matches exactly.

    The seed is the primer's 3'-most 12 nt. On the forward strand the
    seed occurs as-is; on the reverse strand it occurs as its reverse
    complement at the site start. Found with ``str.find`` instead of a
    full sliding-window scan. Only positions where the full primer
    fits are counted, matching :func:`find_sites` semantics.
    """
    seed_length = min(DEFAULT_SEED_LENGTH, len(primer))
    seed = primer[-seed_length:]
    rc_seed = reverse_complement(seed)
    limit = len(template) - len(primer)
    if limit < 0:
        return 0
    count = 0
    for query in (seed, rc_seed):
        start = 0
        while True:
            hit = template.find(query, start)
            if hit < 0 or hit > limit:
                break
            count += 1
            start = hit + 1
    return count


def find_sites(
    primer: str,
    template: str,
    max_mismatches: int = 0,
    seed_length: int = DEFAULT_SEED_LENGTH,
) -> list[dict]:
    """Find binding sites of ``primer`` in ``template``.

    Scans both strands: the primer sequence itself (forward strand)
    and its reverse complement (reverse strand). Returns a list of
    ``{"start": int, "strand": "+" | "-", "mismatches": int,
    "seed_match": bool}`` sorted by (mismatches, start). ``start`` is
    the 0-based position of the site on the template. ``seed_match``
    is True when the primer's 3'-most ``seed_length`` bases match the
    site exactly.
    """
    primer = validate_sequence(primer)
    template = validate_sequence(template)
    if len(primer) > len(template):
        return []
    sites: list[dict] = []
    seed = primer[-seed_length:] if len(primer) >= seed_length else primer
    for strand, query in (("+", primer), ("-", reverse_complement(primer))):
        for start in range(len(template) - len(query) + 1):
            window = template[start : start + len(query)]
            mm = _mismatches(query, window, max_mismatches)
            if mm > max_mismatches:
                continue
            # The 3' end of the primer is the last base of `query` on
            # the forward strand, and the first base of `query` on the
            # reverse strand (query is the reverse complement there,
            # so its 5' end is the primer's 3' end read backwards).
            if strand == "+":
                site_seed = window[-len(seed):]
            else:
                site_seed = reverse_complement(window[: len(seed)])
            sites.append(
                {
                    "start": start,
                    "strand": strand,
                    "mismatches": mm,
                    "seed_match": site_seed == seed,
                }
            )
    sites.sort(key=lambda s: (s["mismatches"], s["start"]))
    return sites


def specificity_report(primer: str, template: str) -> dict:
    """Summarise the specificity of ``primer`` against ``template``.

    Counts exact sites, near sites (<= 2 mismatches), and sites where
    the 3' 12-mer seed matches exactly (including sites found only by
    the seed, i.e. with more than 2 mismatches elsewhere). A primer
    is specific when ``exact_sites == 1`` and ``seed_sites == 1``.
    """
    exact = find_sites(primer, template, max_mismatches=0)
    near = find_sites(primer, template, max_mismatches=2)
    seed_count = _seed_site_count(primer, template)
    return {
        "exact_sites": len(exact),
        "near_sites": len(near),
        "seed_sites": seed_count,
        "exact": exact,
        "near": near,
    }
