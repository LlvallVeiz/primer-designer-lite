"""Quality checks for individual primers and primer pairs.

All checks are simple, deterministic sequence scans -- no alignment
library is needed at primer length (18-25 nt). Each function returns
plain values so the scoring code in :mod:`primer_designer.primer` can
turn them into penalties, and so the ``check`` CLI subcommand can
report them directly.

Conventions
-----------
* A primer is always written 5' -> 3'.
* "3'-anchored" complementarity means the match involves the 3' end
  of the primer: the last base of the primer pairs with a base of the
  other strand. A 3' match matters most because DNA polymerase extends
  from a paired 3' end, so a 3' dimer can prime unwanted extension.
* Hairpin stems are searched with a minimum loop of 3 nt, the smallest
  loop that can close without strain in this simplified model.
"""

from __future__ import annotations

from .thermo import reverse_complement, validate_sequence

__all__ = [
    "longest_homopolymer_run",
    "gc_clamp",
    "longest_complementarity",
    "three_prime_complementarity",
    "hairpin_stem",
    "dimer_report",
    "self_report",
]

_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def longest_homopolymer_run(seq: str) -> int:
    """Length of the longest run of one repeated base in ``seq``.

    A run of 4+ identical bases in a row (e.g. ``AAAA``) is a classic
    primer-design flag: it risks polymerase slippage and synthesis
    errors, and it lowers sequence complexity.
    """
    seq = validate_sequence(seq)
    best = 1
    run = 1
    for prev, cur in zip(seq, seq[1:]):
        run = run + 1 if cur == prev else 1
        best = max(best, run)
    return best


def gc_clamp(seq: str) -> dict:
    """GC-clamp facts for the 3' end of ``seq``.

    Returns the 3' base, whether it is G/C ("clamp present"), and the
    number of G/C bases among the last five. Good practice (Primer3
    manual; Untergasser et al., 2012) is a G or C at the 3' end for
    tighter binding where extension starts, with 1-3 G/C in the last
    five bases -- more than that over-stabilises the 3' end and invites
    mispriming.
    """
    seq = validate_sequence(seq)
    tail = seq[-5:]
    return {
        "three_prime_base": seq[-1],
        "has_clamp": seq[-1] in "GC",
        "gc_in_last_5": tail.count("G") + tail.count("C"),
    }


def _max_run_at_alignment(a: str, b_rc: str, shift: int) -> int:
    """Longest run of paired bases for one alignment shift.

    ``a`` is aligned against ``b_rc`` (the reverse complement of the
    second strand) with ``shift`` offset: position ``i`` of ``a`` pairs
    with position ``i - shift`` of ``b_rc``. Two strands pair exactly
    where the first strand equals the reverse complement of the second,
    so a pair is counted when ``a[i] == b_rc[j]``.
    """
    best = 0
    run = 0
    for i, base in enumerate(a):
        j = i - shift
        if 0 <= j < len(b_rc) and base == b_rc[j]:
            run += 1
            best = max(best, run)
        else:
            run = 0
    return best


def longest_complementarity(a: str, b: str) -> int:
    """Longest run of complementary bases between primers ``a`` and ``b``.

    The maximum over all ungapped alignments of ``a`` against the
    reverse complement of ``b``. For ``a == b`` this is the primer's
    self-complementarity (self-dimer potential).
    """
    a = validate_sequence(a)
    b_rc = reverse_complement(b)
    best = 0
    for shift in range(-len(b_rc) + 1, len(a)):
        best = max(best, _max_run_at_alignment(a, b_rc, shift))
    return best


def three_prime_complementarity(a: str, b: str) -> int:
    """Longest 3'-anchored complementary run of ``a`` against ``b``.

    Only alignments in which the 3' end (last base) of ``a`` is paired
    are counted, and the run must include that last base. This is the
    dimer signal that matters for PCR: a paired 3' end can be extended
    by the polymerase. Returns 0 when the 3' base of ``a`` pairs with
    nothing in ``b``.
    """
    a = validate_sequence(a)
    b_rc = reverse_complement(b)
    best = 0
    last = len(a) - 1
    for j in range(len(b_rc)):
        # Align a[last] with b_rc[j]; count backwards while the bases
        # pair (a[i] == b_rc[k]).
        run = 0
        i, k = last, j
        while i >= 0 and k >= 0 and a[i] == b_rc[k]:
            run += 1
            i -= 1
            k -= 1
        best = max(best, run)
    return best


def hairpin_stem(seq: str, min_loop: int = 3) -> int:
    """Longest hairpin stem (in base pairs) that ``seq`` can form.

    A stem of length ``k`` pairs the 5' arm ``seq[i:i+k]`` antiparallel
    with the 3' arm ``seq[e-k+1:e+1]`` -- base ``i+m`` pairs with base
    ``e-m`` -- leaving a loop of at least ``min_loop`` unpaired bases
    between the arms. The scan is exhaustive over all (start, end)
    pairs and extends each stem while the bases stay complementary;
    trivial at primer length. Returns 0 when no stem of length >= 1
    with a sufficient loop exists.

    Examples (hand-checked): ``GGGGTTTCCCC`` has a 4 bp stem with a
    3 nt loop; ``GCGCGCGCGC`` has a 3 bp stem with a 4 nt loop.
    """
    seq = validate_sequence(seq)
    n = len(seq)
    best = 0
    for i in range(n):
        for e in range(i + 1, n):
            # Pair (i+m) with (e-m); extend while complementary and
            # the arms have not crossed.
            k = 0
            while i + k < e - k and seq[i + k] == _COMPLEMENT[seq[e - k]]:
                k += 1
                # Unpaired bases between the arms after k pairs.
                loop = e - k - (i + k) + 1
                if loop >= min_loop:
                    best = max(best, k)
    return best


def dimer_report(a: str, b: str) -> dict:
    """Dimer report between primers ``a`` and ``b`` (both 5' -> 3').

    Reports the longest complementary run in either orientation, and
    the longest 3'-anchored run for each primer's 3' end. For a
    self-dimer report call with ``a == b``.
    """
    return {
        "max_complementarity": longest_complementarity(a, b),
        "three_prime_a": three_prime_complementarity(a, b),
        "three_prime_b": three_prime_complementarity(b, a),
    }


def self_report(seq: str) -> dict:
    """Self-complementarity and hairpin report for one primer."""
    seq = validate_sequence(seq)
    return {
        "homopolymer_run": longest_homopolymer_run(seq),
        "gc_clamp": gc_clamp(seq),
        "self_complementarity": longest_complementarity(seq, seq),
        "three_prime_self": three_prime_complementarity(seq, seq),
        "hairpin_stem": hairpin_stem(seq),
    }
