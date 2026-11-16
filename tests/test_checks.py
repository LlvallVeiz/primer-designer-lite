"""Tests for primer quality checks (homopolymer, clamp, dimer, hairpin)."""

from __future__ import annotations

from primer_designer.checks import (
    dimer_report,
    gc_clamp,
    hairpin_stem,
    longest_complementarity,
    longest_homopolymer_run,
    self_report,
    three_prime_complementarity,
)

PRIMER_27F = "AGAGTTTGATCATGGCTCAG"
PRIMER_1492R = "GGTTACCTTGTTACGACTT"


def test_longest_homopolymer_run():
    assert longest_homopolymer_run("ACGTAAAAACGT") == 5
    assert longest_homopolymer_run("ACGTACGT") == 1
    assert longest_homopolymer_run("AAATTT") == 3
    assert longest_homopolymer_run(PRIMER_27F) == 3  # TTT


def test_gc_clamp_27f():
    clamp = gc_clamp(PRIMER_27F)
    assert clamp["three_prime_base"] == "G"
    assert clamp["has_clamp"] is True
    assert clamp["gc_in_last_5"] == 3  # 'CTCAG' has C, C, G


def test_gc_clamp_no_clamp():
    clamp = gc_clamp("AAAAATTTTT")
    assert clamp["three_prime_base"] == "T"
    assert clamp["has_clamp"] is False
    assert clamp["gc_in_last_5"] == 0


def test_longest_complementarity_full_dimer():
    # AAAAAA pairs fully with TTTTTT.
    assert longest_complementarity("AAAAAA", "TTTTTT") == 6
    # ACGTACGT is self-complementary: it pairs fully with itself.
    assert longest_complementarity("ACGTACGT", "ACGTACGT") == 8


def test_longest_complementarity_no_match():
    assert longest_complementarity("AAAAAAAA", "CCCCCCCC") == 0
    assert longest_complementarity("AAAAAAAA", "AAAAAAAA") == 0


def test_three_prime_complementarity_full():
    assert three_prime_complementarity("AAAA", "TTTT") == 4
    # ACGTACGT is fully self-complementary: aligning its 3' end with
    # its own reverse complement runs the whole length.
    assert three_prime_complementarity("ACGTACGT", "ACGTACGT") == 8


def test_three_prime_complementarity_partial():
    # GGGGACGT against TTTTACGT: only the 3' ACGT pairs (run of 4).
    assert three_prime_complementarity("GGGGACGT", "TTTTACGT") == 4
    # TTTTAAAA against AAAATTTT: the 3' AAAA pairs, then TTTT, run 4
    # before the alignment runs off the end of the other strand.
    assert three_prime_complementarity("TTTTAAAA", "AAAATTTT") == 4


def test_three_prime_complementarity_none():
    assert three_prime_complementarity("AAAAAAAA", "CCCCCCCC") == 0


def test_hairpin_stem_hand_checked():
    # GGGG + TTT loop + CCCC: 4 bp stem, 3 nt loop.
    assert hairpin_stem("GGGGTTTCCCC") == 4
    # GCGCGCGCGC: best stem is 3 bp with a 4 nt loop.
    assert hairpin_stem("GCGCGCGCGC") == 3
    # ACGT repeat (20 nt): 8 bp stem with a 4 nt loop.
    assert hairpin_stem("ACGTACGTACGTACGTACGT") == 8


def test_hairpin_stem_none():
    # No complementary pair separated by a >= 3 nt loop.
    assert hairpin_stem("AAAAAAAA") == 0
    assert hairpin_stem("ACGT") == 0


def test_dimer_report_universal_primers_are_clean():
    # The published 27F/1492R pair is designed to avoid dimers: the
    # longest cross-complementarity and both 3' runs stay short.
    report = dimer_report(PRIMER_27F, PRIMER_1492R)
    assert report["max_complementarity"] <= 5
    assert report["three_prime_a"] <= 3
    assert report["three_prime_b"] <= 3


def test_self_report_27f():
    report = self_report(PRIMER_27F)
    assert report["homopolymer_run"] == 3
    assert report["gc_clamp"]["has_clamp"] is True
    assert report["three_prime_self"] <= 3
    assert report["hairpin_stem"] <= 4
