"""Tests for primer evaluation, scoring, and candidate generation."""

from __future__ import annotations

import pytest

from primer_designer.primer import DesignConfig, evaluate_primer, generate_candidates
from primer_designer.thermo import reverse_complement

PRIMER_27F = "AGAGTTTGATCATGGCTCAG"
PRIMER_1492R = "GGTTACCTTGTTACGACTT"


def test_evaluate_27f_is_good(template_16s):
    primer = evaluate_primer(PRIMER_27F, 7, "+", template_16s)
    assert primer is not None
    assert primer.seq == PRIMER_27F
    assert primer.tm_nn == pytest.approx(53.48, abs=0.05)
    assert primer.tm_wallace == pytest.approx(58.00)
    assert primer.gc_percent == pytest.approx(45.0)
    # Tm 53.48 is below the accepted 55-65 window, and the primer has
    # a 6 bp self-complementary run (the GATC motif), so penalties are
    # expected -- but no hard rejection and no specificity penalty.
    assert primer.checks["exact_sites"] == 1
    assert primer.checks["seed_sites"] == 1
    assert primer.score == pytest.approx(16.03, abs=0.05)
    assert any("self-complementarity run of 6" in p for p in primer.penalties)


def test_evaluate_1492r_is_good(template_16s):
    primer = evaluate_primer(PRIMER_1492R, 1491, "-", template_16s)
    assert primer is not None
    assert primer.tm_nn == pytest.approx(50.80, abs=0.05)
    assert primer.checks["exact_sites"] == 1


def test_evaluate_rejects_homopolymer_run_of_5(template_16s):
    assert evaluate_primer("ACGTAAAAAACGTACGTACG", 0, "+", template_16s) is None


def test_evaluate_rejects_hairpin_stem_of_7(template_16s):
    # GCGCGCG + TTT loop + CGCGCGC folds into a 7 bp stem.
    assert evaluate_primer("GCGCGCGTTTCGCGCGC", 0, "+", template_16s) is None


def test_evaluate_rejects_strong_3prime_self_dimer(template_16s):
    # ACGT repeat is fully self-complementary with a long 3' run.
    assert evaluate_primer("ACGTACGTACGTACGTACGT", 0, "+", template_16s) is None


def test_evaluate_penalises_missing_gc_clamp(template_16s):
    # Ends with A (no 3' GC clamp); GC 50%, Tm in range.
    primer = evaluate_primer("GCTAGCTAGCTAGGATCGTA", 0, "+", template_16s)
    assert primer is not None
    assert any("no GC clamp" in p for p in primer.penalties)


def test_evaluate_penalises_low_gc_and_low_tm(template_16s):
    primer = evaluate_primer("TATAAATATTAATATAAATA", 0, "+", template_16s)
    assert primer is not None
    assert primer.gc_percent == pytest.approx(0.0)
    assert primer.tm_nn == pytest.approx(31.6, abs=0.5)
    assert any("GC 0.0% below 40%" in p for p in primer.penalties)
    assert any("below 55" in p for p in primer.penalties)
    assert primer.score > 10.0


def test_evaluate_penalises_extra_binding_site():
    template = "TT" + PRIMER_27F + "AAAA" + PRIMER_27F + "TT"
    primer = evaluate_primer(PRIMER_27F, 2, "+", template)
    assert primer is not None
    assert primer.checks["exact_sites"] == 2
    assert any("2 exact binding sites" in p for p in primer.penalties)


def test_evaluate_primer_outside_template_is_allowed(template_16s):
    # A hand-typed primer with no site in the template is evaluated
    # (no specificity penalty), not rejected.
    primer = evaluate_primer("TTGACCTGACCTGACCTGAC", 0, "+", template_16s)
    assert primer is not None
    assert primer.checks["exact_sites"] == 0


def test_generate_candidates_forward(template_16s):
    # Full-length 16S target: forward 3' end near 0-based 26.
    candidates = generate_candidates(template_16s, 26, 1492, "+")
    assert len(candidates) > 50
    assert all(c.strand == "+" for c in candidates)
    assert all(18 <= c.length <= 25 for c in candidates)
    # Sorted best-first by score.
    scores = [c.score for c in candidates]
    assert scores == sorted(scores)
    # Every candidate's 3' end lands within shift=10 of the boundary.
    assert all(16 <= c.end <= 36 for c in candidates)


def test_generate_candidates_reverse(template_16s):
    candidates = generate_candidates(template_16s, 26, 1492, "-")
    assert len(candidates) > 50
    assert all(c.strand == "-" for c in candidates)
    # Reverse candidates are reverse complements of template segments
    # starting within shift=10 of the boundary (0-based 1492).
    assert all(1482 <= c.start <= 1502 for c in candidates)
    for c in candidates:
        segment = template_16s[c.start : c.start + c.length]
        assert c.seq == reverse_complement(segment)


def test_generate_candidates_rejects_bad_region(template_16s):
    with pytest.raises(ValueError):
        generate_candidates(template_16s, 100, 50, "+")
    with pytest.raises(ValueError):
        generate_candidates(template_16s, 0, len(template_16s) + 10, "+")
    with pytest.raises(ValueError):
        generate_candidates(template_16s, 26, 1492, "x")


def test_design_config_defaults():
    config = DesignConfig()
    assert (config.min_len, config.max_len) == (18, 25)
    assert (config.tm_min, config.tm_max) == (55.0, 65.0)
    assert config.na_mM == 50.0
    assert config.primer_nM == 250.0
