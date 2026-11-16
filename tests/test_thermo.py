"""Tests for melting-temperature and base-composition calculations.

Reference values come from two published, independently checkable
sources:

* Biopython documentation example (Bio.SeqUtils.MeltingTemp):
  ``CGTTCCAAAGATGTGGGCATGAGCTTAC`` has Tm_Wallace = 84.00 degC and
  Tm_NN = 60.32 degC with Biopython's defaults (25 nM primer, 25 nM
  template, 50 mM Na+).
* The same sequence and the universal 16S primers 27F / 1492R were
  cross-checked against Biopython 1.8x ``MeltingTemp.Tm_NN`` with
  250 nM primer, 0 nM template, 50 mM Na+ (see docs/RESEARCH.md).
"""

from __future__ import annotations

import math

import pytest

from primer_designer.thermo import (
    _nn_sums,
    gc_content,
    reverse_complement,
    tm_gc_basic,
    tm_nn,
    tm_wallace,
    validate_sequence,
)

SEQ_28 = "CGTTCCAAAGATGTGGGCATGAGCTTAC"  # Biopython docs example
PRIMER_27F = "AGAGTTTGATCATGGCTCAG"  # E. coli 16S positions 8-27
PRIMER_1492R = "GGTTACCTTGTTACGACTT"  # E. coli 16S positions 1492-1510


# --- Wallace rule: hand-computed from base counts --------------------


def test_tm_wallace_biopython_example():
    # 14 A/T + 14 G/C -> 2*14 + 4*14 = 84.00 (Biopython docs: 84.00).
    assert tm_wallace(SEQ_28) == pytest.approx(84.00)


def test_tm_wallace_27f():
    # 11 A/T + 9 G/C -> 2*11 + 4*9 = 58.00.
    assert tm_wallace(PRIMER_27F) == pytest.approx(58.00)


def test_tm_wallace_1492r():
    # 11 A/T + 8 G/C -> 2*11 + 4*8 = 54.00.
    assert tm_wallace(PRIMER_1492R) == pytest.approx(54.00)


def test_tm_wallace_all_gc_and_all_at():
    assert tm_wallace("GCGCGCGC") == pytest.approx(32.00)
    assert tm_wallace("ATATATAT") == pytest.approx(16.00)


# --- Nearest neighbour: published / cross-checked values -------------


def test_tm_nn_biopython_docs_default():
    # Biopython docs: Tm_NN(SEQ_28) = 60.32 with Biopython's defaults
    # (dnac1=25 nM, dnac2=25 nM, Na=50 mM).
    assert tm_nn(SEQ_28, na_mM=50, primer_nM=25, template_nM=25) == pytest.approx(
        60.32, abs=0.05
    )


def test_tm_nn_pcr_conditions():
    # Cross-checked against Biopython Tm_NN (250 nM primer, 0 template,
    # 50 mM Na+).
    assert tm_nn(SEQ_28) == pytest.approx(63.32, abs=0.05)
    assert tm_nn(PRIMER_27F) == pytest.approx(53.48, abs=0.05)
    assert tm_nn(PRIMER_1492R) == pytest.approx(50.80, abs=0.05)
    assert tm_nn("GCGCTCGTCGTCT") == pytest.approx(52.07, abs=0.05)
    assert tm_nn("ATGCATGCATGCATGCATGC") == pytest.approx(59.59, abs=0.05)


def test_tm_nn_owczarzy_example_1m_salt():
    # Owczarzy et al. (2008) worked example: AAGGCGAGTCAGGCTCAGTG
    # has Tm 76.3 degC at 1 M Na+ under their concentration convention;
    # Biopython's convention (250 nM primer, C_eff = primer) gives
    # 75.98 degC -- the ~0.3 degC gap is the convention difference,
    # documented in docs/RESEARCH.md. We pin the Biopython value.
    assert tm_nn("AAGGCGAGTCAGGCTCAGTG", na_mM=1000) == pytest.approx(
        75.98, abs=0.05
    )


def test_tm_nn_salt_dependence():
    # More salt stabilises the duplex: Tm rises with [Na+].
    low = tm_nn(PRIMER_27F, na_mM=10)
    mid = tm_nn(PRIMER_27F, na_mM=50)
    high = tm_nn(PRIMER_27F, na_mM=1000)
    assert low < mid < high


def test_tm_nn_concentration_dependence():
    # More primer stabilises the duplex: Tm rises with concentration.
    dilute = tm_nn(PRIMER_27F, primer_nM=25)
    standard = tm_nn(PRIMER_27F, primer_nM=250)
    concentrated = tm_nn(PRIMER_27F, primer_nM=1000)
    assert dilute < standard < concentrated


def test_tm_nn_template_lowers_effective_concentration():
    # Adding template (dnac2) lowers C_eff = primer - template/2.
    assert tm_nn(SEQ_28, primer_nM=250, template_nM=0) > tm_nn(
        SEQ_28, primer_nM=250, template_nM=250
    )


def test_tm_nn_rejects_bad_concentrations():
    with pytest.raises(ValueError):
        tm_nn(PRIMER_27F, na_mM=0)
    with pytest.raises(ValueError):
        tm_nn(PRIMER_27F, primer_nM=100, template_nM=250)


def test_nn_sums_hand_computed_with_symmetry():
    # ATGCAT is self-complementary: steps AT, TG(=CA/GT), GC, CA, AT.
    # dH = -7.2-8.5-9.8-8.5-7.2 = -41.2; initiation 2x (2.3, 4.1).
    # dS = -20.4-22.7-24.4-22.7-20.4 = -110.6, + 2x4.1, -1.4 symmetry.
    dh, ds = _nn_sums("ATGCAT")
    assert dh == pytest.approx(-36.6)
    assert ds == pytest.approx(-103.8)


def test_tm_nn_order_matters_not_just_composition():
    # Same base counts, different order -> different NN Tm, while
    # Wallace (composition-only) is identical. The two sequences have
    # the same GC count (8/16) but different dinucleotide multisets.
    a = "GCGCGCGCATATATAT"
    b = "GGGGCCCCATATATAT"
    assert tm_wallace(a) == tm_wallace(b)
    assert abs(tm_nn(a) - tm_nn(b)) > 1.0


# --- GC content and basic GC Tm ---------------------------------------


def test_gc_content():
    assert gc_content("GCGCGCGC") == pytest.approx(100.0)
    assert gc_content("ATATATAT") == pytest.approx(0.0)
    assert gc_content(PRIMER_27F) == pytest.approx(45.0)  # 9/20
    assert gc_content("ACGT") == pytest.approx(50.0)


def test_tm_gc_basic_hand_computed():
    # SEQ_28: 14 GC in 28 nt -> 64.9 + 41*(14-16.4)/28 = 61.39.
    assert tm_gc_basic(SEQ_28) == pytest.approx(61.39, abs=0.01)


# --- Helpers ------------------------------------------------------------


def test_reverse_complement():
    assert reverse_complement(PRIMER_27F) == "CTGAGCCATGATCAAACTCT"
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement(reverse_complement(PRIMER_1492R)) == PRIMER_1492R


def test_validate_sequence_accepts_lowercase():
    assert validate_sequence("acgt") == "ACGT"


def test_validate_sequence_rejects_degenerate_and_empty():
    with pytest.raises(ValueError):
        validate_sequence("ACGTN")
    with pytest.raises(ValueError):
        validate_sequence("AGTX")
    with pytest.raises(ValueError):
        validate_sequence("   ")


def test_tm_functions_reject_invalid_sequence():
    with pytest.raises(ValueError):
        tm_wallace("ACGTN")
    with pytest.raises(ValueError):
        gc_content("")
