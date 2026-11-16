"""Tests for binding-site search / specificity against a template."""

from __future__ import annotations

from primer_designer.specificity import find_sites, specificity_report

PRIMER_27F = "AGAGTTTGATCATGGCTCAG"  # E. coli 16S positions 8-27
PRIMER_1492R = "GGTTACCTTGTTACGACTT"  # E. coli 16S positions 1492-1510


def test_find_sites_27f_single_forward_site(template_16s):
    sites = find_sites(PRIMER_27F, template_16s)
    assert len(sites) == 1
    assert sites[0]["start"] == 7  # 0-based; 1-based position 8
    assert sites[0]["strand"] == "+"
    assert sites[0]["mismatches"] == 0
    assert sites[0]["seed_match"] is True


def test_find_sites_1492r_single_reverse_site(template_16s):
    sites = find_sites(PRIMER_1492R, template_16s)
    assert len(sites) == 1
    assert sites[0]["start"] == 1491  # 0-based; 1-based position 1492
    assert sites[0]["strand"] == "-"
    assert sites[0]["mismatches"] == 0
    assert sites[0]["seed_match"] is True


def test_find_sites_absent_primer(template_16s):
    # A primer with no binding site at all (poly-G does not occur).
    assert find_sites("GGGGGGGGGGGGGGGGGGGG", template_16s) == []


def test_find_sites_duplicated_primer():
    template = "TT" + PRIMER_27F + "AAAA" + PRIMER_27F + "TT"
    sites = find_sites(PRIMER_27F, template)
    assert len(sites) == 2
    assert [s["start"] for s in sites] == [2, 26]


def test_find_sites_mismatch_tolerance():
    template = "TT" + PRIMER_27F + "TT"
    # One mismatch outside the 3' 12-mer seed (position 5; the seed
    # covers positions 8-19), so the seed still matches exactly.
    mutated = PRIMER_27F[:5] + ("A" if PRIMER_27F[5] != "A" else "C") + PRIMER_27F[6:]
    assert find_sites(mutated, template, max_mismatches=0) == []
    sites = find_sites(mutated, template, max_mismatches=1)
    assert len(sites) == 1
    assert sites[0]["mismatches"] == 1
    # The mismatch is at position 10 (not in the 3' 12-mer seed),
    # so the seed still matches exactly.
    assert sites[0]["seed_match"] is True


def test_find_sites_seed_mismatch_breaks_seed():
    template = "TT" + PRIMER_27F + "TT"
    # One mismatch inside the 3' 12-mer seed (last base).
    mutated = PRIMER_27F[:-1] + ("A" if PRIMER_27F[-1] != "A" else "C")
    sites = find_sites(mutated, template, max_mismatches=1)
    assert len(sites) == 1
    assert sites[0]["seed_match"] is False


def test_specificity_report_universal_primers(template_16s):
    for primer in (PRIMER_27F, PRIMER_1492R):
        report = specificity_report(primer, template_16s)
        assert report["exact_sites"] == 1
        assert report["near_sites"] == 1
        assert report["seed_sites"] == 1


def test_find_sites_primer_longer_than_template():
    assert find_sites("ACGTACGTACGTACGTACGTACGT", "ACGT") == []
