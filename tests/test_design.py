"""Tests for pair design and batch design.

Coordinates are 0-based half-open here: the full-length 16S target is
[26, 1492) (the 27F 3' end at 1-based 27, the 1492R 3' end at 1-based
1492), and the V4 target is [535, 787) (515F 3' end at 1-based 536,
806R 3' end at 1-based 787).
"""

from __future__ import annotations

import pytest

from primer_designer.design import Target, design_batch, design_pair
from primer_designer.primer import DesignConfig, generate_candidates
from primer_designer.thermo import reverse_complement

PRIMER_27F = "AGAGTTTGATCATGGCTCAG"
PRIMER_1492R = "GGTTACCTTGTTACGACTT"


def test_target_to_zero_based():
    assert Target("full", 27, 1492).to_zero_based() == (26, 1492)
    assert Target("v4", 536, 787).to_zero_based() == (535, 787)


def test_design_full_length_16s(template_16s):
    pairs = design_pair(template_16s, 26, 1492, "16s_full")
    assert len(pairs) == 5
    best = pairs[0]
    # Amplicon: 1466 bp between the primer 3' ends plus both primers
    # (forward 5' end in [0, 18], reverse 5' end in [1500, 1527]).
    assert 1470 <= best.product_length <= 1530
    assert best.delta_tm <= 3.0
    assert best.anneal_temp == pytest.approx(
        min(best.forward.tm_nn, best.reverse.tm_nn) - 5.0, abs=0.01
    )
    # The product really is the template slice between the primers.
    assert best.product_seq == template_16s[
        best.forward.start : best.reverse.end
    ]
    assert best.product_seq.startswith(best.forward.seq)
    assert best.product_seq.endswith(reverse_complement(best.reverse.seq))
    # Pairs are sorted best-first.
    assert [p.score for p in pairs] == sorted(p.score for p in pairs)


def test_design_recovers_universal_sites(template_16s):
    # The published universal primers must be inside the candidate
    # pool: 27F (template[7:27]) among the forward candidates and
    # 1492R (reverse complement of template[1491:1510]) among the
    # reverse candidates. (The top-ranked pairs may use other primers
    # with Tm closer to the 58-62 degC ideal window.)
    forwards = generate_candidates(template_16s, 26, 1492, "+")
    reverses = generate_candidates(template_16s, 26, 1492, "-")
    assert any(c.seq == PRIMER_27F and c.start == 7 for c in forwards)
    assert any(c.seq == PRIMER_1492R and c.start == 1491 for c in reverses)
    # And the designed pairs still land at the target boundaries.
    pairs = design_pair(template_16s, 26, 1492, "16s_full")
    assert all(16 <= p.forward.end <= 36 for p in pairs)
    assert all(1482 <= p.reverse.start <= 1502 for p in pairs)


def test_design_v4_region(template_16s):
    pairs = design_pair(template_16s, 535, 787, "16s_v4")
    assert len(pairs) >= 1
    best = pairs[0]
    # V4 amplicon: 252 bp between the primer 3' ends plus both primers.
    assert 260 <= best.product_length <= 330
    assert best.delta_tm <= 3.0
    assert best.product_seq.startswith(best.forward.seq)
    assert best.product_seq.endswith(reverse_complement(best.reverse.seq))


def test_design_synthetic_template(template_synthetic):
    pairs = design_pair(template_synthetic, 200, 400, "synthetic")
    assert len(pairs) >= 1
    best = pairs[0]
    assert 210 <= best.product_length <= 275
    assert best.delta_tm <= 3.0


def test_design_rejects_impossible_product_window(template_16s):
    config = DesignConfig(min_product=5000, max_product=6000)
    with pytest.raises(ValueError):
        design_pair(template_16s, 26, 1492, "16s_full", config)


def test_design_rejects_tiny_delta_tm(template_16s):
    config = DesignConfig(max_delta_tm=0.0)
    with pytest.raises(ValueError):
        design_pair(template_16s, 26, 1492, "16s_full", config)


def test_design_batch_all_ok(template_16s):
    targets = [
        Target("16s_full", 27, 1492),
        Target("16s_v4", 536, 787),
    ]
    results = design_batch(template_16s, targets)
    assert set(results) == {"16s_full", "16s_v4"}
    assert all(isinstance(v, list) and len(v) >= 1 for v in results.values())


def test_design_batch_one_failure_does_not_abort(template_16s):
    targets = [
        Target("16s_v4", 536, 787),
        Target("impossible", 27, 28),  # 2 bp target: no valid product
    ]
    results = design_batch(template_16s, targets)
    assert isinstance(results["16s_v4"], list)
    assert isinstance(results["impossible"], str)
    assert "no primer pair" in results["impossible"]
