"""primer-designer-lite: PCR primer design from a DNA template (FASTA).

Public API (all functions are also reachable through the CLI)::

    from primer_designer.thermo import tm_nn, tm_wallace, gc_content
    from primer_designer.design import design_pair, design_batch, Target

    pairs = design_pair(template, start, end, "my-target")
"""

from __future__ import annotations

from .checks import dimer_report, gc_clamp, hairpin_stem
from .design import PrimerPair, Target, design_batch, design_pair
from .primer import DesignConfig, Primer, evaluate_primer, generate_candidates
from .specificity import find_sites, specificity_report
from .thermo import gc_content, reverse_complement, tm_gc_basic, tm_nn, tm_wallace

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "DesignConfig",
    "Primer",
    "PrimerPair",
    "Target",
    "design_batch",
    "design_pair",
    "dimer_report",
    "evaluate_primer",
    "find_sites",
    "gc_clamp",
    "gc_content",
    "generate_candidates",
    "hairpin_stem",
    "reverse_complement",
    "specificity_report",
    "tm_gc_basic",
    "tm_nn",
    "tm_wallace",
]
