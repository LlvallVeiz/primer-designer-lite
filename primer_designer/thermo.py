"""Melting-temperature (Tm) and base-composition calculations for DNA primers.

Two Tm estimators are implemented, both for plain A/C/G/T sequences
(5' -> 3'):

* :func:`tm_wallace` -- the Wallace rule of thumb, 2 degC per A/T and
  4 degC per G/C (Wallace et al., 1979). Fast, but it ignores base order
  and is only meant for short oligonucleotides.
* :func:`tm_nn` -- the nearest-neighbour (NN) two-state model with the
  unified DNA/DNA parameters of SantaLucia (1998), including the
  monovalent-salt correction on entropy. This is the primary estimator
  used for primer scoring.

Concentration convention (documented so results are reproducible): the
effective strand concentration is ``C_eff = primer_nM - template_nM / 2``
in nanomolar, following the convention used by Biopython's
``Bio.SeqUtils.MeltingTemp.Tm_NN``. With the defaults (250 nM primer,
no template, 50 mM Na+) the Biopython documentation example
``CGTTCCAAAGATGTGGGCATGAGCTTAC`` gives Tm_NN = 63.32 degC; with
Biopython's own defaults (25 nM primer, 25 nM template, 50 mM Na+) the
same sequence gives 60.32 degC, exactly as published in the Biopython
documentation. See ``docs/RESEARCH.md`` and ``tests/test_thermo.py``.
"""

from __future__ import annotations

import math

__all__ = [
    "gc_content",
    "tm_wallace",
    "tm_gc_basic",
    "tm_nn",
    "validate_sequence",
    "reverse_complement",
]

# Universal gas constant, cal / (mol * K).
R_GAS = 1.987

# Unified DNA/DNA nearest-neighbour parameters (delta_H in kcal/mol,
# delta_S in cal/(mol*K)) at 1 M NaCl, SantaLucia (1998) Table 2,
# as tabulated in Biopython's DNA_NN3 table. Keys are the top-strand
# dinucleotide; symmetry-equivalent dinucleotides map to the same entry:
# TT -> AA/TT, TG -> CA/GT, AC -> GT/CA, AG -> CT/GA, TC -> GA/CT,
# CC -> GG/CC.
NN_PARAMETERS = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}

# Initiation terms (SantaLucia 1998): a terminal A/T pair contributes
# (2.3, 4.1), a terminal G/C pair contributes (0.1, -2.8). Applied once
# per end of the duplex.
INIT_AT = (2.3, 4.1)
INIT_GC = (0.1, -2.8)

# Symmetry correction on entropy (cal/(mol*K)) for self-complementary
# sequences, SantaLucia (1998).
SYMMETRY_DS = -1.4

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def validate_sequence(seq: str) -> str:
    """Return ``seq`` upper-cased and stripped, raising ``ValueError``.

    Only the four standard DNA bases are accepted; degenerate IUPAC
    codes (N, R, Y, ...) are rejected because the NN parameters do not
    cover them.
    """
    cleaned = seq.strip().upper()
    if not cleaned:
        raise ValueError("sequence is empty")
    bad = sorted(set(cleaned) - set("ACGT"))
    if bad:
        raise ValueError(
            f"sequence contains non-ACGT base(s): {', '.join(bad)} "
            "(degenerate IUPAC codes are not supported)"
        )
    return cleaned


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of an A/C/G/T sequence."""
    return validate_sequence(seq).translate(_COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    """Return the GC content of ``seq`` as a percentage (0-100)."""
    seq = validate_sequence(seq)
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def tm_wallace(seq: str) -> float:
    """Wallace-rule Tm: 2 degC per A/T plus 4 degC per G/C.

    Example (hand-computed): ``CGTTCCAAAGATGTGGGCATGAGCTTAC``
    has 14 A/T and 14 G/C, so Tm = 2*14 + 4*14 = 84.00 degC.
    """
    seq = validate_sequence(seq)
    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    return 2.0 * at + 4.0 * gc


def tm_gc_basic(seq: str) -> float:
    """Basic GC-based Tm (Marmur & Doty form).

    Tm = 64.9 + 41 * (n_gc - 16.4) / n, a textbook estimate for
    oligonucleotides that ignores base order. Kept for comparison only;
    :func:`tm_nn` is the estimator used for scoring.
    """
    seq = validate_sequence(seq)
    n_gc = seq.count("G") + seq.count("C")
    return 64.9 + 41.0 * (n_gc - 16.4) / len(seq)


def _nn_sums(seq: str) -> tuple[float, float]:
    """Sum delta_H (kcal/mol) and delta_S (cal/(mol*K)) over a sequence.

    Adds the NN contribution of every dinucleotide, the initiation
    terms for both ends, and the symmetry correction when the sequence
    is self-complementary.
    """
    dh = 0.0
    ds = 0.0
    for i in range(len(seq) - 1):
        step_dh, step_ds = NN_PARAMETERS[seq[i : i + 2]]
        dh += step_dh
        ds += step_ds
    for end_base in (seq[0], seq[-1]):
        term = INIT_GC if end_base in "GC" else INIT_AT
        dh += term[0]
        ds += term[1]
    if seq == reverse_complement(seq):
        ds += SYMMETRY_DS
    return dh, ds


def tm_nn(
    seq: str,
    na_mM: float = 50.0,
    primer_nM: float = 250.0,
    template_nM: float = 0.0,
) -> float:
    """Nearest-neighbour Tm (SantaLucia 1998) in degrees Celsius.

    Parameters
    ----------
    seq:
        Primer sequence, 5' -> 3', A/C/G/T only.
    na_mM:
        Monovalent cation (Na+) concentration in millimolar. The
        parameters are for 1 M NaCl; the salt correction is applied to
        entropy only: ``dS += 0.368 * (N - 1) * ln[Na+]`` with [Na+] in
        molar (SantaLucia 1998). Default 50 mM, a typical PCR buffer.
    primer_nM, template_nM:
        Strand concentrations in nanomolar. The effective concentration
        is ``primer - template / 2`` (Biopython ``Tm_NN`` convention);
        with the defaults this is 250 nM primer in excess over template.
        Both must give a positive effective concentration.

    Returns
    -------
    float
        Tm = 1000 * dH / (dS_salt + R * ln(C_eff)) - 273.15.
    """
    seq = validate_sequence(seq)
    if na_mM <= 0:
        raise ValueError("na_mM must be positive")
    c_eff_nM = primer_nM - template_nM / 2.0
    if c_eff_nM <= 0:
        raise ValueError(
            "effective primer concentration must be positive "
            f"(got {c_eff_nM} nM from primer={primer_nM}, "
            f"template={template_nM})"
        )
    dh, ds = _nn_sums(seq)
    ds_salt = ds + 0.368 * (len(seq) - 1) * math.log(na_mM / 1000.0)
    tm_kelvin = (dh * 1000.0) / (ds_salt + R_GAS * math.log(c_eff_nM * 1e-9))
    return tm_kelvin - 273.15
