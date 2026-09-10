"""
Operon Calculator -- Module F: highly translated internal start codons
(HTISC) (spec section 11). **Status: implemented.**

RBSForge run as a scanner, plus a threshold -- no new physics. This is
the "Translated Open Reading Frames" plot in the web tool: every
candidate start codon's TIR, with annotated CDSs excluded, thresholded.
"""

from typing import Iterable, List

DEFAULT_TAU_HTISC = 1000.0  # au, v1_style_rate scale -- see rbsforge/CLAUDE.md on TIR-scale freezing


def scan_htisc(mrna: str, annotated_starts: Iterable[int], calc, tau_htisc: float = DEFAULT_TAU_HTISC) -> list:
    """Every non-annotated, non-leaderless start codon whose
    ``v1_style_rate`` clears ``tau_htisc`` -- a candidate highly
    translated internal start codon, i.e. a likely truncated-protein
    failure mode.

    ``calc`` is an ``rbsforge.RBSCalculator`` (or anything exposing the
    same ``.predict(mrna) -> PredictResult`` interface). Returns the
    underlying ``StartCodonResult`` objects directly (position via
    ``.index``, rate via ``.v1_style_rate``), ranked strongest first.
    """
    annotated = set(annotated_starts)
    out = calc.predict(mrna)
    hits: List = []
    for r in out.results:
        if r.leaderless:
            continue
        if r.index in annotated:
            continue
        if r.v1_style_rate >= tau_htisc:
            hits.append(r)
    hits.sort(key=lambda r: r.v1_style_rate, reverse=True)
    return hits
