"""
Operon Calculator -- Module F: highly translated internal start codons
(HTISC) (spec section 11).

This module is RBSForge run as a scanner, plus a threshold, and needs no
new physics -- of everything in the spec's "what already works without
new RBS physics" list (section 4.6), this is the cheapest to finish:

    def scan_htisc(mrna, annotated_starts, calc, tau_htisc=1000.0):
        hits = []
        out = calc.predict(mrna)
        for r in out.results:
            if r.leaderless or r.index in annotated_starts:
                continue
            if r.v1_style_rate >= tau_htisc:
                hits.append(r)
        return hits

Status: planned, not yet implemented as a package function (the one-liner
above is usable directly against ``rbsforge.RBSCalculator.predict`` today
and is deliberately reproduced here rather than deferred, since it needs
nothing this repository doesn't already have).
"""

DEFAULT_TAU_HTISC = 1000.0  # au, v1_style_rate scale -- see rbsforge docs on TIR-scale freezing


def scan_htisc(mrna: str, annotated_starts, calc, tau_htisc: float = DEFAULT_TAU_HTISC) -> list:
    """Every non-annotated, non-leaderless start codon whose
    ``v1_style_rate`` clears ``tau_htisc`` -- the "Translated Open Reading
    Frames" plot in the web tool."""
    raise NotImplementedError("Module F (HTISC) is a prepared subsection, not yet implemented")
