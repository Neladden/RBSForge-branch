"""
Operon Calculator -- Module E: the sigma70 Promoter Calculator.

A from-scratch, dependency-free reimplementation of the LaFleur, Hossain &
Salis (2022) 346-parameter linear free-energy model for sigma70 promoter
transcription rate, scanning both strands of an input DNA sequence for
every candidate transcription start site.

    from operon.promoter_calculator import predict, scan_promoters

    result = predict("...promoter DNA...", organism="ecoli")
    best = result.best
    print(best.tss, best.strand, best.tx_rate, best.dg_total)

    cryptic = scan_promoters(assembled_dna, intended_tss=my_tss, tau_tx=1.0)

See ``docs/MODEL.md`` in this package for the full formula set, citations,
and calibration notes; see the top-level ``docs/OPERON_CALCULATOR.md`` for
how this module fits into the rest of the Operon Calculator.
"""

from .calculator import PromoterCalculator, predict, scan_promoters
from .energy_model import PromoterEnergyBreakdown, SigmaFactorElements, score_elements
from .types import PromoterHit, PromoterScanResult

__all__ = [
    "PromoterCalculator",
    "predict",
    "scan_promoters",
    "PromoterHit",
    "PromoterScanResult",
    "PromoterEnergyBreakdown",
    "SigmaFactorElements",
    "score_elements",
]
