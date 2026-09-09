"""Thermodynamics-based ribosome binding site (RBS) strength prediction.

Reimplements the free-energy decomposition of the Salis Lab RBS Calculator
(Salis, Mirsky & Voigt, Nat Biotechnol 2009, and successors -- see
docs/MODEL.md for full version history and citations) in a from-scratch,
dependency-free, HostPack-configurable form:

    delta_G_total = delta_G_standby + delta_G_mRNA:rRNA + delta_G_spacing
                    + delta_G_start + delta_G_stacking - delta_G_mRNA

    r = exp(-beta * delta_G_total)

A lower (more negative) delta_G_total corresponds to a stronger predicted
ribosome binding site / higher translation initiation rate (TIR).
"""
from .calculator import RBSCalculator, StartCodonResult, PredictResult

__all__ = ["RBSCalculator", "StartCodonResult", "PredictResult"]
