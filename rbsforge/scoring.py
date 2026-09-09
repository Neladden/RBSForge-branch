"""Transparent scoring primitives shared by the RBS calculator.

Every score is an additive decomposition of named, physically-
interpretable free-energy terms -- never a single opaque number from a
fitted black-box model. `ScoreBreakdown` is the container the RBS
calculator returns for each start codon: it holds the six terms of

    delta_G_total = delta_G_standby + delta_G_mRNA:rRNA + delta_G_spacing
                    + delta_G_start + delta_G_stacking - delta_G_mRNA

so callers can always see exactly why a sequence scored the way it did
(model summary section 3, "Master equations").
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict

from .constants import BETA_DEFAULT, RT_EFF_DEFAULT, V1_PREFACTOR_K


@dataclass
class ScoreBreakdown:
    """An additive decomposition of delta_G_total.

    Attributes:
        terms: ordered mapping of term name -> contribution (kcal/mol).
            Positive terms weaken/destabilize initiation; negative terms
            strengthen it -- standard free-energy sign convention.
        notes: free-text explanation per term, keyed the same as `terms`.
    """

    terms: Dict[str, float] = field(default_factory=dict)
    notes: Dict[str, str] = field(default_factory=dict)

    def add(self, name: str, value: float, note: str = "") -> None:
        self.terms[name] = value
        if note:
            self.notes[name] = note

    @property
    def total(self) -> float:
        return sum(self.terms.values())

    def to_dict(self) -> dict:
        return {"total": self.total, "terms": dict(self.terms), "notes": dict(self.notes)}

    def __str__(self) -> str:
        lines = [f"delta_G_total = {self.total:+.2f} kcal/mol"]
        for name, value in self.terms.items():
            note = f"  # {self.notes[name]}" if name in self.notes else ""
            lines.append(f"  {name:<16} {value:+7.2f}{note}")
        return "\n".join(lines)


def proportional_rate(delta_g_total: float, beta: float = BETA_DEFAULT) -> float:
    """r = exp(-beta * delta_G_total): the v2.x-style proportional TIR.

    beta is an *empirical* free-energy-to-rate conversion (see
    constants.py and docs/MODEL.md), not 1/RT -- only ratios of this value
    (within one HostPack/temperature) are meaningful, not absolute scale.
    """
    return math.exp(-beta * delta_g_total)


def v1_style_rate(
    delta_g_total: float,
    rt_eff: float = RT_EFF_DEFAULT,
    k_prefactor: float = V1_PREFACTOR_K,
) -> float:
    """r = K * exp(-delta_G_total / RT_eff): the v1.0-style prefactored TIR.

    Provided for comparison with published v1.0 worked examples only. Do
    not mix this with `proportional_rate` values -- they use different
    zeros of energy (model summary section 5).
    """
    return k_prefactor * math.exp(-delta_g_total / rt_eff)
