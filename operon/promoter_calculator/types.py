"""Result types for the sigma70 promoter scan."""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from .energy_model import PromoterEnergyBreakdown

Span = Tuple[int, int]  # [start, end) on the scanned DNA, forward-strand coordinates


@dataclass(frozen=True)
class PromoterHit:
    """One scored candidate sigma70 promoter state: the minimum-dG_total
    element combination found for a given TSS position and strand."""

    tss: int
    strand: int  # +1 or -1
    tx_rate: float
    dg_total: float
    breakdown: PromoterEnergyBreakdown
    promoter_span: Span
    up_span: Span
    hex35_span: Span
    spacer_span: Span
    hex10_span: Span
    disc_span: Span
    itr_span: Span
    up: str
    hex35: str
    spacer: str
    hex10: str
    disc: str
    itr: str

    def is_intended(self, intended_tss: Optional[int], intended_strand: int = 1) -> bool:
        return intended_tss is not None and self.tss == intended_tss and self.strand == intended_strand


@dataclass
class PromoterScanResult:
    """Every TSS position's best (minimum dG_total) hit, both strands, from
    one call to :func:`operon.promoter_calculator.calculator.scan`."""

    sequence: str
    organism: str
    k: float
    beta: float
    forward: Dict[int, PromoterHit]
    reverse: Dict[int, PromoterHit]

    def all_hits(self) -> List[PromoterHit]:
        return list(self.forward.values()) + list(self.reverse.values())

    def ranked(self) -> List[PromoterHit]:
        """All hits, strongest (highest Tx_rate) first."""
        return sorted(self.all_hits(), key=lambda h: h.tx_rate, reverse=True)

    @property
    def best(self) -> Optional[PromoterHit]:
        hits = self.ranked()
        return hits[0] if hits else None

    def internal(self, intended_tss: Optional[int] = None, intended_strand: int = 1,
                 tau_tx: float = 0.0) -> List[PromoterHit]:
        """Hits that are not the intended promoter and clear ``tau_tx``,
        ranked strongest first -- the "cryptic promoter" list (spec
        Module E, section 10.3)."""
        return [
            h for h in self.ranked()
            if not h.is_intended(intended_tss, intended_strand) and h.tx_rate > tau_tx
        ]
