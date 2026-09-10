"""
sigma70 promoter scanner (Operon Calculator spec, Module E / section 10).

For every candidate transcription start site (TSS) in a window, and every
allowed discriminator length (6-10 nt) x spacer length (15-20 nt)
combination, extract the UP element, -35 hexamer, spacer, -10 hexamer,
discriminator and initial transcribed region (ITR) upstream/downstream of
that TSS, score each combination with the linear free-energy model, and
keep the minimum-dG_total (= maximum Tx_rate) combination per TSS. Both
strands are scanned; a circular molecule is the caller's responsibility
(wrap the sequence before calling, per spec section 10.1).

  promoter_span = [up_start, tss + ITR_length)
  min promoter footprint = 24 (UP) + 1 (gap) + 6 (-35) + 15 (spacer)
                          + 6 (-10) + 6 (disc) + 20 (ITR) = 78 nt
  max promoter footprint = ... + 20 (spacer) + 10 (disc) = 89 nt
"""

from dataclasses import replace
from typing import Dict, List, Optional, Tuple

from . import coefficients as coef
from .energy_model import SigmaFactorElements, score_elements, tx_rate
from .types import PromoterHit, PromoterScanResult, Span

UP_LENGTH = 24
UP_HEX35_GAP = 1
HEX35_LENGTH = 6
SPACER_LENGTH_RANGE = (15, 21)   # upper bound exclusive: lengths 15..20
HEX10_LENGTH = 6
DISC_LENGTH_RANGE = (6, 11)      # upper bound exclusive: lengths 6..10
ITR_LENGTH = 20

_REVCOMP = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def _revcomp(seq: str) -> str:
    return seq.translate(_REVCOMP)[::-1]


def _flip_span(span: Span, length: int) -> Span:
    start, end = span
    return (length - end, length - start)


class PromoterCalculator:
    """Stateless sigma70 promoter scorer for one organism/condition.

    Only two Boltzmann constants are calibrated in the published model:
    ``organism="ecoli"`` (in vivo, *E. coli* MG1655) and
    ``organism="in_vitro"``. Every other host inherits the *E. coli*
    hexamer/UP-element/ITR coefficients as a documented prior, exactly as
    RBSForge documents for its own non-*E. coli* HostPacks -- see
    ``docs/MODEL.md``.
    """

    def __init__(self, organism: str = "ecoli"):
        self.organism = organism
        self.k = coef.K
        self.beta = coef.BETA_IN_VIVO if organism == "ecoli" else coef.BETA_IN_VITRO

    def _predict_one_strand(self, sequence: str, tss_range: Tuple[int, int]) -> Dict[int, PromoterHit]:
        hits: Dict[int, PromoterHit] = {}
        seq_len = len(sequence)

        for tss in range(tss_range[0], tss_range[1]):
            if tss + ITR_LENGTH > seq_len:
                continue
            itr = sequence[tss:tss + ITR_LENGTH]
            best: Optional[PromoterHit] = None

            for disc_length in range(*DISC_LENGTH_RANGE):
                if tss - disc_length < 0:
                    continue
                disc = sequence[tss - disc_length:tss]
                disc_span = (tss - disc_length, tss)

                for spacer_length in range(*SPACER_LENGTH_RANGE):
                    hex10_end = tss - disc_length
                    hex10_start = hex10_end - HEX10_LENGTH
                    spacer_end = hex10_start
                    spacer_start = spacer_end - spacer_length
                    hex35_end = spacer_start
                    hex35_start = hex35_end - HEX35_LENGTH
                    up_end = hex35_start - UP_HEX35_GAP
                    up_start = up_end - UP_LENGTH
                    if up_start < 0:
                        continue

                    hex10 = sequence[hex10_start:hex10_end]
                    spacer = sequence[spacer_start:spacer_end]
                    hex35 = sequence[hex35_start:hex35_end]
                    up = sequence[up_start:up_end]

                    elements = SigmaFactorElements(up=up, hex35=hex35, spacer=spacer, hex10=hex10, disc=disc, itr=itr)
                    breakdown = score_elements(elements)
                    dg_total = breakdown.total
                    if best is not None and dg_total >= best.dg_total:
                        continue

                    best = PromoterHit(
                        tss=tss, strand=1,
                        tx_rate=tx_rate(dg_total, organism=self.organism),
                        dg_total=dg_total, breakdown=breakdown,
                        promoter_span=(up_start, tss + ITR_LENGTH),
                        up_span=(up_start, up_end), hex35_span=(hex35_start, hex35_end),
                        spacer_span=(spacer_start, spacer_end), hex10_span=(hex10_start, hex10_end),
                        disc_span=disc_span, itr_span=(tss, tss + ITR_LENGTH),
                        up=up, hex35=hex35, spacer=spacer, hex10=hex10, disc=disc, itr=itr,
                    )

            if best is not None:
                hits[tss] = best

        return hits

    def predict(self, sequence: str, tss_range: Optional[Tuple[int, int]] = None) -> PromoterScanResult:
        """Scan both strands of ``sequence`` for the best sigma70 promoter
        state at every candidate TSS in ``tss_range`` (default: the whole
        sequence)."""
        sequence = sequence.upper().replace("U", "T")
        seq_len = len(sequence)
        if tss_range is None:
            tss_range = (0, seq_len)

        forward = self._predict_one_strand(sequence, tss_range)

        rev_sequence = _revcomp(sequence)
        rev_tss_range = (seq_len - tss_range[1], seq_len - tss_range[0])
        reverse_local = self._predict_one_strand(rev_sequence, rev_tss_range)

        reverse: Dict[int, PromoterHit] = {}
        for local_tss, hit in reverse_local.items():
            fwd_tss = seq_len - local_tss
            reverse[fwd_tss] = replace(
                hit, tss=fwd_tss, strand=-1,
                promoter_span=_flip_span(hit.promoter_span, seq_len),
                up_span=_flip_span(hit.up_span, seq_len),
                hex35_span=_flip_span(hit.hex35_span, seq_len),
                spacer_span=_flip_span(hit.spacer_span, seq_len),
                hex10_span=_flip_span(hit.hex10_span, seq_len),
                disc_span=_flip_span(hit.disc_span, seq_len),
                itr_span=_flip_span(hit.itr_span, seq_len),
            )

        return PromoterScanResult(
            sequence=sequence, organism=self.organism, k=self.k, beta=self.beta,
            forward=forward, reverse=reverse,
        )


def predict(sequence: str, organism: str = "ecoli",
            tss_range: Optional[Tuple[int, int]] = None) -> PromoterScanResult:
    """Convenience wrapper: ``PromoterCalculator(organism).predict(...)``."""
    return PromoterCalculator(organism).predict(sequence, tss_range)


def scan_promoters(dna: str, organism: str = "ecoli", tss_range: Optional[Tuple[int, int]] = None,
                    intended_tss: Optional[int] = None, intended_strand: int = 1,
                    tau_tx: float = 0.0) -> List[PromoterHit]:
    """Operon Calculator Module E entry point (spec section 25):
    ``scan_promoters(dna, host) -> list[PromoterHit]``.

    Returns hits ranked strongest first. If ``intended_tss`` is given,
    only "internal" (cryptic) hits are returned -- those that are not the
    intended promoter and clear ``tau_tx``, matching section 10.3. If
    ``intended_tss`` is ``None`` (Predict/Evaluate mode with no designated
    promoter), every hit above ``tau_tx`` is returned.
    """
    result = predict(dna, organism=organism, tss_range=tss_range)
    if intended_tss is None:
        return [h for h in result.ranked() if h.tx_rate > tau_tx]
    return result.internal(intended_tss=intended_tss, intended_strand=intended_strand, tau_tx=tau_tx)
