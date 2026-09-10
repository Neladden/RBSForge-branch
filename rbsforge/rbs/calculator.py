"""Predict-mode RBS Calculator: mRNA sequence + HostPack + temperature -> TIR.

Reimplements the reconstructed v1.0-class Predict-mode algorithm (model
summary section 7) end to end:

For each start codon in the transcript (AUG/GUG/UUG by default):
  1. Slice a window of mRNA around the start codon (5' UTR up to the
     start, plus `footprint_cds`/`cutoff_post` nt into the CDS).
  2. Fold the window -> delta_G_mRNA (MFE, unconstrained).
  3. Search SD:anti-SD hybridization registers, jointly minimizing
     delta_G_hybrid + delta_G_spacing(aligned_spacing) -> delta_G_mRNA:rRNA
     and delta_G_spacing.
  4. Compute delta_G_standby (v1.0 4-nt forced-unpaired method) on the
     same window.
  5. Look up delta_G_start from the HostPack's start-codon table.
  6. delta_G_stacking: not implemented (v2.1 coefficient is unpublished --
     model summary section 4.6); reported as an explicit 0.0 term.
  7. Sum delta_G_total; r = exp(-beta * delta_G_total).

Start codons with no SD-like sequence found anywhere in the search window
are treated as leaderless and excluded from the ranked results (model
summary section 7, "Leaderless starts ... are skipped with a warning"),
but are still reported in `PredictResult.skipped` rather than silently
dropped.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Optional

from ..constants import celsius_to_kelvin
from ..hostpack import HostPack, get_hostpack
from ..scoring import ScoreBreakdown, proportional_rate, v1_style_rate
from ..seqtools import validate_rna
from ..thermo.duplex import DuplexAlignment, best_hybridization
from ..thermo.fold import BuiltinFolder, best_available_folder
from .spacing import spacing_penalty
from .standby import standby_penalty
from .start_codons import DEFAULT_SCANNED_CODONS, find_start_codons, start_codon_energy

# Safety cap on how much 5' UTR is pulled into the folded/searched window,
# independent of HostPack.cutoff_pre. The model's own default (35 nt) is
# far below this; this cap only guards against pathological runtimes on
# very long input transcripts when cutoff_pre is left unset (meaning "use
# the whole available 5' UTR" -- model summary section 6/13).
MAX_PRE_WINDOW = 200


@dataclass
class StartCodonResult:
    index: int  # 0-based index into the input sequence
    codon: str
    leaderless: bool
    breakdown: Optional[ScoreBreakdown] = None
    translation_initiation_rate: Optional[float] = None
    v1_style_rate: Optional[float] = None
    duplex: Optional[DuplexAlignment] = None
    aligned_spacing: Optional[int] = None
    raw_spacing: Optional[int] = None
    standby_window: Optional[tuple] = None
    window_start: int = 0
    window_end: int = 0
    footprint_span: Optional[tuple] = None  # (start, start + footprint_cds) on the input sequence

    @property
    def delta_g_total(self) -> Optional[float]:
        return self.breakdown.total if self.breakdown is not None else None


@dataclass
class PredictResult:
    sequence: str
    hostpack: HostPack
    temperature_c: float
    results: List[StartCodonResult] = field(default_factory=list)

    @property
    def best(self) -> Optional[StartCodonResult]:
        scored = [r for r in self.results if not r.leaderless]
        if not scored:
            return None
        return max(scored, key=lambda r: r.translation_initiation_rate)

    def ranked(self) -> List[StartCodonResult]:
        scored = [r for r in self.results if not r.leaderless]
        return sorted(scored, key=lambda r: r.translation_initiation_rate, reverse=True)


class RBSCalculator:
    """Predict-mode RBS strength calculator for a given HostPack."""

    def __init__(
        self,
        hostpack: HostPack,
        temperature_c: Optional[float] = None,
        scanned_codons: Iterable[str] = DEFAULT_SCANNED_CODONS,
        min_duplex_length: int = 4,
        max_internal_mismatches: int = 1,
        use_vienna_if_available: bool = True,
    ):
        self.hostpack = hostpack
        self.temperature_c = temperature_c if temperature_c is not None else hostpack.t_growth_c
        self.temperature_k = celsius_to_kelvin(self.temperature_c)
        self.scanned_codons = tuple(scanned_codons)
        self.min_duplex_length = min_duplex_length
        self.max_internal_mismatches = max_internal_mismatches
        self.folder = (
            best_available_folder(self.temperature_k)
            if use_vienna_if_available
            else BuiltinFolder(self.temperature_k)
        )

    @classmethod
    def for_species(cls, species: str, temperature_c: Optional[float] = None, **kwargs) -> "RBSCalculator":
        return cls(get_hostpack(species), temperature_c=temperature_c, **kwargs)

    def predict(self, mrna_sequence: str) -> PredictResult:
        sequence = validate_rna(mrna_sequence)
        hp = self.hostpack
        results: List[StartCodonResult] = []

        for index, codon in find_start_codons(sequence, self.scanned_codons):
            results.append(self.predict_one(sequence, index))

        return PredictResult(
            sequence=sequence,
            hostpack=hp,
            temperature_c=self.temperature_c,
            results=results,
        )

    def predict_one(
        self,
        mrna_sequence: str,
        start: int,
        extra_unpaired_global: Optional[Iterable[int]] = None,
    ) -> StartCodonResult:
        """Score one start codon, with an optional set of additional global
        (whole-`mrna_sequence`-indexed) positions forced unpaired in the
        bound-state fold.

        This is the one constrained-fold hook CDS-footprint unfolding and
        translational coupling both share (model summary section 6.5):
        the default occupied set is `[start, start + 3 + footprint_cds)`
        (the start codon plus HostPack.footprint_cds nt of 30S occupancy
        on the CDS); `extra_unpaired_global` -- coupling's upstream-CDS
        nucleotides that fall inside this window -- is unioned into it.
        Do not add a second, parallel constrained-fold code path; extend
        this one.
        """
        sequence = validate_rna(mrna_sequence)
        codon = sequence[start : start + 3]
        hp = self.hostpack
        pre_available = start
        pre_len = hp.cutoff_pre if hp.cutoff_pre is not None else pre_available
        pre_len = min(pre_len, pre_available, MAX_PRE_WINDOW)
        window_start = start - pre_len
        window_end = min(len(sequence), start + 3 + hp.cutoff_post)
        window = sequence[window_start:window_end]
        start_in_window = start - window_start

        utr_in_window = window[:start_in_window]

        def spacing_fn(aligned_spacing: int) -> float:
            return spacing_penalty(aligned_spacing, hp.s_opt, hp.spacing_push, hp.spacing_pull)

        duplex, _combined = best_hybridization(
            utr_in_window,
            hp.asd,
            start_codon_index=start_in_window,
            spacing_energy_fn=spacing_fn,
            temperature_k=self.temperature_k,
            min_length=self.min_duplex_length,
            max_internal_mismatches=self.max_internal_mismatches,
        )

        if duplex.length == 0:
            return StartCodonResult(
                index=start,
                codon=codon,
                leaderless=True,
                window_start=window_start,
                window_end=window_end,
            )

        n = len(window)
        footprint_end_in_window = min(n, start_in_window + 3 + hp.footprint_cds)
        occupied = set(range(start_in_window, footprint_end_in_window))
        if extra_unpaired_global:
            occupied |= {
                i - window_start for i in extra_unpaired_global if window_start <= i < window_end
            }

        e_initial = self.folder.fold(window).delta_g
        e_bound = self.folder.fold(window, forced_unpaired=frozenset(occupied)).delta_g
        # Constraining a fold can only raise (or leave unchanged) its minimum
        # free energy, so e_bound >= e_initial always; this is the physical
        # cost (>= 0) of reaching the bound-state configuration from the
        # unconstrained equilibrium, not "-1 * unconstrained MFE" renamed.
        unfolding = e_bound - e_initial

        delta_g_spacing = spacing_penalty(
            duplex.aligned_spacing, hp.s_opt, hp.spacing_push, hp.spacing_pull
        )
        delta_g_start = start_codon_energy(codon, hp.start_codon_dg)
        # Standby stays on the plain unconstrained-window baseline (not the
        # footprint-constrained one): occupied is deliberately just the
        # start codon + footprint_cds, so it never overlaps the standby
        # site's 4 nt, and this term does not double-count against it.
        delta_g_standby, standby_start, standby_end = standby_penalty(
            window, duplex.mrna_start, hp.standby.standby_site_nt, self.folder
        )

        breakdown = ScoreBreakdown()
        breakdown.add(
            "standby",
            delta_g_standby,
            f"v1.0 4-nt method, window [{standby_start}:{standby_end}) of the folded region",
        )
        breakdown.add(
            "mRNA:rRNA",
            duplex.delta_g_hybrid,
            f"SD '{duplex.mrna_segment}' : aSD '{duplex.antisd_segment}' "
            f"({duplex.n_watson_crick} WC, {duplex.n_wobble} wobble, {duplex.n_mismatch} mismatch)",
        )
        breakdown.add(
            "spacing",
            delta_g_spacing,
            f"aligned spacing = {duplex.aligned_spacing} nt (s_opt = {hp.s_opt})",
        )
        breakdown.add("start", delta_g_start, f"codon = {codon}")
        breakdown.add("stacking", 0.0, "v2.1 homopolymer term; coefficient unpublished, not implemented")
        breakdown.add(
            "mRNA",
            unfolding,
            f"bound-state unfolding cost: E_bound({e_bound:.2f}) - E_initial({e_initial:.2f}), "
            f"footprint {hp.footprint_cds} nt + start codon constrained unpaired",
        )

        rate = proportional_rate(breakdown.total, hp.beta)
        rate_v1 = v1_style_rate(breakdown.total, rt_eff=1.0 / hp.beta)

        return StartCodonResult(
            index=start,
            codon=codon,
            leaderless=False,
            breakdown=breakdown,
            translation_initiation_rate=rate,
            v1_style_rate=rate_v1,
            duplex=duplex,
            aligned_spacing=duplex.aligned_spacing,
            raw_spacing=duplex.raw_spacing,
            standby_window=(standby_start, standby_end),
            window_start=window_start,
            window_end=window_end,
            footprint_span=(start, start + hp.footprint_cds),
        )
