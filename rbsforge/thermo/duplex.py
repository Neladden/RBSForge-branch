"""SD:anti-SD duplex hybridization, jointly optimized with spacing.

Reimplements the core of the RBS Calculator v1.0 algorithm (model summary
section 4.2): search every plausible SD:anti-SD alignment (register and
extent), and for each compute both its hybridization free energy and its
*aligned spacing* to the start codon (the Chen/Salis-corrected spacer
length, not a raw nucleotide count -- see below). The alignment kept is
the one that minimizes delta_G_hybrid + delta_G_spacing(aligned_spacing)
jointly, exactly as the source algorithm does ("keep the duplex that
minimizes ΔG_mRNA-rRNA + ΔG_spacing").

Aligned spacing. The 16S rRNA's anti-SD tail is a fixed-length, fixed
sequence (`anti_sd`, given 5'->3', typically the last 9 nt of 16S rRNA).
Real SD sequences frequently pair with only part of it, and not always
flush against its 3' terminus. If the alignment's most-upstream paired
base (the mRNA window's 5'-most engaged position) pairs an anti-SD
position short of the tail's true 3' terminus, part of the tail dangles
unused. That unused length ("overhang") is subtracted from the raw
nucleotide spacer count to get the aligned spacing: an alignment that
uses the tail right up to its 3' end (overhang = 0) has aligned spacing
== raw spacing; one that stops short is treated as effectively more
compressed, since less of the tail needed to be "spent" reaching the
duplex.

Algorithm (deliberately simplified relative to a full NuPACK/ViennaRNA
cofold, in the interest of a transparent, dependency-free, from-scratch
reimplementation -- see docs/MODEL.md): for every relative register of
`anti_sd` against `mrna_window`, classify each aligned base pair as
Watson-Crick, G-U wobble, or mismatch, then brute-force every contiguous
sub-window (bounded by a minimum length and a maximum number of internal
mismatches/bulges), scoring each by nearest-neighbor stacking sums plus
helix-initiation and terminal-AU penalties. Because both strands are
short (anti_sd <= ~15 nt), this is cheap.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional, Tuple

from ..constants import REFERENCE_TEMPERATURE_K
from ..seqtools import is_watson_crick, is_wobble, validate_rna
from .nn_params import (
    helix_initiation,
    internal_mismatch_penalty,
    terminal_au_penalty,
    watson_crick_stack,
    wobble_penalty,
)

WC, WOBBLE, MISMATCH = "wc", "wobble", "mismatch"


def _pair_type(mrna_base: str, antisd_base: str) -> str:
    if is_watson_crick(mrna_base, antisd_base):
        return WC
    if is_wobble(mrna_base, antisd_base):
        return WOBBLE
    return MISMATCH


def _is_au_or_gu(mrna_base: str, antisd_base: str) -> bool:
    pair = mrna_base + antisd_base
    return pair in ("AU", "UA", "GU", "UG")


@dataclass
class DuplexAlignment:
    """One candidate SD:anti-SD hybridization, with its spacing context.

    `mrna_start`/`mrna_end` are 0-based, end-exclusive indices into the
    mRNA window that was searched.
    """

    delta_g_hybrid: float
    mrna_start: int
    mrna_end: int
    mrna_segment: str
    antisd_segment: str
    n_watson_crick: int
    n_wobble: int
    n_mismatch: int
    raw_spacing: int
    overhang: int
    aligned_spacing: int

    @property
    def length(self) -> int:
        return self.mrna_end - self.mrna_start


_NO_ALIGNMENT = DuplexAlignment(
    delta_g_hybrid=0.0,
    mrna_start=0,
    mrna_end=0,
    mrna_segment="",
    antisd_segment="",
    n_watson_crick=0,
    n_wobble=0,
    n_mismatch=0,
    raw_spacing=0,
    overhang=0,
    aligned_spacing=0,
)


def best_hybridization(
    mrna_window: str,
    anti_sd: str,
    start_codon_index: int,
    spacing_energy_fn: Callable[[int], float],
    temperature_k: float = REFERENCE_TEMPERATURE_K,
    min_length: int = 4,
    max_internal_mismatches: int = 1,
) -> Tuple[DuplexAlignment, float]:
    """Search all registers/sub-windows for the alignment minimizing
    delta_G_hybrid + spacing_energy_fn(aligned_spacing).

    Returns (best_alignment, combined_delta_g). If no window of at least
    `min_length` paired bases can be formed anywhere, returns
    (a zero/empty DuplexAlignment, spacing_energy_fn(no SD available)) --
    callers should treat that as "no SD-like sequence found" (a leaderless
    or SD-independent start).
    """
    mrna_window = validate_rna(mrna_window)
    anti_sd = validate_rna(anti_sd)
    n, k = len(mrna_window), len(anti_sd)
    if n == 0 or k == 0:
        return _NO_ALIGNMENT, spacing_energy_fn(start_codon_index)

    best: Optional[DuplexAlignment] = None
    best_total = float("inf")

    for offset in range(-(k - 1), n):
        lo = max(0, offset)
        hi = min(n, offset + k) - 1
        if hi - lo + 1 < min_length:
            continue

        antisd_index = lambda i: (k - 1) - (i - offset)  # noqa: E731
        pair_types: List[str] = [
            _pair_type(mrna_window[i], anti_sd[antisd_index(i)]) for i in range(lo, hi + 1)
        ]

        for start in range(len(pair_types)):
            if pair_types[start] == MISMATCH:
                continue
            for end in range(start, len(pair_types)):
                if pair_types[end] == MISMATCH:
                    continue
                length = end - start + 1
                if length < min_length:
                    continue
                mismatches = pair_types[start : end + 1].count(MISMATCH)
                if mismatches > max_internal_mismatches:
                    break  # extending `end` further only adds more bases

                candidate = _build_alignment(
                    mrna_window,
                    anti_sd,
                    offset,
                    lo + start,
                    lo + end,
                    pair_types[start : end + 1],
                    start_codon_index,
                    temperature_k,
                )
                total = candidate.delta_g_hybrid + spacing_energy_fn(candidate.aligned_spacing)
                if total < best_total:
                    best, best_total = candidate, total

    if best is None:
        return _NO_ALIGNMENT, spacing_energy_fn(start_codon_index)
    return best, best_total


def _build_alignment(
    mrna_window: str,
    anti_sd: str,
    offset: int,
    i_start: int,
    i_end: int,
    pair_types: List[str],
    start_codon_index: int,
    temperature_k: float,
) -> DuplexAlignment:
    k = len(anti_sd)
    antisd_index = lambda i: (k - 1) - (i - offset)  # noqa: E731

    energy = helix_initiation(temperature_k)
    for i in range(i_start, i_end):
        t0, t1 = pair_types[i - i_start], pair_types[i - i_start + 1]
        if t0 == WC and t1 == WC:
            m0, m1 = mrna_window[i], mrna_window[i + 1]
            a0, a1 = anti_sd[antisd_index(i)], anti_sd[antisd_index(i + 1)]
            energy += watson_crick_stack(f"{m0}{m1}/{a0}{a1}", temperature_k)
        elif t0 == MISMATCH or t1 == MISMATCH:
            energy += internal_mismatch_penalty(temperature_k)
        else:
            energy += wobble_penalty(temperature_k)

    if _is_au_or_gu(mrna_window[i_start], anti_sd[antisd_index(i_start)]):
        energy += terminal_au_penalty(temperature_k)
    if _is_au_or_gu(mrna_window[i_end], anti_sd[antisd_index(i_end)]):
        energy += terminal_au_penalty(temperature_k)

    antisd_segment = "".join(anti_sd[antisd_index(i)] for i in range(i_start, i_end + 1))

    # "how far from the 3' end of 16S that pair sits": the aligned window's
    # most-upstream (mRNA 5'-most) engaged base pairs an anti-SD position
    # `overhang` nucleotides short of the tail's true 3' terminus.
    overhang = i_start - offset
    raw_spacing = start_codon_index - (i_end + 1)
    aligned_spacing = raw_spacing - overhang

    return DuplexAlignment(
        delta_g_hybrid=energy,
        mrna_start=i_start,
        mrna_end=i_end + 1,
        mrna_segment=mrna_window[i_start : i_end + 1],
        antisd_segment=antisd_segment,
        n_watson_crick=pair_types.count(WC),
        n_wobble=pair_types.count(WOBBLE),
        n_mismatch=pair_types.count(MISMATCH),
        raw_spacing=raw_spacing,
        overhang=overhang,
        aligned_spacing=aligned_spacing,
    )
