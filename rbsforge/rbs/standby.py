"""delta_G_standby: the v1.0 four-nucleotide standby-site penalty.

Model summary section 4.5 (v1.0): the standby site is the 4 nt
immediately upstream (5') of the SD:aSD duplex, where the 30S subunit
first contacts the mRNA before sliding to the start codon. The penalty is
the energetic cost of forcing those 4 nt to stay single-stranded (i.e. not
sequestered in a competing structure) in the folded mRNA window:

    delta_G_standby_v1 = min(0, E_before - E_after_forcing_unpaired)
    delta_G_total += -delta_G_standby_v1   # a negative unfolding cost becomes a positive penalty

`E_before` is the unconstrained MFE of the window; `E_after` is the MFE
of the same window with the standby-site positions forbidden from
pairing. If forcing them unpaired doesn't change (or improves) the
energy, there is no penalty (the clipping to <= 0 also means this term
can never award a bonus for a well-exposed standby site).

This is the v1.0 formula, chosen as the minimal-viable default (model
summary section 17); the richer v2.0 module geometry (proximal/distal
ssRNA binding site lengths, hairpin height, sliding penalty) is a
documented extension point (`StandbyParams` in hostpack.py already
carries its coefficients) not implemented here.
"""
from __future__ import annotations

from typing import Tuple


def standby_penalty(
    utr_window: str,
    sd_duplex_start: int,
    standby_site_nt: int,
    folder,
) -> Tuple[float, int, int]:
    """Returns (penalty_kcal_mol >= 0, standby_start, standby_end) for the
    `standby_site_nt` positions immediately upstream of `sd_duplex_start`
    (0-based index into `utr_window` of the duplex's 5'-most paired base).

    If there is no room for a standby site (the duplex starts at or before
    position 0), returns (0.0, 0, 0) -- no penalty is charged for a
    situation this simplified model cannot represent, rather than guessing.
    """
    standby_end = sd_duplex_start
    standby_start = max(0, standby_end - standby_site_nt)
    if standby_end <= 0 or standby_start >= standby_end:
        return 0.0, standby_start, standby_end

    e_before = folder.fold(utr_window).delta_g
    forced = frozenset(range(standby_start, standby_end))
    e_after = folder.fold(utr_window, forced_unpaired=forced).delta_g

    delta_v1 = min(0.0, e_before - e_after)
    penalty = -delta_v1
    return penalty, standby_start, standby_end
