"""Nearest-neighbor RNA duplex free-energy parameters, temperature-scaled.

The RBS Calculator's published equations (Reis & Salis 2020, section
"Computational stack") scale nearest-neighbor free energies to a target
temperature via the standard van't Hoff decomposition:

    delta_G(T) = delta_H - T * delta_S
    delta_S    = (delta_H - delta_G_37) / 310.15

This module stores (delta_G_37, delta_H) pairs -- the standard way
nearest-neighbor tables are published -- and derives delta_S/delta_G(T)
from them, exactly as described above. Values are the commonly-cited
unified RNA nearest-neighbor parameters (Xia et al. 1998, Biochemistry;
consistent with Turner 2004 rules for Watson-Crick stacks). This is the
"turner1999"-class table referenced as the v1.0/v2.0 default in the model
summary; swapping in a different published table (e.g. Andronescu 2007,
used by v2.1) is a matter of replacing the numbers below, not the code
structure.

Caveat (stated explicitly in the model documentation): heat-capacity
(delta_Cp) effects are neglected, i.e. delta_H and delta_S are treated as
temperature-independent. This is a good approximation up to ~60 C and a
progressively worse one above that -- errors of 1-3 kcal/mol per short
helix have been reported by 75-100 C (see docs/MODEL.md).
"""
from __future__ import annotations

import math

from ..constants import REFERENCE_TEMPERATURE_K, rt as _rt

# Watson-Crick nearest-neighbor stacking parameters, keyed as "XY/ZW"
# meaning:
#   5'-X Y-3'
#   3'-Z W-5'
# i.e. X pairs Z, Y pairs W. Values are (delta_G_37, delta_H) in kcal/mol.
_WATSON_CRICK_STACKS_37 = {
    "AA/UU": (-0.93, -6.82),
    "AU/UA": (-1.10, -9.38),
    "UA/AU": (-1.33, -7.69),
    "CU/GA": (-2.08, -10.48),
    "CA/GU": (-2.11, -10.44),
    "GU/CA": (-2.24, -11.40),
    "GA/CU": (-2.35, -12.44),
    "CG/GC": (-2.36, -10.64),
    "GG/CC": (-3.26, -13.39),
    "GC/CG": (-3.42, -14.88),
}


def _add_reverse_duplicates(table: dict) -> dict:
    out = dict(table)
    for key, value in table.items():
        top, bottom = key.split("/")
        rev_key = f"{bottom[::-1]}/{top[::-1]}"
        out.setdefault(rev_key, value)
    return out


_WATSON_CRICK_STACKS_37 = _add_reverse_duplicates(_WATSON_CRICK_STACKS_37)

# Helix initiation and terminal A-U/G-U penalty (entropy-dominated;
# treated as delta_H ~= 0, i.e. scaling with T only through the RT term
# used elsewhere, not through this module's van't Hoff machinery).
HELIX_INITIATION_37 = 4.09
TERMINAL_AU_PENALTY_37 = 0.45

# Simplified stand-ins (see thermo/duplex.py and thermo/fold.py docstrings)
# for interior-loop/bulge and G-U-wobble contributions that a full nearest-
# neighbor table would otherwise supply. Treated as delta_H ~= 0.
INTERNAL_MISMATCH_PENALTY_37 = 1.8
WOBBLE_PENALTY_37 = 0.5


def _scaled(dg37: float, dh: float, temperature_k: float) -> float:
    if dh == 0.0:
        return dg37
    ds = (dh - dg37) / REFERENCE_TEMPERATURE_K
    return dh - temperature_k * ds


def watson_crick_stack(key: str, temperature_k: float) -> float:
    """Temperature-scaled delta_G (kcal/mol) for one NN stack, e.g. 'GC/CG'."""
    dg37, dh = _WATSON_CRICK_STACKS_37.get(key, (0.0, 0.0))
    return _scaled(dg37, dh, temperature_k)


def helix_initiation(temperature_k: float) -> float:
    return HELIX_INITIATION_37  # delta_H ~= 0 approximation


def terminal_au_penalty(temperature_k: float) -> float:
    return TERMINAL_AU_PENALTY_37


def internal_mismatch_penalty(temperature_k: float) -> float:
    return INTERNAL_MISMATCH_PENALTY_37


def wobble_penalty(temperature_k: float) -> float:
    return WOBBLE_PENALTY_37


# --- Single-strand hairpin loop initiation -------------------------------
#
# Hairpin loop free energies are overwhelmingly entropic (loop closure
# restricts conformational freedom), so -- absent a full per-size
# delta_H/delta_S table -- they are scaled proportionally to T (i.e.
# delta_H ~= 0, delta_G(T) = delta_G_37 * T/310.15), which is the standard
# simplifying assumption for loop-initiation terms when only delta_G_37 is
# tabulated.
HAIRPIN_LOOP_INITIATION_37 = {
    3: 5.4,
    4: 5.6,
    5: 5.7,
    6: 5.4,
    7: 6.0,
    8: 5.5,
    9: 6.4,
}
_HAIRPIN_LOOP_REFERENCE_SIZE = 9
_HAIRPIN_LOOP_REFERENCE_DG37 = HAIRPIN_LOOP_INITIATION_37[_HAIRPIN_LOOP_REFERENCE_SIZE]


def hairpin_loop_energy(loop_size: int, temperature_k: float) -> float:
    """Free energy (kcal/mol) to initiate a hairpin loop of `loop_size` nt,
    at 37 C from the tabulated/extrapolated value, then entropically
    scaled to `temperature_k` (see module docstring)."""
    if loop_size < 3:
        return float("inf")  # sterically impossible
    if loop_size in HAIRPIN_LOOP_INITIATION_37:
        dg37 = HAIRPIN_LOOP_INITIATION_37[loop_size]
    else:
        dg37 = _HAIRPIN_LOOP_REFERENCE_DG37 + 1.75 * _rt(REFERENCE_TEMPERATURE_K) * math.log(
            loop_size / _HAIRPIN_LOOP_REFERENCE_SIZE
        )
    return dg37 * (temperature_k / REFERENCE_TEMPERATURE_K)
