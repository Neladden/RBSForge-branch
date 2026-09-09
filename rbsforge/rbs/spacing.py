"""delta_G_spacing: the free-energy penalty for non-optimal aligned spacing.

Exact v1.0 formulas (model summary section 4.3). Let `ds = s - s_opt`
where `s` is the *aligned* spacing (see thermo/duplex.py), not a raw
nucleotide count.

Compression (s < s_opt), a logistic "push" penalty:

    delta_G_spacing = push[0] / (1 + exp(push[1] * (ds + push[2])))**push[3]

Stretch (s >= s_opt), a quadratic "pull" penalty:

    delta_G_spacing = pull[0] * ds**2 + pull[1] * ds + pull[2]

At s == s_opt (ds == 0), both branches agree to ~0. `push`/`pull` are
HostPack fields (`spacing_push`, `spacing_pull`); the E. coli-fit defaults
are push=[12.2, 2.5, 2.0, 3.0], pull=[0.048, 0.24, 0.0].
"""
from __future__ import annotations

import math
from typing import Tuple


def spacing_penalty(
    aligned_spacing: int,
    s_opt: int,
    push: Tuple[float, float, float, float],
    pull: Tuple[float, float, float],
) -> float:
    ds = aligned_spacing - s_opt
    if ds < 0:
        p0, p1, p2, p3 = push
        return p0 / (1 + math.exp(p1 * (ds + p2))) ** p3
    q0, q1, q2 = pull
    return q0 * ds ** 2 + q1 * ds + q2
