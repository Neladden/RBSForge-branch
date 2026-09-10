"""
Operon Calculator -- Module B: translational coupling (spec section 7).

The only translation physics not already covered by RBSForge. Source:
Tian & Salis, NAR 43:7137 (2015), https://doi.org/10.1093/nar/gkv635.

Requires ``rbsforge`` to expose ``predict_one(mrna, start,
extra_unpaired_global=None)`` (spec section 6.5) and a Vienna-backed
folder -- the builtin folder cannot represent bulge-containing intergenic
hairpins (see ``docs/OPERON_CALCULATOR.md``, known failure mode 11).

Status: planned, not yet implemented. RBSForge does not yet expose
``predict_one``; this module is blocked on that first.
"""


def k_reinitiation(d: int) -> float:
    """Re-initiation coefficient as a function of intergenic distance
    ``d`` (Tian & Salis 2015, section 7.2). This piecewise fit is a ratio
    of rates and transfers across RBS-calculator au scales unchanged --
    it does not need re-fitting on RBSForge.

        d == -4          -> 0.022      (AUGA, tRNAfMet 4-nt pairing, no scanning)
        0 <= d <= 25      -> 0.0072    (forward scan, ~flat)
        d > 25            -> 0.0072 * 0.6 ** ((d - 25) / 100)
        d <= -25          -> 0.0072 / 11.6   (reverse-scan collision)
        -25 < d < 0, d!=-4 -> linear interpolation between the -25 and -4 points
    """
    if d == -4:
        return 0.022
    if 0 <= d <= 25:
        return 0.0072
    if d > 25:
        return 0.0072 * (0.6 ** ((d - 25) / 100.0))
    if d <= -25:
        return 0.0072 / 11.6
    return 0.0072 / 11.6 + (0.0072 - 0.0072 / 11.6) * (d + 25) / 25.0


def coupled_tirs(assembled, host, calc) -> list:
    """Nested 5'->3' iteration computing each CDS's coupled TIR as
    ``r_reinit + r_denovo`` (spec section 7.3):

        r_reinit_i = host.k_p * k_reinitiation(d) * r_{i-1}
        f          = min(1, host.c_unfold * to_physical(r_{i-1}, host))
        r_denovo_i = (1-f) * r_folded + f * r_unfolded
        r_i        = r_reinit_i + r_denovo_i

    where ``r_folded``/``r_unfolded`` come from two calls to
    ``calc.predict_one`` on the assembled mRNA -- the constrained-fold
    hook this module shares with CDS-footprint unfolding.

    Not a closed-form inverse: do not try to encode a target TIR ratio by
    inverting this function (section 7.5).
    """
    raise NotImplementedError(
        "Module B (coupling) is a prepared subsection, not yet implemented; "
        "blocked on rbsforge.RBSCalculator.predict_one"
    )
