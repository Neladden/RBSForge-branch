"""
The sigma70 promoter linear free-energy model (Appendix D / Module E of the
Operon Calculator spec):

    dG_10     = HEX10_LEFT[hex10[0:3]] + HEX10_RIGHT[hex10[3:6]]
    dG_35     = HEX35_LEFT[hex35[0:3]] + HEX35_RIGHT[hex35[3:6]]
    dG_ext10  = EXT10_2MER[spacer[-3:-1]]
    dG_disc   = DISC_FIRST3[disc[0:3]]
    dG_spacer = 0.1463*s**2 - 4.9113*s + 41.119          # s = len(spacer)
    dG_ITR    = coef_itr * (dg_hybrid / NORM_ITR)
    dG_UP     = coef_dist*(groove_dist/NORM_DIST_UP) + coef_prox*(groove_prox/NORM_PROX_UP)
              + coef_rigidity*(rigidity/NORM_RIGIDITY)
    dG_total  = dG_10 + dG_35 + dG_disc + dG_ITR + dG_ext10 + dG_spacer + dG_UP + INTERCEPT

    Tx_rate   = K * exp(-BETA * dG_total)

This is a direct, from-scratch re-expression of the published model: only
the fitted numeric coefficients (``coefficients.py``) and the published
dinucleotide parameter tables (``tables.py``) are reused; the control flow
below is independently written.
"""

from dataclasses import dataclass

from . import coefficients as coef
from .features import dna_rna_hybrid_energy, groove_width, rigidity


@dataclass(frozen=True)
class SigmaFactorElements:
    """The six sequence elements of one candidate sigma70 promoter state,
    in 5' to 3' order. ``spacer`` and ``disc`` (the discriminator) are
    variable length; every other element is fixed length."""

    up: str          # 24 nt UP element (distal 12 + proximal 12)
    hex35: str        # 6 nt -35 hexamer
    spacer: str        # 15-20 nt spacer (last 2 nt are the extended -10 / TGn motif)
    hex10: str         # 6 nt -10 hexamer
    disc: str          # 6-10 nt discriminator
    itr: str           # >=20 nt initial transcribed region


@dataclass(frozen=True)
class PromoterEnergyBreakdown:
    """Per-term free-energy contributions (kcal/mol-equivalent, on the
    model's own fitted scale) for one scored promoter state."""

    dg_10: float
    dg_35: float
    dg_disc: float
    dg_ext10: float
    dg_spacer: float
    dg_itr: float
    dg_up: float
    intercept: float

    @property
    def total(self) -> float:
        return (
            self.dg_10 + self.dg_35 + self.dg_disc + self.dg_ext10
            + self.dg_spacer + self.dg_itr + self.dg_up + self.intercept
        )

    def terms(self) -> dict[str, float]:
        return {
            "hex10": self.dg_10,
            "hex35": self.dg_35,
            "discriminator": self.dg_disc,
            "ext10": self.dg_ext10,
            "spacer": self.dg_spacer,
            "itr": self.dg_itr,
            "up": self.dg_up,
            "intercept": self.intercept,
        }


def score_elements(elements: SigmaFactorElements) -> PromoterEnergyBreakdown:
    """Score one fully-specified candidate promoter state.

    Raises ``KeyError`` if any fixed-length categorical slice contains a
    base outside ``ACGT`` (uppercase DNA is required; lowercase or
    ambiguity codes are not accepted by the fitted hexamer tables).
    """
    hex10, hex35, spacer, disc = elements.hex10, elements.hex35, elements.spacer, elements.disc
    up = elements.up

    dg_10 = coef.HEX10_LEFT[hex10[0:3]] + coef.HEX10_RIGHT[hex10[3:6]]
    dg_35 = coef.HEX35_LEFT[hex35[0:3]] + coef.HEX35_RIGHT[hex35[3:6]]
    dg_disc = coef.DISC_FIRST3[disc[0:3]]
    dg_ext10 = coef.EXT10_2MER[spacer[-3:-1]]

    s = float(len(spacer))
    dg_spacer = 0.1463 * s ** 2 - 4.9113 * s + 41.119

    _, _, dg_hybrid = dna_rna_hybrid_energy(elements.itr)
    dg_itr = coef.NUMERICAL_COEFS["itr"] * (dg_hybrid / coef.NORM_ITR)

    prox_up = up[-(len(up) // 2):]
    dist_up = up[0:len(up) // 2]
    width_dist = groove_width(dist_up)
    width_prox = groove_width(prox_up)
    up_rigidity = rigidity(up + hex35 + spacer[0:14])

    dg_up = (
        coef.NUMERICAL_COEFS["dist_up"] * (width_dist / coef.NORM_DIST_UP)
        + coef.NUMERICAL_COEFS["prox_up"] * (width_prox / coef.NORM_PROX_UP)
        + coef.NUMERICAL_COEFS["rigidity"] * (up_rigidity / coef.NORM_RIGIDITY)
    )

    return PromoterEnergyBreakdown(
        dg_10=dg_10, dg_35=dg_35, dg_disc=dg_disc, dg_ext10=dg_ext10,
        dg_spacer=dg_spacer, dg_itr=dg_itr, dg_up=dg_up, intercept=coef.INTERCEPT,
    )


def tx_rate(dg_total: float, *, organism: str = "ecoli") -> float:
    """``K * exp(-BETA * dG_total)``. ``organism`` selects the fitted
    Boltzmann constant: ``"ecoli"`` (in vivo, MG1655) or ``"in_vitro"``.
    Only these two are calibrated -- see ``docs/MODEL.md``."""
    import math

    beta = coef.BETA_IN_VIVO if organism == "ecoli" else coef.BETA_IN_VITRO
    return coef.K * math.exp(-beta * dg_total)
