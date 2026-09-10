"""
Operon Calculator -- Module D: mRNA stability (spec section 9).

Independent of coupling except that it consumes coupled TIRs as
ribosome-protection inputs. Implement the 2021 biophysical model first
(Cetnar & Salis, ACS Synth. Biol. 2021,
https://doi.org/10.1021/acssynbio.0c00471); the 2024 isoform-aware LightGBM
expansion (https://doi.org/10.1038/s41467-024-54059-7) is a later layer.

Decay in *E. coli* for these synthetic operons is 5'-end / ribosome-
protection dominated -- do not put terminator strength into the half-life
formula (section 9.1, known failure mode 6).

Status: planned, not yet implemented. RBSForge explicitly does not
include this; it is entirely this module's responsibility.
"""


def mrna_stability(assembled, tirs: list, host):
    """Ribosome-protection + RNase-accessibility half-life estimate
    (section 9.2)::

        dx = r_elong / r_init - L_footprint          # nt of naked mRNA between ribosomes
        RNase_score = sum(ss[i] * motif_w[nt[i]] for i in 5' UTR)   # unpaired AU-richness
        log(mRNA) = a0 + a1*(-RNase_score) + a2*(-dx) + a3*1[5p_hairpin]
        k_decay   = k0 * exp(b1*RNase_score + b2*dx)
        t_half    = ln(2) / k_decay

    {a, b} must be fit on a paired sequence/qPCR set (the 2021 paper is
    the calibration target); this function is not usable until they are.
    """
    raise NotImplementedError("Module D (stability) is a prepared subsection, not yet implemented")
