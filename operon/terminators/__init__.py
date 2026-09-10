"""
Operon Calculator -- Module H: internal transcriptional terminators
(spec section 13).

There is no Salis Terminator Calculator; this is a scan plus a parts-
database append, not a fitted rate model of its own. Intrinsic
(rho-independent) detection follows Kingsford, Ayanbule & Salzberg,
Genome Biology 8:R22 (2007), https://doi.org/10.1186/gb-2007-8-2-r22, with
strength optionally scored per Chen et al., Nat. Methods 10:659 (2013),
https://doi.org/10.1038/nmeth.2515. Rho-dependent detection follows
RhoTermPredict (Di Salvo et al. 2019,
https://doi.org/10.1186/s12859-019-2704-x).

Do not put terminator strength into the *E. coli* mRNA half-life formula
(Module D) -- Cetnar & Salis 2021 found no effect on upstream mRNA.

Status: planned, not yet implemented.
"""


def scan_intrinsic_terminators(dna: str) -> list:
    """Rho-independent (hairpin + U-tract) terminator hits, both strands
    (section 13.2). A hit is "internal" if it is not the user-supplied 3'
    terminator (allow a ~50 bp window around the annotated terminator).
    """
    raise NotImplementedError("Module H (terminators) is a prepared subsection, not yet implemented")


def scan_rho_terminators(dna: str) -> list:
    """Rho-dependent terminator hits via Rut-site + pause detection,
    both strands (section 13.3)."""
    raise NotImplementedError("Module H (terminators) is a prepared subsection, not yet implemented")
