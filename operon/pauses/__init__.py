"""
Operon Calculator -- Module G: ribosomal pause sites (spec section 12).

Independent motif + structure scanner inside CDSs. The lab never
published a standalone Pause Calculator; the operational definition is
the union of internal Shine-Dalgarno / anti-SD hybridization (reuse
``rbsforge.thermo.duplex.best_hybridization``), slow-codon runs,
polyproline/stall motifs, and (optionally) stable hairpins in the
ribosomal E/P site window.

Status: planned, not yet implemented.
"""


def scan_pauses(cds_nt: str, host, calc) -> list:
    """Positions and scores of candidate ribosomal pause sites inside one
    CDS (section 12). Host-dependent: use ``host.pack.asd``, never a
    hardcoded *E. coli* anti-SD tail.
    """
    raise NotImplementedError("Module G (pauses) is a prepared subsection, not yet implemented")
