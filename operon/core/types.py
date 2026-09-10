"""
Shared types for the Operon Calculator (spec section 2.1).

These are the coordinate system and record types every module operates
on. They intentionally carry no behavior -- each module (assembly,
coupling, promoter_calculator, ...) is free to read and populate them
without depending on any other module's internals.

Coordinate convention: assembled DNA is 5' to 3' on the coding strand.
mRNA is the same sequence with T to U from the promoter's TSS to the
terminator's 3' end. All scanners report ``[start, end)`` on that
molecule plus strand.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class StandbyParams:
    """v2.0 standby-site module-geometry coefficients (Espah Borujeni,
    Channarasappa & Salis, NAR 2014). Scaffolded on HostPack; RBSForge
    wires the v1.0 4-nt form only -- see rbsforge/docs/MODEL.md."""

    c1: float = 0.038
    c2: float = -1.629
    c3: float = 17.359
    c_slide: float = 0.20
    h_cap: float = 15.0
    a0: float = 15.0


@dataclass
class OperonHost:
    """Wraps an ``rbsforge.HostPack`` with the Operon-level, organism-
    specific fields the RBS calculator itself has no reason to know about
    (spec section 2.1, Appendix F). Do not fork HostPack; do not put these
    fields into it.
    """

    pack: object  # an rbsforge.HostPack -- typed loosely to avoid a hard rbsforge import here
    codon_table: Optional[dict] = None
    codon_weights_high: Optional[dict] = None
    codon_weights_balanced: Optional[dict] = None
    r_elong_nominal_nt_s: float = 60.0
    ribosome_footprint_nt: int = 30       # 70S coupling occupancy -- NOT HostPack.footprint_cds (13)
    c_unfold: float = 0.81                # Tian 2015 C; re-fit per host/au-scale before trusting
    k_p: float = 10.0                     # Tian 2015 k_P; re-fit per host/au-scale before trusting
    is_motifs: List[str] = field(default_factory=list)
    att_motifs: List[str] = field(default_factory=list)
    forbidden_re_sites: List[str] = field(default_factory=list)
    rnase_e_like: bool = True             # True for Gram-negative hosts (RNase E); False -> RNase J/Y


@dataclass
class CDS:
    id: str
    aa_or_nt: str
    target_tir: Optional[float] = None    # Design mode only, au (v1_style_rate scale)
    rbs_constraint: Optional[str] = None  # IUPAC, same length as the designed RBS
    annotated_start: Optional[int] = None  # 0-based on the assembled mRNA
    annotated_stop: Optional[int] = None


@dataclass
class Junction:
    between: tuple            # (upstream CDS id, downstream CDS id)
    stop_i: int
    start_i_plus_1: int
    d_nt: int                 # start_{i+1} - stop_i_end; negative = overlap


@dataclass
class Operon:
    promoter_dna: str
    cds_list: List[CDS]
    rbs_list: List[str]
    terminator_dna: str
    host: OperonHost
    intergenic_policy: str = "free"  # 'overlap-4' | 'overlap-1' | 'abut' | 'spacer-N' | 'free'


@dataclass
class Assembled:
    dna: str
    mrna: str                 # T -> U from TSS to terminator end
    features: List[Dict]      # [{id, type, start, end, strand}]
    junctions: List[Junction]
