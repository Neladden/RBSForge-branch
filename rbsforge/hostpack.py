"""HostPack: the organism-specific parameter set for RBS strength prediction.

This is the extension point for scoring sequences from organisms other
than E. coli, following the "replication architecture" laid out in the
model summary (docs/MODEL.md section 13): the thermodynamic machinery
(nearest-neighbor RNA folding, SD:aSD duplex search, the six-term free
energy sum) is organism-agnostic. Only a small, explicit set of numbers is
organism-specific, and they are all collected here in one `HostPack`.

**What is well-supported by evidence for a new organism** (model summary
section 10 and 16): the 16S rRNA anti-SD tail sequence (`asd`), the growth
temperature (`t_growth_c`), and whether the ribosome uses a Gram-positive-
or Gram-negative-style standby platform. Swapping these and re-running the
same equations at the right temperature ("Layer 1" in the model summary)
is a one-day, defensible step.

**What requires real calibration data before you should trust absolute
numbers** (section 12, 14 "Layer 2"): beta (the free-energy-to-rate
conversion), the spacing optimum/curvature, and the standby-site
coefficients. Every non-E.-coli HostPack below inherits the E. coli
values for these as an explicit, documented *prior*, not a fitted result.
Do not report absolute translation initiation rates from an uncalibrated
HostPack as if they were validated; relative rankings *within* one
HostPack are far more defensible than absolute rates across HostPacks.

To build a HostPack for a new organism:

1. Determine the mature 16S rRNA 3' end. Do not trust a raw NCBI/GenBank
   16S annotation -- many are truncated before the anti-SD (Nakov et al.
   2018 RNA). Prefer RNA-seq 3'-end mapping, or the "helix-45 -> CCUCCU ->
   tail" rule cross-checked against a close relative with a mapped tail.
   Use the *whole* mature tail (often 9-15 nt), not a blind 9-mer crop.
2. Record the growth/expression temperature and Gram stain.
3. Leave beta, spacing, and standby coefficients at the E. coli prior
   until you have a calibration library (model summary section 14, Layer
   2) -- a few dozen insulated reporter constructs spanning SD strength,
   spacer length, start codon, and structure, in that host, at that
   temperature.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

from .seqtools import validate_rna

# Start codon initiation free energies (kcal/mol), v1.0 source values:
# apparent hybridization energy of each codon to the initiator fMet-tRNA
# anticodon (CAU). These are *apparent* values -- they implicitly include
# IF3 discrimination in vivo, not a pure RNA:RNA duplex energy -- fit on
# E. coli. AUG is the reference (most favorable); alternative starts are
# progressively less favorable, in the well-established relative order
# AUG > GUG > UUG (CUG is tabulated but excluded from v1.0's start-codon
# scan; included here for completeness).
DEFAULT_START_CODON_DG: Dict[str, float] = {
    "AUG": -1.194,
    "GUG": -0.0748,
    "UUG": -0.0435,
    "CUG": -0.03406,
}


@dataclass(frozen=True)
class StandbyParams:
    """v2.0-style standby-site distortion coefficients (model summary 4.5).

    delta_G_distortion = C1*A_s^2 + C2*A_s + C3, A_s = A0 + P + D - min(H, H_cap)
    """

    c1: float = 0.038
    c2: float = -1.629
    c3: float = 17.359
    c_slide: float = 0.20  # kcal/mol per nt of downstream hairpin height
    h_cap: int = 15
    a0: int = 15
    standby_site_nt: int = 4  # v1.0 4-nt forced-unpaired window, used by rbs/standby.py


@dataclass(frozen=True)
class HostPack:
    """Complete organism-specific parameter set for one prediction run."""

    name: str
    phylogeny: str
    gram_stain: str  # "negative", "positive", or "unknown"
    t_growth_c: float
    asd: str  # mature 16S rRNA 3' tail, RNA alphabet, 9-15 nt, 5'->3'
    asd_modifications: Dict[int, str] = field(default_factory=dict)  # {0-based pos: mod name}
    s_opt: int = 5
    spacing_push: Tuple[float, float, float, float] = (12.2, 2.5, 2.0, 3.0)
    spacing_pull: Tuple[float, float, float] = (0.048, 0.24, 0.0)
    standby: StandbyParams = field(default_factory=StandbyParams)
    start_codon_dg: Dict[str, float] = field(default_factory=lambda: dict(DEFAULT_START_CODON_DG))
    beta: float = 0.45  # mol/kcal
    nn_param_set: str = "xia1998/turner"
    footprint_cds: int = 13  # nt past the start codon that must stay unfolded
    cutoff_post: int = 35  # nt of CDS included in the folded window
    cutoff_pre: Optional[int] = None  # nt of 5' UTR included; None = entire UTR given
    rbs_max_len: int = 35
    notes: str = ""

    def __post_init__(self):
        object.__setattr__(self, "asd", validate_rna(self.asd))


# --- Built-in HostPacks ---------------------------------------------------

ECOLI_K12 = HostPack(
    name="Escherichia coli K-12",
    phylogeny="Gammaproteobacteria",
    gram_stain="negative",
    t_growth_c=37.0,
    asd="ACCUCCUUA",  # v1.0 source: rRNA = "acctcctta"
    s_opt=5,
    spacing_push=(12.2, 2.5, 2.0, 3.0),
    spacing_pull=(0.048, 0.24, 0.0),
    standby=StandbyParams(),
    start_codon_dg=dict(DEFAULT_START_CODON_DG),
    beta=0.45,
    nn_param_set="xia1998/turner (v1.0/v2.0 default; v2.1 uses Andronescu 2007)",
    notes=(
        "Reference organism: every coefficient here is the published v1.0/"
        "v2.0 fitted value (Salis, Mirsky & Voigt 2009; Espah Borujeni, "
        "Channarasappa & Salis 2013). All other HostPacks in this module "
        "inherit beta/spacing/standby from this one as an unfit prior."
    ),
)

BACILLUS_SUBTILIS = HostPack(
    name="Bacillus subtilis",
    phylogeny="Firmicutes (Bacilli)",
    gram_stain="positive",
    t_growth_c=37.0,
    # RNA-seq-mapped mature tail (model summary section 10): 15 nt, longer
    # than E. coli's -- using only a cropped 9-mer would discard real
    # pairing to the extra U-rich 3' extension.
    asd="GAUCACCUCCUUUCU",
    # Classical Gram-positive spacing optimum (Vellanoweth & Rabinowitz
    # 1992); public, but not the same as a v2.1-refit value (that
    # coefficient table is unpublished -- model summary section 4.3/16).
    s_opt=8,
    spacing_push=(12.2, 2.5, 2.0, 3.0),  # E. coli prior; curvature unrefit
    spacing_pull=(0.048, 0.24, 0.0),  # E. coli prior; curvature unrefit
    standby=StandbyParams(),  # E. coli S1-platform prior; Gram-pos coefficients unpublished
    start_codon_dg=dict(DEFAULT_START_CODON_DG),
    beta=0.45,
    nn_param_set="xia1998/turner",
    notes=(
        "aSD tail length/sequence is real (RNA-seq mapped, per model "
        "summary section 10); s_opt reflects the documented Gram-positive "
        "spacing literature. beta, spacing curvature, and standby "
        "coefficients are carried over from the E. coli fit as a "
        "documented prior, per the model summary's Layer-1 recommendation "
        "('keep E. coli beta, spacing, standby, start table as the prior; "
        "swap only aSD and T') -- refit before trusting absolute TIR."
    ),
)

GEOBACILLUS_THERMOGLUCOSIDASIUS = HostPack(
    name="Geobacillus/Parageobacillus thermoglucosidasius",
    phylogeny="Firmicutes (Bacilli), thermophile",
    gram_stain="positive",
    t_growth_c=60.0,
    # No confirmed RNA-seq-mapped tail for this genus at authoring time;
    # using the B. subtilis (related Firmicute) tail as the documented
    # starting point the model summary recommends (section 15), NOT a
    # verified Geobacillus sequence. Confirm by RNA-seq before real use.
    asd="GAUCACCUCCUUUCU",
    s_opt=8,
    spacing_push=(12.2, 2.5, 2.0, 3.0),
    spacing_pull=(0.048, 0.24, 0.0),
    standby=StandbyParams(),
    start_codon_dg=dict(DEFAULT_START_CODON_DG),
    beta=0.45,
    nn_param_set="xia1998/turner",
    notes=(
        "Placeholder profile for a thermophilic Firmicute (T_opt ~55-65C). "
        "asd is borrowed from B. subtilis (unconfirmed for this genus -- "
        "verify by RNA-seq 3'-end mapping before real designs). All "
        "kinetic/energetic coefficients besides T are the E. coli prior. "
        "This profile only performs the model summary's 'Layer 1' "
        "temperature-aware physics step, not a real calibration (Layer 2)."
    ),
)

THERMUS_THERMOPHILUS = HostPack(
    name="Thermus thermophilus",
    phylogeny="Deinococcus-Thermus, thermophile",
    gram_stain="negative",
    t_growth_c=70.0,
    asd="ACCUCCUUA",  # CCUCCU-core placeholder; confirm mature tail by RNA-seq
    # Pseudouridines at 16S positions 1540/1541 (Guymon et al. 2006) fall
    # within the anti-SD; flagged here (0-based positions from the 5' end
    # of `asd`) as a documented, unapplied caveat -- Psi stacks more
    # strongly than U, so scoring with unmodified U underestimates SD:aSD
    # stability for this organism (model summary section 12.3).
    asd_modifications={7: "pseudouridine", 8: "pseudouridine"},
    s_opt=5,  # unrefit; Thermus 30S lacks S21 and its platform differs from E. coli's
    spacing_push=(12.2, 2.5, 2.0, 3.0),
    spacing_pull=(0.048, 0.24, 0.0),
    standby=StandbyParams(),
    start_codon_dg=dict(DEFAULT_START_CODON_DG),
    beta=0.45,
    nn_param_set="xia1998/turner",
    notes=(
        "Highest-uncertainty built-in profile. asd sequence is an "
        "unconfirmed CCUCCU-core placeholder (native SDs in this species "
        "can be short and non-canonical, e.g. AGAGGG in thrS); the "
        "pseudouridine modifications at the tail's 3' end are noted but "
        "not numerically applied by the built-in folder/duplex search. "
        "s_opt/standby are the E. coli prior despite a documented "
        "different 30S platform (no S21). Do not trust absolute TIR "
        "without a Layer-2 calibration in this species."
    ),
)

GENERIC_UNCHARACTERIZED = HostPack(
    name="generic/uncharacterized",
    phylogeny="unknown",
    gram_stain="unknown",
    t_growth_c=37.0,
    asd="ACCUCCUUA",
    s_opt=5,
    standby=StandbyParams(),
    start_codon_dg=dict(DEFAULT_START_CODON_DG),
    beta=0.45,
    nn_param_set="xia1998/turner",
    notes=(
        "Placeholder using the E. coli prior throughout. Replace `asd` "
        "with the organism's actual mature 16S rRNA 3' tail at minimum "
        "before drawing any conclusions; treat every other field as "
        "unvalidated until a calibration library is run (model summary "
        "section 14)."
    ),
)

BUILTIN_HOSTPACKS: Dict[str, HostPack] = {
    "ecoli": ECOLI_K12,
    "e_coli": ECOLI_K12,
    "e_coli_k12": ECOLI_K12,
    "b_subtilis": BACILLUS_SUBTILIS,
    "bacillus_subtilis": BACILLUS_SUBTILIS,
    "geobacillus": GEOBACILLUS_THERMOGLUCOSIDASIUS,
    "parageobacillus": GEOBACILLUS_THERMOGLUCOSIDASIUS,
    "thermus": THERMUS_THERMOPHILUS,
    "thermus_thermophilus": THERMUS_THERMOPHILUS,
    "generic": GENERIC_UNCHARACTERIZED,
}


def get_hostpack(species: str) -> HostPack:
    """Look up a built-in HostPack by short name (case-insensitive)."""
    key = species.strip().lower().replace(" ", "_").replace(".", "").replace("-", "_")
    if key not in BUILTIN_HOSTPACKS:
        available = ", ".join(sorted(BUILTIN_HOSTPACKS))
        raise KeyError(f"Unknown species/HostPack {species!r}. Available: {available}")
    return BUILTIN_HOSTPACKS[key]
