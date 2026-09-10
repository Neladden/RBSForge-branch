"""
Operon Calculator -- Module B: translational coupling (spec section 7).
**Status: implemented.**

The only translation physics not already covered by RBSForge. Source:
Tian & Salis, NAR 43:7137 (2015), https://doi.org/10.1093/nar/gkv635.

Requires ``rbsforge.RBSCalculator.predict_one(mrna, start,
extra_unpaired_global=None)`` (spec section 6.5) -- implemented -- and, by
default, a Vienna-backed folder: the builtin folder cannot represent
bulge-containing intergenic hairpins (known failure mode 11), so
``coupled_tirs`` raises unless ``require_vienna=False`` is passed
explicitly (development/testing only -- results from the builtin folder
under-stabilize coupling-relevant hairpins and should not be trusted).
"""

import warnings
from typing import List, Optional

from operon.core import Assembled, OperonHost


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


def to_physical(rate_au: float, host: OperonHost) -> float:
    """Convert a v1_style_rate (au) to the "physical" rate Tian & Salis's
    ``f_unfold = min(1, C * r_phys)`` formula expects.

    The spec is explicit that this conversion is *absorbed into C* when
    using Tian's published C=0.81: "the conversion between RBS-calculator
    au and s^-1 is absorbed into C. If you reuse 0.81 you must use their
    au." I.e. with a properly re-fit ``host.c_unfold``, ``rate_au`` is
    already the right input -- there is no independently-published
    au -> s^-1 formula to implement here, and inventing one would
    silently double-convert. This function exists only as the single
    seam the spec asks for ("put the conversion in one function so
    HTISC, coupling, and Design never mix units"); it is the identity
    map until a real reporter-library refit of ``c_unfold`` says
    otherwise (spec section 7.2, 23).
    """
    return rate_au


class LeaderlessCouplingError(Exception):
    """Raised when a downstream CDS's start is leaderless and coupling
    cannot compute a de novo rate for it (spec section 7.2: this is a
    real architecture -- SD-less overlapping starts -- not an error to
    paper over; callers that hit this should route around ``coupled_tirs``
    for that junction, not force a numeric TIR)."""


def coupled_tirs(
    assembled: Assembled,
    host: OperonHost,
    calc,
    require_vienna: bool = True,
) -> List[float]:
    """Nested 5'->3' iteration computing each CDS's coupled TIR
    (``v1_style_rate`` scale throughout -- spec section 4.5) as
    ``r_reinit + r_denovo`` (section 7.3).

    ``assembled.starts``/``assembled.cds_end`` must be populated (see
    ``operon.assembly.assemble``); ``assembled.junctions[i-1].d_nt`` is
    the intergenic distance between CDS_{i-1} and CDS_i.

    Each step uses the already-coupled upstream rate, not a mono-cistronic
    prediction -- this is why changing an upstream RBS can move a
    downstream CDS's coupled TIR even though its own sequence didn't
    change (known failure mode 1's general form).
    """
    folder_name = getattr(getattr(calc, "folder", None), "name", None)
    if folder_name != "vienna":
        message = (
            "operon.coupling.coupled_tirs: calc.folder is not Vienna-backed "
            f"(got {folder_name!r}). The builtin folder cannot represent "
            "bulge-containing intergenic hairpins and will under-stabilize "
            "coupling-relevant structure -- install ViennaRNA (RNAfold on "
            "PATH) before trusting these numbers."
        )
        if require_vienna:
            raise RuntimeError(message + " Pass require_vienna=False to proceed anyway (not recommended).")
        warnings.warn(message, stacklevel=2)

    n = len(assembled.starts)
    if n == 0:
        return []

    r: List[Optional[float]] = [None] * n
    first = calc.predict_one(assembled.mrna, assembled.starts[0])
    if first.leaderless:
        raise LeaderlessCouplingError(
            f"CDS 0 (start {assembled.starts[0]}) is leaderless; coupling needs a scored first CDS."
        )
    r[0] = first.v1_style_rate

    for i in range(1, n):
        d = assembled.junctions[i - 1].d_nt
        upstream_cds_span = range(assembled.starts[i - 1], assembled.cds_end[i - 1])

        tir_folded = calc.predict_one(assembled.mrna, assembled.starts[i])
        r_reinit = host.k_p * k_reinitiation(d) * r[i - 1]

        if tir_folded.leaderless:
            # No SD at all for this start: a real SD-less overlapping
            # architecture. De novo initiation is 0; keep only re-initiation.
            r[i] = r_reinit
            continue

        tir_unfolded = calc.predict_one(
            assembled.mrna, assembled.starts[i], extra_unpaired_global=upstream_cds_span
        )
        f_unfold = min(1.0, host.c_unfold * to_physical(r[i - 1], host))
        r_denovo = (1 - f_unfold) * tir_folded.v1_style_rate + f_unfold * tir_unfolded.v1_style_rate
        r[i] = r_reinit + r_denovo

    return r
