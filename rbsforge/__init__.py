"""RBSForge: a from-scratch, dependency-free, organism-configurable
reimplementation of the Salis Lab RBS Calculator's free-energy model for
predicting bacterial translation initiation rates -- including for
organisms outside its E. coli-tuned defaults (see docs/MODEL.md).

Quick start::

    from rbsforge import predict

    result = predict(
        "AGGAGGACAACTAAATGAAACGCATTAGCACCACC...",
        species="ecoli",
        temperature_c=37.0,
    )
    best = result.best
    print(best.codon, best.delta_g_total, best.translation_initiation_rate)
    print(best.breakdown)
"""
from typing import Optional

from .hostpack import BUILTIN_HOSTPACKS, HostPack, get_hostpack
from .rbs.calculator import PredictResult, RBSCalculator, StartCodonResult

__all__ = [
    "predict",
    "RBSCalculator",
    "PredictResult",
    "StartCodonResult",
    "HostPack",
    "get_hostpack",
    "BUILTIN_HOSTPACKS",
]


def predict(
    mrna_sequence: str,
    species: str = "ecoli",
    temperature_c: Optional[float] = None,
) -> PredictResult:
    """Predict translation initiation rates for every start codon in
    `mrna_sequence`.

    Args:
        mrna_sequence: the mRNA (or DNA, T is accepted) sequence, ideally
            from the transcription start site through some of the CDS.
        species: a built-in HostPack name (see `BUILTIN_HOSTPACKS`) or
            pass a `HostPack` instance directly to `RBSCalculator` for a
            custom/uncharacterized organism.
        temperature_c: expression/growth temperature in Celsius. Defaults
            to the HostPack's own `t_growth_c` if omitted.
    """
    calculator = RBSCalculator.for_species(species, temperature_c=temperature_c)
    return calculator.predict(mrna_sequence)
