"""
Inverse RBS design (spec section 5, Appendix A): a v1.0-style simulated
annealer calling ``rbsforge.predict`` as the energy oracle. Mono-
cistronic -- see ``operon/design/CLAUDE.md`` for why its returned TIR is
only a proposal for any CDS after the first.
"""

import math
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

from rbsforge.constants import RT_EFF_DEFAULT, V1_LOG_K
from rbsforge.rbs.calculator import RBSCalculator
from rbsforge.rbs.start_codons import DEFAULT_SCANNED_CODONS

# v1.0 constants (Salis 2009/2011; RBS_MC_Design.py), copied verbatim -- do not retune
TOL_KCAL_MOL = 0.25
RT_INIT = 0.6
ANNEAL_EVERY_N_MOVES = 50
ANNEAL_ACCEPT_LOW = 0.01
ANNEAL_ACCEPT_HIGH = 0.20
MAX_INIT_ENERGY = 10.0
DG_RANGE_HIGH = 25.0
DG_RANGE_LOW = -18.0
MOVE_WEIGHTS = {"replace": 0.80, "insert": 0.10, "delete": 0.10}
MIN_CORE_LEN = 4
MAX_CORE_LEN = 9

IUPAC_ALPHABET: Dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "U",
    "N": "ACGU", "R": "AG", "Y": "CU", "S": "GC", "W": "AU",
    "K": "GU", "M": "AC", "B": "CGU", "D": "AGU", "H": "ACU", "V": "ACG",
}

_COMPLEMENT = {"A": "U", "U": "A", "T": "A", "G": "C", "C": "G"}


@dataclass
class DesignResult:
    rbs: str
    predicted_tir: Optional[float]
    dg_total: Optional[float]
    breakdown: Optional[object]
    n_evals: int
    hit_tol: bool
    objective: float = field(default=float("inf"))


def dg_target_from_tir(target_tir: float, beta: float) -> float:
    """``dG_target = RT_eff * (logK - ln(TIR_target))`` (spec section
    5.1), using RBSForge's own v1.0 K=2500/logK=7.824 constants and the
    host's own beta (RT_eff = 1/beta) rather than a hardcoded 2.222, so
    this works for any HostPack, not just the calibrated *E. coli* one.
    """
    rt_eff = 1.0 / beta if beta else RT_EFF_DEFAULT
    return rt_eff * (V1_LOG_K - math.log(target_tir))


def _has_start_codon(rbs: str, scanned_codons=DEFAULT_SCANNED_CODONS) -> bool:
    seq = rbs.upper().replace("T", "U")
    codons = {c.upper().replace("T", "U") for c in scanned_codons}
    return any(seq[i : i + 3] in codons for i in range(len(seq) - 2))


def _random_base(code: str, rng) -> str:
    options = IUPAC_ALPHABET.get(code.upper(), "ACGU")
    return options[int(rng.random() * len(options))]


def _matches_code(base: str, code: str) -> bool:
    return base.upper() in IUPAC_ALPHABET.get(code.upper(), "ACGU")


def _plant_sd_core(length: int, asd: str, s_opt: int, dg_target: float, rng) -> str:
    """SD-biased random init (spec section 5.4): plant a core
    complementary to the host's anti-SD tail, length biased by how
    strong ``dg_target`` is, at aligned spacing ``s_opt`` from the 3' end
    (where the start codon will follow); randomize the rest.
    """
    frac_strength = min(1.0, max(0.0, (DG_RANGE_HIGH - dg_target) / (DG_RANGE_HIGH - DG_RANGE_LOW)))
    core_len = min(len(asd), max(MIN_CORE_LEN, min(MAX_CORE_LEN, round(MIN_CORE_LEN + frac_strength * (MAX_CORE_LEN - MIN_CORE_LEN)))))

    core = "".join(_COMPLEMENT[b] for b in reversed(asd[-core_len:]))
    core_start = max(0, length - s_opt - core_len)
    core_end = min(length, core_start + len(core))
    core = core[: core_end - core_start]

    bases = [_random_base("N", rng) for _ in range(length)]
    for k, b in enumerate(core):
        bases[core_start + k] = b
    return "".join(bases)


def _strip_start_codons(bases, rng, mask=None):
    """Replace any AUG/GUG/UUG formed within the RBS with a random
    (mask-respecting) base at its middle position, repeating until none
    remain or attempts are exhausted."""
    for _ in range(10):
        seq = "".join(bases)
        if not _has_start_codon(seq):
            return bases
        for i in range(len(seq) - 2):
            if seq[i : i + 3] in {"AUG", "GUG", "UUG"}:
                code = mask[i + 1] if mask else "N"
                bases[i + 1] = _random_base(code, rng)
                break
    return bases


def _score(calc: RBSCalculator, constant_upstream: str, rbs: str, cds_nt: str):
    mrna = constant_upstream + rbs + cds_nt
    start = len(constant_upstream) + len(rbs)
    result = calc.predict_one(mrna, start)
    return result


def design_rbs(
    cds_nt: str,
    target_tir: float,
    hostpack,
    constraint_iupac: Optional[str] = None,
    constant_upstream: Optional[str] = None,
    rbs_init: Optional[str] = None,
    rbs_len_range: Tuple[int, int] = (20, 35),
    max_iter: int = 10000,
    temperature_c: Optional[float] = None,
    rng=None,
) -> DesignResult:
    """Simulated-annealing inverse RBS design against ``rbsforge.predict``
    (Appendix A). See module docstring and ``operon/design/CLAUDE.md``.
    """
    import random as _random_module

    rng = rng or _random_module.Random()
    constant_upstream = constant_upstream or ""
    calc = RBSCalculator(hostpack, temperature_c=temperature_c)
    max_len = min(rbs_len_range[1], hostpack.rbs_max_len)
    min_len = max(1, rbs_len_range[0])

    dg_target = dg_target_from_tir(target_tir, hostpack.beta)

    def objective(rbs: str) -> Tuple[float, Optional[object]]:
        result = _score(calc, constant_upstream, rbs, cds_nt)
        if result.leaderless:
            return float("inf"), result
        return abs(result.delta_g_total - dg_target), result

    n_evals = 0
    mask = constraint_iupac
    fixed_length = mask is not None

    if rbs_init is not None:
        current_str = rbs_init
        current_o, current_result = objective(current_str)
        n_evals += 1
    else:
        length = len(mask) if fixed_length else max_len
        current_str, current_o, current_result = None, None, None
        for _ in range(50):
            candidate = list(_plant_sd_core(length, hostpack.asd, hostpack.s_opt, dg_target, rng))
            if mask:
                candidate = [b if _matches_code(b, mask[i]) else _random_base(mask[i], rng) for i, b in enumerate(candidate)]
            candidate = _strip_start_codons(candidate, rng, mask)
            candidate_str = "".join(candidate)
            o, result = objective(candidate_str)
            n_evals += 1
            current_str, current_o, current_result = candidate_str, o, result  # keep the latest draw regardless
            if o <= MAX_INIT_ENERGY:
                break

    best_str, best_o, best_result = current_str, current_o, current_result

    rt = RT_INIT
    accepted_in_block = 0

    for move_i in range(1, max_iter + 1):
        candidate = list(current_str)
        if fixed_length:
            move = "replace"
        else:
            r = rng.random()
            move = "replace" if r < MOVE_WEIGHTS["replace"] else ("insert" if r < MOVE_WEIGHTS["replace"] + MOVE_WEIGHTS["insert"] else "delete")

        if move == "replace" and candidate:
            pos = int(rng.random() * len(candidate))
            code = mask[pos] if mask else "N"
            candidate[pos] = _random_base(code, rng)
        elif move == "insert" and len(candidate) < max_len:
            pos = int(rng.random() * (len(candidate) + 1))
            candidate.insert(pos, _random_base("N", rng))
        elif move == "delete" and len(candidate) > min_len:
            pos = int(rng.random() * len(candidate))
            candidate.pop(pos)
        else:
            continue  # move not applicable this iteration; try the next one

        candidate = _strip_start_codons(candidate, rng, mask if fixed_length else None)
        candidate_str = "".join(candidate)

        if mask and (len(candidate_str) != len(mask) or not all(_matches_code(b, c) for b, c in zip(candidate_str, mask))):
            continue  # hard reject: IUPAC violated
        if not (min_len <= len(candidate_str) <= max_len):
            continue  # hard reject: length out of range

        candidate_o, candidate_result = objective(candidate_str)
        n_evals += 1
        if candidate_result.leaderless:
            continue  # hard reject: intended start must not be leaderless

        accept = candidate_o < current_o or rng.random() < math.exp((current_o - candidate_o) / rt)
        if accept:
            current_str, current_o, current_result = candidate_str, candidate_o, candidate_result
            accepted_in_block += 1
            if candidate_o < best_o:
                best_str, best_o, best_result = candidate_str, candidate_o, candidate_result

        if best_o <= TOL_KCAL_MOL:
            return DesignResult(
                rbs=best_str, predicted_tir=best_result.v1_style_rate,
                dg_total=best_result.delta_g_total, breakdown=best_result.breakdown,
                n_evals=n_evals, hit_tol=True, objective=best_o,
            )

        if move_i % ANNEAL_EVERY_N_MOVES == 0:
            accept_ratio = accepted_in_block / ANNEAL_EVERY_N_MOVES
            if accept_ratio > ANNEAL_ACCEPT_HIGH:
                rt /= 2
            elif accept_ratio < ANNEAL_ACCEPT_LOW:
                rt *= 2
            accepted_in_block = 0

    return DesignResult(
        rbs=best_str, predicted_tir=best_result.v1_style_rate,
        dg_total=best_result.delta_g_total, breakdown=best_result.breakdown,
        n_evals=n_evals, hit_tol=False, objective=best_o,
    )
