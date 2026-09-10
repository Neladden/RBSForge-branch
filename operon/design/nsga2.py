"""
NSGA-II Pareto search over RBS nucleotides and synonymous codon choices
(spec section 18). Objective-agnostic: the caller supplies the objective
vector and hard constraints (spec section 18.2's `f1..f12` and its hard-
constraint list) as plain functions; ``operon/design/objectives.py`` has
a set of ready-made ones built from the now-implemented modules.

Not a weighted-sum scalarizer -- returns a de-duplicated, non-dominated
set of designs (spec section 18.5), never a single winner.
"""

from dataclasses import dataclass, replace
from typing import Callable, List, Optional

from operon.core import Assembled, Operon
from operon.elongation import ECOLI_HIGHLY_TRANSLATED, mutate_codon
from operon.assembly import assemble

from .design_rbs import _random_base

ObjectiveFn = Callable[[Operon, Assembled], float]  # lower is better, always
ConstraintFn = Callable[[Operon, Assembled], bool]  # True = satisfied


@dataclass
class Design:
    operon: Operon
    assembled: Assembled
    objective_values: List[float]


@dataclass
class _Individual:
    operon: Operon
    assembled: Optional[Assembled] = None
    objective_values: Optional[List[float]] = None
    n_violations: int = 0
    feasible: bool = True
    rank: int = 0
    crowding: float = 0.0

    def genome_key(self):
        return (tuple(c.id + "|" + c.aa_or_nt for c in self.operon.cds_list), tuple(self.operon.rbs_list))


def _dominates(a: List[float], b: List[float]) -> bool:
    return all(x <= y for x, y in zip(a, b)) and any(x < y for x, y in zip(a, b))


def _constrained_dominates(a: _Individual, b: _Individual) -> bool:
    if a.feasible and not b.feasible:
        return True
    if not a.feasible and b.feasible:
        return False
    if not a.feasible and not b.feasible:
        return a.n_violations < b.n_violations
    return _dominates(a.objective_values, b.objective_values)


def _fast_nondominated_sort(pop: List[_Individual]) -> List[List[int]]:
    n = len(pop)
    dominated_by = [set() for _ in range(n)]
    domination_count = [0] * n
    fronts: List[List[int]] = [[]]

    for p in range(n):
        for q in range(n):
            if p == q:
                continue
            if _constrained_dominates(pop[p], pop[q]):
                dominated_by[p].add(q)
            elif _constrained_dominates(pop[q], pop[p]):
                domination_count[p] += 1
        if domination_count[p] == 0:
            pop[p].rank = 0
            fronts[0].append(p)

    i = 0
    while fronts[i]:
        next_front = []
        for p in fronts[i]:
            for q in dominated_by[p]:
                domination_count[q] -= 1
                if domination_count[q] == 0:
                    pop[q].rank = i + 1
                    next_front.append(q)
        i += 1
        fronts.append(next_front)
    return fronts[:-1]  # drop the trailing empty front


def _assign_crowding(pop: List[_Individual], front: List[int]) -> None:
    if not front:
        return
    for i in front:
        pop[i].crowding = 0.0
    if not pop[front[0]].objective_values:
        return
    n_obj = len(pop[front[0]].objective_values)
    for m in range(n_obj):
        front_sorted = sorted(front, key=lambda i: pop[i].objective_values[m])
        lo = pop[front_sorted[0]].objective_values[m]
        hi = pop[front_sorted[-1]].objective_values[m]
        pop[front_sorted[0]].crowding = float("inf")
        pop[front_sorted[-1]].crowding = float("inf")
        spread = hi - lo
        if spread == 0:
            continue
        for k in range(1, len(front_sorted) - 1):
            prev_v = pop[front_sorted[k - 1]].objective_values[m]
            next_v = pop[front_sorted[k + 1]].objective_values[m]
            pop[front_sorted[k]].crowding += (next_v - prev_v) / spread


def _tournament_select(pop: List[_Individual], rng) -> _Individual:
    a, b = pop[int(rng.random() * len(pop))], pop[int(rng.random() * len(pop))]
    if a.rank != b.rank:
        return a if a.rank < b.rank else b
    return a if a.crowding >= b.crowding else b


def _mutate_operon(operon: Operon, rng, codon_table=None) -> Operon:
    """One local mutation (section 18.4): an RBS nucleotide flip inside
    its IUPAC mask, or a protein-preserving codon swap -- picked
    uniformly among every mutable target in the operon."""
    codon_table = codon_table or ECOLI_HIGHLY_TRANSLATED
    targets = []
    for i, rbs in enumerate(operon.rbs_list):
        if rbs:
            targets.append(("rbs", i))
    for i in range(len(operon.cds_list)):
        targets.append(("codon", i))
    if not targets:
        return operon

    kind, i = targets[int(rng.random() * len(targets))]
    new_cds_list = list(operon.cds_list)
    new_rbs_list = list(operon.rbs_list)

    if kind == "rbs":
        rbs = new_rbs_list[i]
        mask = operon.cds_list[i].rbs_constraint if i < len(operon.cds_list) else None
        mask = mask or ("N" * len(rbs))
        pos = int(rng.random() * len(rbs))
        code = mask[pos] if pos < len(mask) else "N"
        new_base = _random_base(code, rng)
        new_rbs_list[i] = rbs[:pos] + new_base + rbs[pos + 1 :]
    else:
        cds = new_cds_list[i]
        try:
            new_nt = mutate_codon(cds.aa_or_nt, codon_table, rng)
        except ValueError:
            new_nt = cds.aa_or_nt  # not a valid in-frame CDS (e.g. constant filler) -- leave untouched
        new_cds_list[i] = replace(cds, aa_or_nt=new_nt)

    return replace(operon, cds_list=new_cds_list, rbs_list=new_rbs_list)


def _evaluate(operon: Operon, objectives: List[ObjectiveFn], constraints: List[ConstraintFn]) -> _Individual:
    ind = _Individual(operon=operon)
    try:
        assembled = assemble(operon)
    except ValueError:
        ind.feasible = False
        ind.n_violations = len(constraints) + 1
        ind.objective_values = [float("inf")] * len(objectives)
        return ind

    ind.assembled = assembled
    violations = sum(0 if c(operon, assembled) else 1 for c in constraints)
    ind.n_violations = violations
    ind.feasible = violations == 0
    ind.objective_values = [f(operon, assembled) for f in objectives]
    return ind


def design(
    operon: Operon,
    calc,
    objectives: List[ObjectiveFn],
    constraints: List[ConstraintFn],
    population_size: int = 30,
    n_generations: int = 20,
    max_results: int = 20,
    rng=None,
) -> List[Design]:
    """NSGA-II Pareto search (section 18.3's inner loop, section 18.4's
    mutation operators). Returns a de-duplicated (by exact RBS/CDS
    genome), non-dominated set of designs -- never a single winner.

    ``calc`` is accepted (matching spec section 25's signature) for
    objective/constraint functions that need it (an `rbsforge.RBSCalculator`);
    the engine itself doesn't call it directly, since scoring is entirely
    delegated to the supplied ``objectives``/``constraints``.
    """
    import random as _random_module

    rng = rng or _random_module.Random()

    population = [_evaluate(operon, objectives, constraints)]
    while len(population) < population_size:
        candidate = _mutate_operon(population[0].operon, rng)
        population.append(_evaluate(candidate, objectives, constraints))

    for _ in range(n_generations):
        offspring = []
        while len(offspring) < population_size:
            parent = _tournament_select(population, rng)
            child_operon = _mutate_operon(parent.operon, rng)
            offspring.append(_evaluate(child_operon, objectives, constraints))

        combined = population + offspring
        fronts = _fast_nondominated_sort(combined)
        next_population: List[_Individual] = []
        for front in fronts:
            _assign_crowding(combined, front)
            if len(next_population) + len(front) <= population_size:
                next_population.extend(combined[i] for i in front)
            else:
                remaining = population_size - len(next_population)
                front_sorted = sorted(front, key=lambda i: combined[i].crowding, reverse=True)
                next_population.extend(combined[i] for i in front_sorted[:remaining])
                break
        population = next_population

    fronts = _fast_nondominated_sort(population)
    pareto = [population[i] for i in fronts[0]] if fronts else []

    seen = set()
    results: List[Design] = []
    for ind in pareto:
        key = ind.genome_key()
        if key in seen:
            continue
        seen.add(key)
        results.append(Design(operon=ind.operon, assembled=ind.assembled, objective_values=ind.objective_values))
        if len(results) >= max_results:
            break
    return results
