# operon.design

Two related but separable pieces: inverse RBS design (`design_rbs`, spec
section 5, Appendix A) and the multi-objective operon-level search
(`design`, spec section 18). **Status: both implemented.** This was the
last module in the spec's implementation order (section 20, steps 6 and
14) for good reason — everything else was built first.

## `design_rbs` — implemented in `design_rbs.py`

Verified across the full realistic TIR range (10 to 50,000 au) in
`tests/test_design_rbs.py`: every case there converges (`hit_tol=True`,
objective <= 0.25 kcal/mol) within a few thousand evaluations at most,
usually far fewer. `dg_target_from_tir` uses RBSForge's own
`V1_LOG_K`/`RT_EFF_DEFAULT` constants and the *host's own* `beta`
(`RT_eff = 1/beta`) rather than hardcoding 2.222, so it isn't silently
*E. coli*-only the way a literal `2.222` in the formula would be.

A v1.0-style simulated annealer (Salis, Mirsky & Voigt 2009,
https://doi.org/10.1038/nbt.1568; Salis 2011 section 4.4 for IUPAC) that
calls `rbsforge.predict` as a black-box energy oracle. No dependency on
coupling, footprint unfolding, or Vienna.

```
dG_target = 2.222 * (7.824 - ln(target_tir))
O(seq)    = |predict(constant_upstream + rbs + cds_nt).dG_total - dG_target|
```

Constants to copy exactly from v1.0 (spec section 5.3, Appendix A) — do
not retune any of these without a reason from the spec: `tol = 0.25`
kcal/mol, `RT_init = 0.6`, anneal every 50 moves (`RT /= 2` if accept
ratio `> 0.20`, `RT *= 2` if `< 0.01`), `max_iter = 10000` default (cap
lower, ~2000, for the inner-loop use from `design()` below, since it will
call this many times), move weights `{replace: 0.80, insert: 0.10, delete:
0.10}` (replace-only if the IUPAC mask is fixed-length — Operon
constraints usually are), `max_init_energy = 10.0` kcal/mol (redraw if
worse), `dG_range_high = 25.0` / `dG_range_low = -18.0` for the SD-biased
init, `rbs_max_len = 35` (from `HostPack.rbs_max_len`).

**Reject, don't Metropolis, on**: IUPAC violation, a start codon
appearing inside the RBS, the intended start scoring leaderless, length
outside `rbs_len_range`/`rbs_max_len`.

`DesignResult.predicted_tir` is `v1_style_rate`, not
`translation_initiation_rate` — `dG_target` was derived by inverting the
`v1_style_rate` formula (`K * exp(-dG/RT_eff)`), so that's the scale that
actually matches `target_tir` apples-to-apples. Don't switch it to the
proportional scale without also changing `dg_target_from_tir`.

**Explicitly skip** (spec section 5.3's table): `kinetic_score`/
`three_state_indicator` guards (RBSForge doesn't compute the folding-
kinetics quantities they need — don't fake them) and helical-loop repair
(Vienna's `dG_mRNA` already penalizes pathological folds when Vienna is
the folder).

**This function is mono-cistronic — it does not know about coupling.**
For every CDS after the first, its returned TIR is a proposal, not
ground truth: `design()` must re-score against the assembled mRNA via
`operon.coupling` afterward. Never surface `design_rbs`'s TIR for CDS_2+
as if it were the final coupled TIR (spec section 5.6, known failure
mode 2). Pass `constant_upstream` = the actual upstream CDS tail so the
folded window at least sees the junction, but treat the result as a
starting point for Stage B, not an answer.

If `design_rbs` isn't built yet and you need a placeholder proposal
inside `design()`'s inner loop, spec section 5.8 gives one: plant
`AGGAGG` at `s_opt`, randomize the rest of a 35-mer inside the IUPAC
mask, and let the outer GA move it from there.

## `design` (NSGA-II) — implemented in `nsga2.py`; ready-made objectives in `objectives.py`

The engine (`design()`) is objective-agnostic — it takes
`objectives: List[Callable[[Operon, Assembled], float]]` (minimize) and
`constraints: List[Callable[[Operon, Assembled], bool]]` (True =
satisfied) and runs standard NSGA-II (fast non-dominated sort +
crowding-distance selection + constrained dominance: feasible always
beats infeasible, fewer constraint violations wins among infeasible
individuals) over them. `objectives.py` has ready-made ones built from
the now-implemented modules — `tir_error_objective`,
`repeat_length_objective`, `htisc_count_objective`,
`synthesis_feasibility_constraint`, `protein_sequence_constraint` — so
`design()` is demonstrably usable end to end against real scanners, not
just an abstract shell (see `tests/test_design_nsga2.py`).

**A seed `Operon` that already violates a hard constraint you're
searching under poisons the whole run.** Elitism (`combined = population
+ offspring`, then non-dominated sort) means a feasible individual, once
found, persists across generations — but if the *seed* and every early
mutation of it are infeasible (e.g. a homopolymer run already in the
seed RBS), the search may never find a feasible individual in a short
run at all, and `design()` will then return infeasible designs in "front
0" (they're still non-dominated *among infeasible individuals*). This
bit an early version of this module's own test suite — the seed's RBS
had a 15-nt run of A's, silently failing `synthesis_feasibility_constraint`
from generation 0. Start from a feasible seed.

`tir_error_objective` is mono-cistronic (`calc.predict_one` per CDS, not
`operon.coupling`) by default — cheap and Vienna-independent, which
matters inside a GA loop that may call it hundreds of times. A caller
who needs a coupling-aware TIR objective should build one around
`operon.coupling.coupled_tirs` themselves; wiring that in as the default
would make every `design()` call require Vienna.

## `design` — the outer Pareto search, general design notes

**NSGA-II (or another Pareto GA), never a weighted-sum scalarization.**
The spec is explicit about why (section 18.2, known failure mode 9): a
weighted sum hides exactly the tradeoffs — TIR error vs. repeat count vs.
HTISC — that the real Operon Calculator's web tool surfaces by returning
several "equally optimal" designs. If you're tempted to combine
objectives into one scalar for convenience, don't; keep the objective
vector (section 18.2, `f1`...`f12`) and search over it directly.

Hard constraints get repaired or the candidate is rejected — they never
sit on the Pareto front: `n_forbidden_re == 0`, protein sequence
unchanged, IUPAC RBS constraints satisfied, homopolymer limits from
`operon.synthesis`.

Two-stage practical scheme (section 18.3): **Stage A** — decoupled
`design_rbs` proposal per CDS plus recoding, choosing junctions per
`operon.assembly`'s policy. **Stage B** — assemble, run full Predict
including `operon.coupling`, and NSGA-II-polish whatever drifted. Stage A
alone is already a usable, much cheaper "Design mode"; don't skip
straight to Stage B as the only mode.

Mutation operators stay **local** (section 18.4): RBS nucleotide flip
inside the IUPAC mask, single codon swap to a synonym (repaired for
RE/12-mer/homopolymer via `operon.repeats`/`operon.synthesis`),
N-terminal-recode-plus-RBS-redesign as one paired move, junction recode
targeted at a specific `d`/`ΔG_coupling`, and *targeted* repair (a
cryptic-promoter hit at position `p` mutates only that hexamer, not the
whole CDS). Crossover is uniform **over cistrons**, never over raw
nucleotides — nucleotide-level crossover wrecks Shine-Dalgarno sequences.

Return a de-duplicated (by 12-mer identity), non-dominated set of 5-20
designs — not a single winner.

## Interface to build toward

```python
def design_rbs(cds_nt, target_tir, hostpack, constraint_iupac=None,
                constant_upstream=None, rbs_init=None, rbs_len_range=(20, 35),
                max_iter=10000, temperature_c=None): ...
def design(operon, calc, objectives, constraints) -> list: ...
```
