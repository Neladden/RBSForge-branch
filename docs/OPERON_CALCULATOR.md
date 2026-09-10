# Operon Calculator — architecture and status

This is the short index. The full design spec every module below
implements against is `OPERON_CALCULATOR_SPEC.md` in this directory —
read that for algorithms, citations, and worked examples; this page is
just the map.

## Why RBSForge is a subset

RBSForge (`rbsforge/`) reimplements the Salis Lab **RBS Calculator**:
given an assembled mRNA, it predicts the translation initiation rate
(TIR) of every candidate start codon. That is exactly one box in the
Operon Calculator's pipeline — the rest of this project is the
assembler that builds the mRNA RBSForge scores, the coupling model that
makes CDS-after-CDS TIR not just a mono-cistronic RBSForge call, and the
independent scanners/recoders (codon optimization, mRNA stability,
cryptic promoters, terminators, repeats, synthesis constraints) that
turn "predict a TIR" into "predict and design a whole operon."

## Pipeline

```
                    +-----------------------------------------+
                    |           OperonAssembler                |
                    |  promoter + [RBS_i + CDS_i] + terminator |
                    +------------------+------------------------+
                                       | assembled DNA / mRNA
          +---------------+------------+------------+-----------------+
          v               v            v            v                 v
   Translation      mRNAStability   CrypticTX    CrypticTerm     Instability
   (RBSForge +      (RNase +        (Promoter     (intrinsic +    (repeats,
    coupling +       ribosome        Calculator    rho + pause)    IS, att, RE,
    TIR + TER +      protection)     both strands)                 synthesis)
    HTISC)
          |               |            |            |                 |
          +---------------+------------+------------+-----------------+
                                       |
                                       v
                            ObjectiveVector / Scores
                                       |
                    Design mode only:  MultiObjectiveSearch
                    mutates RBS IUPAC windows and synonymous
                    codon choices, re-assembles, re-scores
```

Every box is its own function with a closed interface (spec section 25).
They share a host-organism record (`operon.core.OperonHost`) and a
sequence coordinate system, and nothing else — swapping one module (e.g.
a different Promoter Calculator for a non-sigma70 host) never requires
touching another.

## Module status

| Module | Subpackage | What | Status |
|---|---|---|---|
| — | `rbsforge` (sibling package) | TIR engine, Predict mode, incl. `predict_one` + CDS-footprint unfolding | **implemented** |
| — | `operon.core` | shared types: `OperonHost`, `CDS`, `Operon`, `Assembled`, `Junction` | **implemented** |
| A | `operon.assembly` | operon assembly, intergenic distance `d` | **implemented** |
| B | `operon.coupling` | translational coupling (Tian & Salis 2015) | **implemented** |
| C | `operon.elongation` | TER + synonymous codon recoding | **implemented** (placeholder codon-usage table) |
| D | `operon.stability` | mRNA stability (Cetnar & Salis 2021/2024) | **features implemented**; final rate needs real coefficients |
| E | `operon.promoter_calculator` | sigma70 promoter / cryptic-promoter scan | **implemented** |
| F | `operon.htisc` | highly translated internal start codons | **implemented** |
| G | `operon.pauses` | ribosomal pause sites | **implemented** (3 of 4 signals) |
| H | `operon.terminators` | intrinsic + rho-dependent terminators | **implemented** (presence layer) |
| I | `operon.repeats` | repeats (>=12 bp), IS/att sites | **implemented** (seed-length, not extended) |
| J | `operon.synthesis` | synthesis complexity, restriction sites | **implemented** (numeric hard rules; no hairpin/G-quad) |
| K | — (bookkeeping over other modules) | system-level RNAP/ribosome load | not started |
| — | `operon.design` | `design_rbs` (inverse RBS) + NSGA-II operon design | **implemented** |

Every module above has real tests in `tests/` — see each subpackage's
own `CLAUDE.md` for exactly what's implemented vs. deliberately deferred
(most modules have at least one documented scope boundary; "implemented"
here means "has a real, tested algorithm behind it," not "matches the
full spec with nothing left to extend").

## Implementation order (history)

The spec (section 20) gave a dependency-ordered build sequence; this is
the order it actually happened in, since it diverged from the spec's own
suggested order in two places:

1. Assembler + intergenic distance (Module A)
2. TIR scale freeze + `predict_tir` adapter — already true of `rbsforge`
   from the start
3. **Promoter Calculator (Module E)** — pulled forward out of spec order
   (originally step 10) because it has no dependency on any other
   unimplemented module; it only needs a DNA sequence
4. HTISC (Module F) — cheapest step; needs no new physics, only wraps
   `rbsforge.RBSCalculator.predict`
5. `predict_one` + forced-unpaired + CDS-footprint unfolding — an
   `rbsforge` change, the prerequisite for Module B
6. Coupling (Module B) — built against the builtin folder with a
   `require_vienna` refuse/warn gate (this environment has no `RNAfold`
   on `PATH`; the tests exercise the warn path explicitly)
7. Codon recoding + TER (Module C)
8. Repeat finder (Module I)
9. RE sites + homopolymer/GC (Module J)
10. Intrinsic + rho-dependent terminator scan (Module H)
11. Pause / internal-SD scan (Module G)
12. mRNA stability features (Module D) — final rate deliberately left
    gated on real coefficients, not shipped with an invented default
13. `design_rbs` (inverse RBS design) — verified across the full 10-
    50,000 au target-TIR range
14. NSGA-II design search (`operon.design.design`, plus ready-made
    objectives in `operon.design.objectives`)

Remaining, later-layer expansions per spec section 20's step 15 (2024
GBDT stability, trained RhoTermPredict/OPLS-DA classifier, IS tables,
system load bookkeeping, v2.0 standby, SSC's random forest) are not
implemented — see each module's own `CLAUDE.md` for what's deferred and
why.

## What not to do

The spec's "known failure modes" (section 24) and "what to copy vs.
re-fit" (section 23) apply across every module built here. The two that
matter most while modules are still being filled in:

- Don't put OperonHost-only fields (codon tables, IS motifs, coupling
  constants) into RBSForge's `HostPack`. Wrap it (`operon.core.OperonHost`
  does this); don't fork it.
- Don't copy a paper's fitted constant across a boundary it wasn't fit
  for: Tian 2015's `C` and `k_P` are on their own RBS-calculator au scale
  and need re-fitting against RBSForge's `v1_style_rate`; the Promoter
  Calculator's 343 coefficients, by contrast, are reused as-is here
  because they're already a from-scratch trained model over raw DNA
  sequence, not something layered on an au scale that changed underneath
  them.
