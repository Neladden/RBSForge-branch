# operon.synthesis — Module J

DNA synthesis complexity and restriction-site constraints. **Status:
prepared subsection, not implemented.** Spec section 15. Source: Halper,
Hossain & Salis 2020 (https://doi.org/10.1021/acssynbio.9b00460,
Synthesis Success Calculator).

## Start with the hard rules, not the trained forest

The section 15.1 hard-rule table (homopolymer C/G <=9, homopolymer A/T
<=13, no 6x period-3 trinucleotide tandem, no 10x dinucleotide tandem,
windowed GC% bounds at 20/100 bp and the terminal 30 bp, local repeat-
density caps) already defines a feasible set on its own. The SSC random
forest (F1 0.928) is a refinement for designing >1 kb fragments — build
it later, as a smooth `P(synthesis success)` objective layered on top,
not as a prerequisite.

**Use these as GA constraints, not just scores.** A candidate that
violates a hard rule should be rejected or repaired by `operon.design`'s
search, the same way a forbidden restriction site is (below) — don't
fold them into a soft objective that a Pareto front could "trade off"
against TIR.

## Restriction sites

User-supplied enzyme list. Exclude recognition sequences on both
strands (both orientations if the enzyme is non-palindromic). This is a
**hard constraint**: `n_forbidden_re > 0` is invalid, full stop — not a
score component. If a user wants an RE site at a part boundary on
purpose, that's expressed by making those nucleotides constant in the
IUPAC RBS/intergenic constraint fed to `operon.design`, not by this
module's logic.

## This module overlaps `operon.repeats` — don't duplicate work

Halper et al. found repeats were the SSC forest's most important feature
class, and killing 12-mer repeats (`operon.repeats`) does most of the
synthesis-complexity work already. If you're about to reimplement
repeat-density scanning here, check `operon.repeats` first — call into
it instead.

## One thing to leave alone

RBSForge's `delta_G_stacking = 0.0` means polyA/polyU spacers (which this
module's homopolymer rules will often push a design toward, especially
combined with mRNA-stability or RNase-avoidance pressure from
`operon.stability`) are mis-ranked on TIR by the translation engine.
That's a known, accepted gap — treat TIR as ordinal rather than absolute
when homopolymer spacers dominate a Pareto set; don't try to compensate
for it inside this module.

## Interface to build toward

```python
def synthesis_score(dna: str): ...
```
