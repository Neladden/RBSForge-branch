# operon.synthesis — Module J

DNA synthesis complexity and restriction-site constraints. **Status:
implemented** for the section 15.1 numeric hard-rule table and the RE-
site hard constraint. Source: Halper, Hossain & Salis 2020
(https://doi.org/10.1021/acssynbio.9b00460, Synthesis Success
Calculator).

## Not implemented — deliberately out of scope so far

- Hairpin/GC-rich/terminal-hairpin detection, G-quadruplex, and i-motif
  motifs (section 15.1's last few rows) are not implemented. They need
  either folding (RBSForge's folder) or dedicated motif detectors this
  pass didn't build — don't assume `synthesis_score` catches them.
- The trained random-forest `P(synthesis success)` refinement (SSC's own
  headline result, F1 0.928) is not implemented; see below.
- `_repeat_density_hits` reuses `operon.repeats.find_repeats`, whose
  direct/inverted hits are seed(k)-length, not extended (see
  `operon/repeats/CLAUDE.md`) — density estimates here are therefore
  conservative (a real, longer repeat is undercounted), not exact.

## A gotcha already fixed once — don't reintroduce it

A run of a single repeated base trivially satisfies *any* period's
tandem recurrence (`seq[j] == seq[j-period]` holds for constant regions
regardless of `period`). `_tandem_hits` explicitly skips a period-2/3
"unit" that's actually just one repeated character (`len(set(unit)) ==
1`) so a long homopolymer isn't double-reported as a dinucleotide *and*
trinucleotide tandem violation on top of its own homopolymer hit. If you
touch `_tandem_hits`, keep that guard.

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
