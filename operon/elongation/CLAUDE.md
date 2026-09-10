# operon.elongation — Module C

Translation elongation rate (TER) prediction and synonymous codon
recoding. **Status: prepared subsection, not implemented.** Spec
section 8. Independent of coupling; operates on one CDS at a time.

## Before implementing

- **Harmonic mean, not arithmetic**, for `TER_cds` — one slow codon
  dominates dwell time:
  `TER_cds = 3 / mean(1/r_elong_codon[c] for c in cds)`.
- Two codon tables are Design-mode inputs, not two implementations to
  pick one of: **Highly Translated** (always the max-`w` codon; use for a
  small number of overexpressed CDSs) and **Balanced** (mixes toward the
  target TIR via `softmax(mix * log(table_high) + (1-mix) *
  log(table_genomic))`, `mix = clip(log10(target_tir)/log10(tir_max), 0,
  1)`; use when many CDSs are expressed at varied levels). Fast
  initiation + slow elongation on the wrong table causes pile-ups that
  silently invalidate the TIR you designed for.
- Host-specific tables only — never reuse *E. coli* codon usage for
  another organism (`table_high` = usage of highly-expressed genes,
  `table_genomic` = whole-genome usage, both host-specific).
- **`recode_cds` is a neighborhood operator for the Design-mode GA, not a
  one-shot rewrite.** A single greedy CAI pass fights every other
  objective (repeats, cryptic promoters, HTISC, RNase sites) — implement
  it as "pick a random synonym-eligible codon, swap to a synonym biased
  by the active table, repair for RE/12-mer/homopolymer," callable
  repeatedly, not as a whole-CDS transform.
- **Until CDS-footprint unfolding is wired into `rbsforge`** (see
  `rbsforge/CLAUDE.md`), freeze the first ~10-15 codons during recoding —
  they're TIR-sensitive and a naive recode there silently invalidates the
  upstream RBS design (known failure mode 3). Once footprint unfolding
  lands, the correct behavior is to re-run `predict_tir` after every
  N-terminal recode instead of freezing — switch to that, don't keep
  both paths.
- Internal start-codon suppression (removing/weakening an in-CDS
  AUG/GUG/UUG or its SD) is a recoding move that serves `operon.htisc`'s
  objective — implement it here as a mutator, but the HTISC *count* is
  `operon.htisc`'s objective, not folded into this module's TIR error
  (spec section 8.4).

## Interface to build toward

```python
def recode_cds(aa: str, table, rng) -> str: ...
def ter_nt_s(cds_nt: str, host) -> float: ...
```
