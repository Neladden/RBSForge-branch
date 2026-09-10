# operon.elongation — Module C

Translation elongation rate (TER) prediction and synonymous codon
recoding. **Status: implemented, with a placeholder codon-usage table.**
Spec section 8. Independent of coupling; operates on one CDS at a time.

## What's real vs. placeholder here

- `genetic_code.py`'s `GENETIC_CODE`/`translate` are the universal
  genetic code — 100% factual, not organism-specific, safe to use as-is.
- `codon_tables.py`'s `ECOLI_HIGHLY_TRANSLATED`/`ECOLI_GENOMIC` are
  **not** a measured codon-usage table. Each amino acid's single most
  commonly cited *E. coli* preferred codon gets weight 1.0; every other
  synonym gets a flat placeholder weight (0.3) — a two-tier
  approximation, not a real relative-frequency ranking. Both "high" and
  "genomic" tables are currently identical (built the same way) because
  there is no real whole-genome usage data wired in yet. **Before this
  module's numbers should be trusted for anything beyond "does the
  pipeline run," replace `codon_tables.py` with real data** (Codon Usage
  Database, https://www.kazusa.or.jp/codon/, or a table computed from the
  actual expression host's own genes) — see the module's own docstring
  for where the real `table_high`/`table_genomic` split (spec section
  8.2) is supposed to come from.
- `recode_cds(aa, table, rng)` matches the spec section 25 signature
  exactly: a full initial protein->RNA recode, sampling each position
  independently from `table[aa]`. `mutate_codon(cds_nt, table, rng)` is
  the separate single-position GA neighborhood operator spec section 8.3
  describes — don't conflate the two; `recode_cds` is not meant to be
  called repeatedly as a mutator.
- `ter_nt_s`'s `r_max` divides `host.r_elong_nominal_nt_s` (nt/s) by 3 to
  get codons/s before applying the harmonic-mean formula — the spec's
  own "~20 aa/s ~= 60 nt/s" note is the tell that `r_max` in
  `r_elong_codon[c] = r_max * w[c]/w_max` must be in codons/s, not nt/s;
  get this wrong and TER comes out ~3x too high.

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
