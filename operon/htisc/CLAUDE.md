# operon.htisc — Module F

Highly translated internal start codons: RBSForge run as a scanner, plus
a threshold. **Status: implemented.** Spec section 11.

## Why this one is easy

Of everything in the spec's "what already works without new RBS physics"
list (section 4.6), this needs *nothing* the repo doesn't already have.
The whole implementation is:

```python
def scan_htisc(mrna, annotated_starts, calc, tau_htisc=1000.0):
    hits = []
    out = calc.predict(mrna)
    for r in out.results:
        if r.leaderless or r.index in annotated_starts:
            continue
        if r.v1_style_rate >= tau_htisc:
            hits.append(r)
    return hits
```

against `rbsforge.RBSCalculator.predict` as it exists today — no
`predict_one`, no coupling, no Vienna requirement. If you're looking for
something to implement first in `operon.*` beyond the promoter
calculator, this is it.

## Gotchas

- `tau_htisc` defaults to `1000.0` on the `v1_style_rate` scale (already
  defined as `DEFAULT_TAU_HTISC` in `__init__.py`) — **not** on
  `proportional_rate`. Mixing scales here is exactly known failure mode
  10; RBSForge's own TIR-scale freeze (spec section 4.5) exists partly so
  this threshold means one fixed thing.
- Leaderless starts have `tir=None` and are never HTISC hits — skip them,
  don't coerce to `0`.
- `annotated_starts` must be checked by `r.index`, not by codon identity
  — an internal AUG in a different frame is a different `index` even if
  it's also "AUG".
- This module only *counts* hits; the Design-mode kill (synonym swaps
  that remove the start or its SD, or introduce a nearby stop) is
  `operon.elongation`'s job (spec section 8.4) — don't implement mutation
  logic here.

## Interface to build toward

```python
def scan_htisc(mrna: str, annotated_starts, calc, tau_htisc: float = DEFAULT_TAU_HTISC) -> list: ...
```
