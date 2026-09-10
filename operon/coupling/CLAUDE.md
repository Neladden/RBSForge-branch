# operon.coupling — Module B

Translational coupling: the only translation physics not already in
RBSForge. **Status: prepared subsection, not implemented — and currently
blocked.** Spec section 7, source Tian & Salis, NAR 2015
(https://doi.org/10.1093/nar/gkv635).

## Blocked on two prerequisites — check both before starting

1. **`rbsforge.RBSCalculator` must expose `predict_one(mrna, start,
   extra_unpaired_global=None)`.** It does not yet (see
   `rbsforge/CLAUDE.md`). Coupling's de novo term needs two calls to it
   per downstream CDS — one folded, one with the upstream CDS's
   nucleotides forced unpaired — and both must go through the *same*
   constrained-fold code path CDS-footprint unfolding uses (spec section
   6.5). Do not build a second, coupling-only constrained-fold mechanism.
2. **Vienna (`RNAfold` on `PATH`) is required, not optional, when this
   module is on.** The builtin folder cannot represent bulge-containing
   intergenic hairpins; running coupling on it under-stabilizes
   inhibitory helices and flattens the sigmoid (known failure mode 11).
   Refuse or warn if `folder.backend != "vienna"` rather than silently
   producing a wrong answer.

## The model, once unblocked

```
r_1 = predict_one(mrna, start_1).tir                      # no coupling
r_reinit_i = host.k_p * k_reinitiation(d) * r_{i-1}        # k_reinitiation already implemented, see below
f = min(1, host.c_unfold * to_physical(r_{i-1}, host))
r_denovo_i = (1-f) * r_folded + f * r_unfolded
r_i = r_reinit_i + r_denovo_i
```

- `k_reinitiation(d)` **is already implemented** in `__init__.py` — it's
  a ratio and transfers across au scales unchanged, unlike `c_unfold`/
  `k_p`, which do not (see below). Don't touch it without a reason from
  the spec.
- `host.c_unfold` (`C`) and `host.k_p` default to Tian's published
  values (0.81, 10) as a prior. **Re-fit both before trusting an absolute
  coupled TIR** — they mix physical ribosome occupancy with whichever au
  scale RBSForge outputs, and Tian's values were fit on a different
  scale. `k_reinitiation(d)` is the one piece that transfers as-is.
- If the downstream start is leaderless, `r_folded` is undefined: treat
  de novo as `0` and keep only re-initiation (a real SD-less overlapping
  architecture) — never substitute `TIR = 0` for a leaderless start
  (known failure mode 14).
- There is no closed-form inverse of this module. Do not build one for
  `operon.design` — coupling-aware design is sample-and-score (spec
  section 18), not inversion (spec section 7.5, "do not invert coupling").
- Iterate 5' to 3', using each step's *already-coupled* `r_{i-1}`, not
  the mono-cistronic prediction — this is why changing an upstream RBS
  can make a downstream RBS "disappear" in Evaluate mode (spec section
  7.3's note, known failure mode 1 is the promoter-DNA-vs-mRNA version of
  the same mistake).

## Interface to build toward

```python
def k_reinitiation(d: int) -> float: ...   # already implemented
def coupled_tirs(assembled, host, calc) -> list[float]: ...
```
