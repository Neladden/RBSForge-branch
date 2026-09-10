# operon.stability — Module D

mRNA stability / decay rate prediction. **Status: features implemented
(dx, RNase-accessibility score, 5' hairpin detection); the final fitted
rate is deliberately gated on real coefficients.** Spec section 9.
RBSForge explicitly does not include this — it is entirely this module's
responsibility. Consumes coupled TIRs from `operon.coupling` as
ribosome-protection inputs, but is otherwise independent.

## Why `mrna_stability` can return `k_decay=None`

```
dx = r_elong / r_init - L_footprint          # nt of naked mRNA between ribosomes
RNase_score = sum(ss[i] * motif_w[nt[i]] for i in 5' UTR)   # unpaired AU-richness
log(mRNA) = a0 + a1*(-RNase_score) + a2*(-dx) + a3*1[5p_hairpin]
k_decay   = k0 * exp(b1*RNase_score + b2*dx)
t_half    = ln(2) / k_decay
```

`dx` (`ribosome_protection_dx`) and `RNase_score`
(`rnase_accessibility_score`) are implemented for real — both are
physically-grounded and need no fitted data beyond `host.r_elong_nominal_nt_s`/
`host.ribosome_footprint_nt` (already on `OperonHost`) and a fold. The
final `log(mRNA)`/`k_decay` combination is **not** wired to a default
`{a, b}` — the spec is explicit that those must be fit on a real paired
sequence/qPCR set (the 2021 paper's 82-operon set is the calibration
target), and this repo's convention is not to invent coefficients that
were never published. `mrna_stability(..., coefficients=None)` (the
default) returns only the features; pass a real fit as
`coefficients={"a0":..., "a1":..., "a2":..., "a3":..., "k0":..., "b1":...,
"b2":...}` to get `k_decay`/`t_half` back. **Never add a
`DEFAULT_COEFFICIENTS` dict to this module** — `tests/test_stability.py`
has a regression test (`test_never_invents_coefficients_of_its_own`)
checking exactly that name doesn't appear in the source; if you have a
real fit, wire it through a `HostPack`/`OperonHost` field or a caller-
supplied argument, not a hardcoded module-level default.

- `ribosome_protection_dx` converts au TIR via `operon.coupling.to_physical`
  — the same documented identity seam, same reasoning (no published
  au -> s^-1 formula exists). Keep both in sync if that seam ever changes.
- `rnase_accessibility_score`'s per-base weights (`_RNASE_MOTIF_WEIGHT`:
  A/U=1.0, G/C=0.3) are directionally motivated by the paper's own
  finding (AU-rich unpaired tracts destabilize, GC protects), not a
  fitted table — same honesty tier as `operon.elongation`'s placeholder
  codon weights, not a source of real absolute scores.
- Multi-cistronic aggregation uses the **minimum** `dx` across cistrons
  (`StabilityResult.per_cistron_dx` has every cistron's value if you need
  the length-weighted-average alternative instead) — the least-
  translated cistron is the RNase entry point, spec section 9.2's stated
  preference.
- Spec section 9.3's 2024 isoform-aware LightGBM layer (needs
  `operon.promoter_calculator` for TSS/isoform selection) is a later,
  separate expansion, not a prerequisite for this module.

## The one rule that matters most here

**Do not put terminator strength or 3' UTR ssRNA into this module's
half-life formula.** Cetnar & Salis 2021 found no effect on upstream
mRNA for either, in *E. coli*. This is spec section 9.1's explicit
finding and known failure mode 6 — it will be tempting to add a
terminator-strength term because `operon.terminators` exists and
produces a convenient number; don't.

Design-mode actions this module should expose as mutators (section 9.4):
raise CDS1 TIR (don't accidentally under-express the protective ribosome
train), insert a small 5' hairpin (6-8 bp stem, 4-6 nt loop) that doesn't
occlude the standby site, shorten unpaired AU tracts in the 5' UTR. Do
not add "massage the terminator" as an option.

## Interface to build toward

```python
def mrna_stability(assembled, tirs: list[float], host) -> StabilityResult: ...
```
