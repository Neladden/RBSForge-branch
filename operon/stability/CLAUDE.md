# operon.stability — Module D

mRNA stability / decay rate prediction. **Status: prepared subsection,
not implemented.** Spec section 9. RBSForge explicitly does not include
this — it is entirely this module's responsibility. Consumes coupled
TIRs from `operon.coupling` as ribosome-protection inputs, but is
otherwise independent.

## Implement the 2021 biophysical model first, not the 2024 GBDT

Spec section 9.2 (2021, Cetnar & Salis,
https://doi.org/10.1021/acssynbio.0c00471) is the model to build first —
two interpretable features, one linear fit. Spec section 9.3 (2024,
isoform-aware LightGBM, https://doi.org/10.1038/s41467-024-54059-7) is a
later layer that additionally depends on `operon.promoter_calculator`
(for isoform/TSS selection) and is not a prerequisite for a working
Module D.

```
dx = r_elong / r_init - L_footprint          # nt of naked mRNA between ribosomes
RNase_score = sum(ss[i] * motif_w[nt[i]] for i in 5' UTR)   # unpaired AU-richness
log(mRNA) = a0 + a1*(-RNase_score) + a2*(-dx) + a3*1[5p_hairpin]
k_decay   = k0 * exp(b1*RNase_score + b2*dx)
t_half    = ln(2) / k_decay
```

- `{a, b}` **must be fit on a paired sequence/qPCR set** — the 2021
  paper's 82-operon set is the calibration target. This function is not
  usable with made-up coefficients; don't ship a plausible-looking
  default.
- For multi-cistronic operons, weight `dx` toward the 5' cistron (or use
  the minimum protection across cistrons) — the least-translated cistron
  is the RNase entry point, not a length-weighted average across all of
  them equally.

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
