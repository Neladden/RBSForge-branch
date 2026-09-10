# operon.terminators — Module H

Internal transcriptional terminator scanning: intrinsic (rho-independent)
and rho-dependent. **Status: prepared subsection, not implemented.** Spec
section 13. There is no Salis Terminator Calculator to port — this is a
scan plus a parts-database append, not a fitted rate model of its own.

## Two independent detectors, don't conflate them

**Intrinsic** (spec section 13.2): GC-rich hairpin + U-tract.
- Presence layer (enough for the Operon count): TransTermHP-style —
  stem 4-15 bp, loop 3-8, U-tract ~15 nt, `score = hairpin MFE + U-tract
  RNA:DNA hybrid energy`. Kingsford, Ayanbule & Salzberg, Genome Biology
  2007 (https://doi.org/10.1186/gb-2007-8-2-r22).
- Strength layer (optional, better as a weight than a binary): Chen et
  al. 2013 (https://doi.org/10.1038/nmeth.2515) — `dG_U` (8-nt U-tract
  RNA:DNA hybrid) and `dG_L` (loop-closure energy) are what tracked
  termination strength empirically; **hairpin MFE alone did not**. Don't
  substitute MFE for `dG_L` if you implement the strength layer.

**Rho-dependent** (spec section 13.3): no Salis model exists; use a
published motif detector — RhoTermPredict (Di Salvo et al. 2019,
https://doi.org/10.1186/s12859-019-2704-x: 78-nt Rut, C/G > 1, C every
11-13 nt, max C/G in a 128-nt neighborhood, then a pause within
~50-100 nt) or Nadiras et al. 2018
(https://doi.org/10.1093/nar/gky705: C>G skew + regular CC/UC spacing +
low structure). Rho hits are common in high-CAI, C-rich recoding — this
is one of the reasons `operon.elongation` should default to the Balanced
codon table.

## The one rule that matters most here

**A hit at the user-supplied 3' terminator is not "internal."** Allow a
~50 bp window around the annotated terminator position before counting a
hit as a real internal terminator (spec section 13.2). This module also
feeds `operon.stability`, and Cetnar & Salis 2021 found terminator
strength has **no effect on upstream mRNA half-life** — don't let this
module's output leak into `operon.stability`'s formula (known failure
mode 6; see `operon/stability/CLAUDE.md`).

## Gotchas

- Scan both strands. Forgetting the reverse strand is explicitly called
  out as a known failure mode (spec section 24, item 8).
- If the operon has no user-supplied terminator, appending one is a
  Design-mode action (spec section 13.4: a Chen-class intrinsic
  terminator from a non-repetitive toolbox, `L_max=12` vs. the rest of
  the operon, checked against `operon.repeats` for 12-mer clashes) — a
  part-selection decision, not this module's own detector logic.

## Interface to build toward

```python
def scan_intrinsic_terminators(dna: str) -> list: ...
def scan_rho_terminators(dna: str) -> list: ...
```
