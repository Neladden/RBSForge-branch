# operon.terminators — Module H

Internal transcriptional terminator scanning: intrinsic (rho-independent)
and rho-dependent. **Status: implemented** (presence layer only for
both). Spec section 13. There is no Salis Terminator Calculator to port
— this is a scan plus a parts-database append, not a fitted rate model
of its own.

## What's not implemented

- The Chen-class strength layer (dG_U 8-nt U-tract hybrid energy, dG_L
  loop-closure energy) is not implemented — `scan_intrinsic_terminators`
  uses the crude presence gate instead (stem/loop geometry + a >=5-in-8
  U-tract count + hairpin MFE <= -7 kcal/mol via `BuiltinFolder`), which
  the spec itself calls "enough for the Operon count."
- `scan_rho_terminators` is a simplified Rut-site presence check (C/G >
  1, C count consistent with "roughly every 11-13 nt" in a 78-nt window,
  plus a downstream pause), not a trained RhoTermPredict/OPLS-DA
  classifier. Expect it to report several overlapping hits around one
  real Rut-like region (consecutive shifted 78-nt windows that all still
  qualify) rather than one merged hit — no de-duplication/merging pass
  exists yet.
- Neither scanner is fast on genome-scale input: `_find_hairpin_at` tries
  every stem/loop combination at every position (~0.1-0.2s per 1000 nt
  for the intrinsic scan in this implementation). Fine for operon-scale
  DNA (hundreds to a few thousand nt); don't reach for this on a whole
  plasmid or genome without profiling first.

## Things worth knowing if you touch this

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
