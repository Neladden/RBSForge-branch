# operon.assembly — Module A

Concatenates promoter + `[RBS_i + CDS_i]` + terminator into one assembled
DNA/mRNA molecule and computes each junction's intergenic distance `d`.
**Status: implemented.** Spec section 3.

## Design decisions worth knowing before touching this

- `intergenic_policy` is uniform across the whole `Operon`, not per-
  junction (`Operon` only carries one field for it). If a design ever
  needs different policies at different junctions, that's a breaking
  change to `operon.core.Operon`, not a workaround in `assemble()`.
- **`"overlap-N"` does not try to reconcile two overlapping reading
  frames.** It trims the last N nt already placed and splices in CDS_i's
  own leading N nt as the physically shared span — CDS_i's sequence is
  authoritative for that span, full stop. If you need both frames to
  encode specific, correct proteins across the overlap (the usual case
  for a real AUGA-style junction), that reconciliation is a Design-mode
  joint recode of "the last 5-10 codons of CDS_i plus the first 5 of
  CDS_{i+1}" (spec section 18.3/18.4) — do it before calling `assemble`,
  not inside it. Don't add frame-reconciliation logic here.
- `Assembled.features` spans are on the full DNA (0-based from
  `promoter_dna`'s start) — that's the only coordinate system in which
  the promoter and terminator exist at all. `Assembled.starts` /
  `Assembled.cds_end` are on `Assembled.mrna` instead (0-based from the
  TSS), matching `CDS.annotated_start`'s documented convention and what
  `operon.coupling`'s section 7.3 pseudocode (`assembled.starts`,
  `assembled.cds_end`) expects. Don't mix the two coordinate systems when
  extending this module.
- `Operon.tss` (added alongside this implementation) is `None` by
  default, meaning "TSS is exactly where `promoter_dna` ends." Once
  `operon.promoter_calculator` is wired into an end-to-end pipeline, a
  real TSS call should set this explicitly rather than assuming the
  promoter/UTR boundary is the transcription start.
- Cryptic-promoter and cryptic-terminator scans need the *full* DNA
  (`Assembled.dna`), not `Assembled.mrna` — the promoter and terminator
  regions are in scope for those scans (spec section 3.3). Don't slice
  `Assembled.mrna` for anything that's supposed to see the whole
  construct.

## Tests

`tests/test_assembly.py`. Covers: single-CDS layout and feature spans,
default vs. explicit TSS, mRNA-relative `starts`/`cds_end`, all four
`intergenic_policy` kinds (including the canonical `overlap-4` AUGA-style
junction and `overlap-25` insulation), a 3-cistron operon's two
junctions, and the reject cases (empty `cds_list`, unknown policy, wrong-
length spacer, overlap longer than the downstream CDS).
