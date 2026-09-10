# operon.assembly — Module A

Concatenates promoter + `[RBS_i + CDS_i]` + terminator into one assembled
DNA/mRNA molecule and computes each junction's intergenic distance `d`.
**Status: prepared subsection, not implemented.** Spec section 3.

## Before implementing

- This is the first module to build (spec section 20, step 1) — nothing
  else in `operon.*` has a hard dependency on it being done first, but
  every other module's function signatures take an `Assembled` (from
  `operon.core`), so get its shape right early.
- `d = start_{i+1} - stop_i_end`. Negative means overlap. Compute it once
  per junction at assembly time and store it on `Assembled.junctions`;
  don't recompute it downstream (`operon.coupling` needs the exact same
  value `operon.stability`'s protection-window math would derive
  independently otherwise).
- The mRNA submitted to translation/coupling/stability starts at the TSS
  and ends at the terminator's U-tract — not the full assembled DNA.
  Cryptic-promoter and cryptic-terminator scans run on the *full DNA*,
  both strands, including the promoter and terminator regions themselves
  (spec section 3.3) — don't trim those before handing DNA to
  `operon.promoter_calculator` or `operon.terminators`.
- Default junction policies (spec section 3.2): `d <= -25` for
  independent/insulated expression; `d == -4` (AUGA) for hairpin-gated
  coupling; `d` in `{-4, -1, +3}` with a strong SD and no inhibitory
  structure for a simple linear leak. Don't default to a long
  unstructured spacer for insulation — it kills coupling *and* creates an
  RNase E site *and* may create a cryptic promoter (known failure mode 4).

## Interface to build toward

```python
def assemble(operon: Operon) -> Assembled: ...
```
(signature and field-level docstring already in `__init__.py`).

## Tests to write alongside the implementation

- Round-trip a known bi-cistronic construct: AUGA junction should give
  `d == -4` (spec section 20's own suggested test for this step).
- `Assembled.mrna` starts at the TSS, not at the start of `Assembled.dna`.
- Overlap, abutting, and spacer junction configurations all produce the
  `d` the table in spec section 3.1 predicts.
