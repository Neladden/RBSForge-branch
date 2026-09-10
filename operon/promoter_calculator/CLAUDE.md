# operon.promoter_calculator — Module E

sigma70 promoter transcription-rate scanner: both strands, every
candidate TSS. **Status: implemented.** Spec section 10, Appendix D.

Read `docs/MODEL.md` (in this directory) before changing the energy
model — it has the full formula set, the coefficient provenance, and why
`dG_spacer` is an analytic quadratic rather than the categorical lookup
the reference model's coefficient file technically has room for.

## File map

- `tables.py` — raw dinucleotide physical parameter tables (groove
  width, rigidity, DNA:DNA / RNA:DNA hybrid energies). Static published
  data; don't hand-edit without a citation.
- `coefficients.py` — the 343 fitted weights + intercept, extracted
  verbatim from the public reference release. **Never hand-edit these
  values or their key ordering** (lexicographic 3-mer/2-mer, matching
  `sklearn.OneHotEncoder` fit on `itertools.product('ACGT', repeat=k)` —
  see `docs/MODEL.md` for how they were pulled and how to re-pull them if
  the upstream model is ever updated).
- `features.py` — the four numerical features (groove width x2,
  rigidity, DNA:RNA hybrid) as functions over `tables.py`.
- `energy_model.py` — `score_elements(SigmaFactorElements) ->
  PromoterEnergyBreakdown`. Pure function, no scanning logic.
- `calculator.py` — `PromoterCalculator`, the TSS x
  discriminator-length x spacer-length scan, both-strand handling, and
  the top-level `predict`/`scan_promoters` convenience functions.
- `types.py` — `PromoterHit`, `PromoterScanResult`.
- `web/promoter-calculator-console.html` — self-contained, dependency-free
  in-browser reimplementation of this package (a direct JS port,
  numerically verified against it: same random/consensus sequences give
  identical `dG_total`/`Tx_rate` to floating-point noise, across both
  organisms and both strands — see the file's own model-section comment
  for the parity-check pattern). Sibling to `rbsforge/web/rbsforge-console.html`,
  same design system. Opens directly in any browser, no server or build
  step. **If you change a formula in this package, the console's JS port
  drifts out of sync silently — there is no shared source between them.**
  Re-verify the port before merging a model change.

## Load-bearing constraints

- **Only two conditions are calibrated**: `organism="ecoli"` (in vivo,
  MG1655, `BETA_IN_VIVO`) and `organism="in_vitro"`. There is no
  principled prior-transfer story for a third host the way RBSForge has
  for non-*E. coli* HostPacks — this model was trained end-to-end on raw
  DNA sequence, not composed from separately-measured physical terms.
  Don't add a third organism preset without a real refit.
- **`Tx_rate` is on this model's own scale.** Never compare or threshold
  it against RBSForge's `v1_style_rate` or `proportional_rate` in the
  same objective/threshold (the general form of known failure mode 10).
- **sigma70 only.** Other sigma factors and T7 RNAP need a different,
  unfitted model — that's a module swap (a new subpackage), not an
  extension of this one (spec section 10.4, Appendix G).
- **Circular molecules are the caller's problem.** This scanner does not
  wrap the origin; `operon.assembly` (once built) or the caller must
  concatenate/wrap before calling `predict`.
- **Reverse-strand coordinates are flipped back to the original,
  forward-strand numbering** (`_flip_span` in `calculator.py`) so every
  `PromoterHit`, on either strand, reports spans in one consistent
  coordinate system matching `operon.core`'s convention. If you touch the
  reverse-strand path, re-run
  `test_reverse_strand_hit_has_identical_energy_to_forward` — it asserts
  exact energy equality between a forward scan and the same promoter
  embedded as a reverse complement, which is the real invariant here (not
  just "doesn't crash").
- `scan_promoters(dna, intended_tss=None, ...)` returns *every* hit above
  `tau_tx` when `intended_tss` is `None` (Predict/Evaluate mode listing
  everything found); with `intended_tss` set, it returns only the
  cryptic/internal hits (spec section 10.3). Don't change the default
  behavior of the `intended_tss=None` case — other modules
  (`operon.stability`'s 2024 isoform layer, eventually) will call this in
  "list everything" mode.

## Tests

`tests/test_promoter_calculator.py`, from the repo root:
`python -m unittest tests.test_promoter_calculator -v`. Covers
consensus-vs-scrambled hexamer discrimination, breakdown-terms-sum
invariant, span-recovers-hexamer sanity, the reverse-strand exact-energy
round-trip, and `scan_promoters`'s intended-TSS exclusion.
