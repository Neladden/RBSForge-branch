# rbsforge

The translation-initiation-rate (TIR) engine: reimplements Predict mode of
the Salis Lab RBS Calculator from the published v1.0-v2.1 papers plus the
public v1.0 GPL source. **Status: implemented.** This is the one module
every other `operon.*` subpackage that touches translation calls into.

Read `docs/MODEL.md` (in this directory) for the full formula set and
citations before changing anything here. Read the top-level
`docs/OPERON_CALCULATOR_SPEC.md` section 4 ("RBSForge as the TIR engine")
for how the rest of the Operon Calculator depends on this package's exact
API and TIR scale.

## Public API

```python
from rbsforge import predict, RBSCalculator, HostPack, get_hostpack
result = predict(mrna_sequence, species="ecoli", temperature_c=37.0)
```

`PredictResult.results` is a `List[StartCodonResult]`; `.best` and
`.ranked()` exclude leaderless starts. See `docs/MODEL.md` for the full
field list.

## Live console

`web/rbsforge-console.html` is a self-contained, dependency-free
in-browser reimplementation of this package (a direct JS port of
`rbsforge/`, numerically verified against it — same test sequence gives
identical `delta_G_total`/TIR, including `v1_style_rate`, across every
built-in HostPack and temperature) with a live UI: sequence input, an
organism selector showing which HostPack fields are `calibrated` vs. an
unfit `prior`, editable-but-locked HostPack parameters gated behind a
confirm step, and a temperature slider. Opens directly in any browser, no
server or build step. **If you change a formula in this package, the
console's JS port drifts out of sync silently — there is no shared
source between them.** Re-verify the port (see the file's own comments
for the parity-check pattern) before merging a model change.

## HostPack Decoder

`docs/HOSTPACK_DECODER.md` is a draft design document (not implemented)
for a network that would predict `HostPack` parameters de novo from a
genome embedding, instead of by hand curation. Read it before adding a
new built-in HostPack by hand if the organism in question might be in
scope for that project instead — it lays out which fields are honestly
learnable from sequence alone and which need real calibration data, a
distinction worth checking before hand-curating another `prior`-tier
guess.

## Load-bearing constraints — do not change without reading the spec first

- **`beta` is an empirical fit, not `1/RT`.** Only the `ecoli` HostPack at
  ~37 C is calibrated. Every other built-in HostPack inherits
  `beta`/spacing/standby from the *E. coli* fit as a documented *prior*
  (rankings are defensible; absolute TIR is not). Don't "fix" a non-ecoli
  HostPack's absolute output without a real calibration library (spec
  section 4.1, `docs/MODEL.md`'s "Beta is not portable" section).
- **TIR scale is frozen at `v1_style_rate`** (`2500 * exp(-dG/2.222)`) for
  everything downstream that thresholds or mixes rates (HTISC, coupling,
  design). `translation_initiation_rate` (the `exp(-beta*dG)` proportional
  scale) exists but must never be compared against `v1_style_rate` or
  against published occupancy constants like Tian's `C=0.81` — see spec
  section 4.5 and known failure mode 10/12.
- **`delta_G_stacking = 0.0` is deliberate.** The v2.1 coefficient was
  never published. Do not invent one.
- **`HostPack.footprint_cds = 13` is wired in** via
  `RBSCalculator.predict_one(mrna, start, extra_unpaired_global=None)`
  (spec section 6.2/6.5/Appendix B): the `mRNA` breakdown term is now the
  bound-state unfolding cost `E_bound - E_initial` (start codon +
  `footprint_cds` nt forced unpaired), not the old blanket
  `-1 * unconstrained MFE`. `predict()` calls this for every start codon
  it finds — there is exactly one constrained-fold code path. `operon.coupling`
  reuses this same method via `extra_unpaired_global` (the upstream CDS's
  nucleotides that fall in the downstream start's window) — do not add a
  second, parallel constrained-fold mechanism there.
- **v2.0 standby (`StandbyParams` on `HostPack`) is scaffolded, not
  wired.** Only the v1.0 4-nt forced-unpaired form is implemented. Do not
  ship both v1.0 and v2.0 additively if you wire this in (spec section
  19, 23).
- **No dependency on ViennaRNA is required.** The builtin folder
  (`thermo/fold.py`) is a from-scratch nested-helix folder (no bulges,
  internal loops, or multiloops) — conservative, fine for ranking isolated
  RBS strength, *not* fine for coupling-relevant intergenic hairpins.
  `RNAfold` is used automatically if on `PATH`. Any module that needs
  accurate structure for bulge-containing hairpins (coupling) must check
  `folder.backend == "vienna"` and refuse or warn otherwise, not silently
  fall back to the builtin folder.

## What this package will never do (by design, not by gap)

No mRNA-stability sidecar, no promoter prediction, no terminator model, no
codon recoding, no Design mode (inverse RBS design). Those belong in
`operon.*` subpackages and call into this one; do not add them here.
`design_rbs` in particular is planned in `operon/design/`, calling
`rbsforge.predict` as its energy oracle — see `operon/design/CLAUDE.md`.

## Tests

```
python -m unittest discover -s tests -p 'test_*.py'
```
(from the repo root — `tests/` is shared across the whole repo, not
per-subpackage; files here are named after what they cover:
`test_calculator.py`, `test_duplex.py`, `test_hostpack.py`,
`test_nn_params.py`, `test_spacing.py`, `test_standby.py`,
`test_start_codons.py`.)
