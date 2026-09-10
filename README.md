# Operon Calculator

A from-scratch recreation of the Salis Lab Operon Calculator: an assembler,
a translation-initiation-rate (TIR) engine, a translational-coupling model,
and a set of independent scanners/recoders for everything else a multi-CDS
bacterial operon needs — codon optimization, mRNA stability, cryptic
promoters, internal terminators, ribosomal pause sites, repeats, and
synthesis/restriction-site constraints — composed behind one assembler and
(eventually) one multi-objective design search.

**RBSForge is one module of this project, not the whole of it.** It is the
Predict-mode TIR engine (`rbsforge/`) that the rest of the stack calls as
its translation-initiation oracle. Everything else — assembly, coupling,
codon recoding, mRNA stability, the promoter calculator, terminator/repeat/
synthesis scanners, and inverse RBS + operon design — lives under
`operon/`, one subpackage per module.

See `docs/OPERON_CALCULATOR.md` for the architecture and module status,
and `docs/OPERON_CALCULATOR_SPEC.md` for the full design spec these
subpackages implement against.

## Layout

```
rbsforge/                       RBS Calculator recreation (TIR engine; Predict mode)
  docs/MODEL.md                 RBSForge's own free-energy model doc
  docs/HOSTPACK_DECODER.md      draft design doc: predicting HostPacks from a genome embedding
  web/rbsforge-console.html     live, in-browser demo (self-contained, no build step)

operon/                         everything else in the Operon Calculator
  core/                         shared types (OperonHost, CDS, Operon, Assembled)
  assembly/                     Module A — operon assembly, intergenic distance d
  coupling/                     Module B — translational coupling
  elongation/                   Module C — TER + synonymous codon recoding
  stability/                    Module D — mRNA stability
  promoter_calculator/          Module E — sigma70 promoter / cryptic-promoter scan
    docs/MODEL.md                promoter calculator's own free-energy model doc
    web/promoter-calculator-console.html   live, in-browser demo (self-contained, no build step)
  htisc/                        Module F — highly translated internal start codons
  pauses/                       Module G — ribosomal pause sites
  terminators/                  Module H — intrinsic + rho-dependent terminators
  repeats/                      Module I — repeats, IS/att sites
  synthesis/                    Module J — synthesis complexity, restriction sites
  design/                       inverse RBS design (design_rbs) + NSGA-II operon design

docs/
  OPERON_CALCULATOR.md          architecture overview + module status table
  OPERON_CALCULATOR_SPEC.md     full design spec (source of truth for every module)

examples/
tests/
```

## Status

| Module | What | Status |
|---|---|---|
| RBSForge | TIR engine (Predict mode) | implemented |
| `operon.core` | shared types | implemented |
| `operon.promoter_calculator` (Module E) | sigma70 promoter scan | implemented |
| everything else under `operon/` | see table above | prepared subsection, not yet implemented |

The promoter calculator is the current focus — a from-scratch, dependency-
free reimplementation of the 346-parameter LaFleur, Hossain & Salis (2022)
linear free-energy model for sigma70 transcription rate. See
`operon/promoter_calculator/docs/MODEL.md`.

## RBSForge

A thermodynamics-based (not ML) ribosome binding site (RBS) strength
predictor, reimplementing the Salis Lab RBS Calculator's free-energy model
from scratch for organisms outside its *E. coli*-tuned defaults.

```python
from rbsforge import predict

result = predict(
    "AGGAGGACAACTAAATGAAACGCATTAGCACCACC...",
    species="ecoli",       # or b_subtilis, geobacillus, thermus, generic
    temperature_c=37.0,
)
best = result.best
print(best.codon, best.delta_g_total, best.translation_initiation_rate)
print(best.breakdown)
```

Command line:

```
python -m rbsforge.cli "AGGAGG...ATGAAACGC..." --species ecoli --temperature 37
```

Live demo: open `rbsforge/web/rbsforge-console.html` directly in a browser
for an interactive version — no server, no build step.

See `rbsforge/docs/MODEL.md` for the full formula set, citations, and —
important — the limitations around what's calibrated (*E. coli*, ~37 C)
versus what's a documented but unfit prior (every other built-in
organism/temperature).

## Promoter Calculator

A sigma70 promoter scanner: given a DNA sequence, scores every candidate
transcription start site (TSS), on both strands, for its UP element, -35
hexamer, spacer, -10 hexamer, discriminator, and initial transcribed
region, and reports a transcription rate.

```python
from operon.promoter_calculator import predict, scan_promoters

result = predict("...promoter DNA...", organism="ecoli")
best = result.best
print(best.tss, best.strand, best.tx_rate, best.dg_total)

# cryptic promoters relative to a designed/annotated TSS
cryptic = scan_promoters(assembled_dna, intended_tss=my_tss, tau_tx=1.0)
```

Live demo: open `operon/promoter_calculator/web/promoter-calculator-console.html`
directly in a browser for an interactive version — no server, no build
step; same design and JS-port approach as the RBSForge console above.

See `operon/promoter_calculator/docs/MODEL.md` for the formula set,
citations, and calibration notes.

## Install (dev)

No required third-party dependencies for either `rbsforge` or
`operon.promoter_calculator`. Optional: install
[ViennaRNA](https://www.tbi.univie.ac.at/RNA/) and have `RNAfold` on PATH
for a higher-accuracy folding backend in RBSForge; it falls back to a
built-in simplified folder otherwise.

```
python -m unittest discover -s tests
```
