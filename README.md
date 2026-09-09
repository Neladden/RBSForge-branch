# RBSForge

A thermodynamics-based (not ML) ribosome binding site (RBS) strength
predictor, reimplementing the Salis Lab RBS Calculator's free-energy model
from scratch for organisms outside its *E. coli*-tuned defaults.

Given an mRNA sequence, a species, and a temperature, it predicts a
translation initiation rate for every candidate start codon, with a full
breakdown of the underlying free-energy terms — no black box.

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

## Live demo

`web/rbsforge-console.html` is a self-contained, dependency-free
in-browser reimplementation of the same model (a direct JS port of
`rbsforge/`, numerically verified against the Python package) with a live
UI: type a sequence, pick a host organism, and drag the temperature
slider to see the six-term ΔG breakdown and the SD:anti-SD alignment
recompute instantly. Open it directly in any browser — no server, no
build step.

See `docs/MODEL.md` for the full formula set, citations, and — important —
the limitations around what's calibrated (E. coli, ~37 C) versus what's a
documented but unfit prior (every other built-in organism/temperature).

## Install (dev)

No required third-party dependencies. Optional: install
[ViennaRNA](https://www.tbi.univie.ac.at/RNA/) and have `RNAfold` on PATH
for a higher-accuracy folding backend; the package falls back to a
built-in simplified folder otherwise.

```
python -m unittest discover -s tests
```
