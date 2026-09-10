# The Promoter Calculator model

`operon.promoter_calculator` reimplements the sigma70 promoter free-energy
model from:

> La Fleur, T., Hossain, A. & Salis, H. M. Automated model-predictive
> design of synthetic promoters to control transcriptional profiles in
> bacteria. *Nat. Commun.* **13**, 5159 (2022).
> https://doi.org/10.1038/s41467-022-32829-5

The reference implementation and its trained coefficients are public:
https://github.com/hsalis/SalisLabCode/tree/master/Promoter_Calculator
(GPL-3.0, La Fleur/Hossain/Salis). This package does not vendor or import
that code. It re-extracts the trained model's **numeric coefficients**
(343 fitted weights + 1 intercept — published, trained parameters, not
this package's own invention) and re-expresses the scoring algorithm as
independently written Python with no runtime dependencies, the same
approach RBSForge takes with the published Xia et al. 1998 nearest-
neighbor parameters (see `rbsforge/docs/MODEL.md`).

## The six elements

For every candidate transcription start site (TSS), walking upstream:

| Element | Length | Role |
|---|---|---|
| ITR (initial transcribed region) | 20 nt, starts at TSS | R-loop / DNA:RNA hybrid strength |
| discriminator | 6–10 nt | contacts sigma70 region 1.2 |
| -10 hexamer | 6 nt | contacts sigma70 region 2.3/2.4 |
| spacer | 15–20 nt | its length alone sets a quadratic penalty; its last 2 nt are the extended -10 / TGn motif |
| -35 hexamer | 6 nt | contacts sigma70 region 4.2 |
| 1 nt gap | 1 nt | between -35 and UP |
| UP element | 24 nt (12 distal + 12 proximal) | minor-groove width + DNA bending rigidity |

Minimum promoter footprint: 78 bp (shortest spacer + discriminator).
Maximum: 89 bp. Both strands are scanned; a circular molecule must be
wrapped by the caller before scanning (matching RBSForge's own stance on
windowing — this package does not implement origin-wrapping).

## Free-energy model

```
dG_10     = HEX10_LEFT[hex10[0:3]] + HEX10_RIGHT[hex10[3:6]]
dG_35     = HEX35_LEFT[hex35[0:3]] + HEX35_RIGHT[hex35[3:6]]
dG_ext10  = EXT10_2MER[spacer[-3:-1]]
dG_disc   = DISC_FIRST3[disc[0:3]]
dG_spacer = 0.1463*s**2 - 4.9113*s + 41.119          # s = len(spacer); analytic, not a lookup
dG_ITR    = coef_itr * (dg_hybrid / NORM_ITR)
dG_UP     = coef_dist*(width_dist/NORM_DIST_UP) + coef_prox*(width_prox/NORM_PROX_UP)
          + coef_rigidity*(rigidity/NORM_RIGIDITY)
dG_total  = dG_10 + dG_35 + dG_disc + dG_ITR + dG_ext10 + dG_spacer + dG_UP + INTERCEPT

Tx_rate   = K * exp(-BETA * dG_total)
```

`dg_hybrid`, `width_dist`/`width_prox`, and `rigidity` are the four
numerical features (`operon/promoter_calculator/features.py`), each a sum
or mean of a dinucleotide-indexed physical parameter table
(`tables.py`): DNA:RNA hybrid energies (Sugimoto et al. 1995) and DNA:DNA
duplex energies for the ITR term, minor-groove accessibility and
persistence-length proxies (Geggier & Vologodskii 2010) for the UP-element
term.

`dG_spacer` is an **analytic quadratic in spacer length**, not a lookup
table — the reference model does compute a 3-parameter categorical
spacer-length term during training, but it is dead code in the reference
implementation (never wired into the final score); we reproduce the
behavior that is actually used, matching the spec's Appendix D, and
confirmed against the reference `.npy` coefficients (see below).

## Constants

```
K              = 42.0                  # both organism presets
BETA_IN_VIVO   = 1.636217004872062     # E. coli MG1655, in vivo
BETA_IN_VITRO  = 0.81632623            # in vitro
```

Only these two conditions are calibrated. Every other host would need its
own MPRA-scale training set to refit the 343 coefficients and beta — there
is no principled "prior" transfer here the way RBSForge documents for
non-*E. coli* HostPacks, because this model was trained end-to-end on raw
DNA sequence rather than composed from separately-measured physical terms.
Do not present `Tx_rate` from a non-`ecoli`/`in_vitro` organism as
meaningful.

## How the coefficients were obtained

`coefficients.py`'s 343 values were extracted directly from the public
release's `free_energy_coeffs.npy` (339 categorical dG values across six
motif tables) and `model_intercept.npy` (1 value), plus 4 coefficients on
the normalized numerical features that live at the end of the same array.
The categorical tables are keyed by the sequence motif they score, in the
same lexicographic 3-mer/2-mer ordering `sklearn.preprocessing.
OneHotEncoder` produces when fit on `itertools.product('ACGT', repeat=k)`
(alphabetical — this package reproduces that ordering directly, with no
`sklearn` dependency at runtime). This was verified end-to-end: a
consensus sigma70 promoter (TTGACA … 17 bp spacer … TATAAT) scores
strongly negative `dG_total` / high `Tx_rate` relative to a scrambled
control, and a sequence's reverse complement produces the exact same
`dG_total` when scanned on the minus strand as the original scores on the
plus strand (see `tests/test_promoter_calculator.py`).

## Known limitations (inherited from the reference model)

- sigma70 only. Other sigma factors (sigma19/24/28/32/38/54) and T7 RNAP
  are different, unfitted models — swapping this module, not extending
  it (spec section 10.4, Appendix G).
- No isoform/TSS-selection layer — every TSS above threshold is reported
  independently; the 2024 mRNA-stability model's "top 5 isoforms" step
  (spec section 9.3) is Module D's responsibility, not this one's.
- Circular molecules: wrap the sequence before scanning; this package
  does not wrap the origin for you.
- `Tx_rate` is a rate on this model's own fitted scale — do not mix it
  with RBSForge's `v1_style_rate` or `proportional_rate` in a single
  threshold or objective (spec's "known failure mode" 10, restated for a
  different pair of scales).
