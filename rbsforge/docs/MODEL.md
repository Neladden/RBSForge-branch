# The RBSForge model

RBSForge reimplements the free-energy model behind the Salis Lab RBS
Calculator (Predict mode only) from scratch, in pure Python with no
required external dependencies, parameterized so it can be pointed at
organisms other than *E. coli*.

## Primary sources

1. Salis, Mirsky & Voigt, *Nat. Biotechnol.* **27**, 946 (2009) — v1.0, the
   five-term equilibrium free-energy model and its `K * exp(-dG/RT_eff)`
   Boltzmann map.
2. Salis, *Methods Enzymol.* **498**, 19 (2011) — v1.1, final-state
   spacing minimization and CDS-footprint unfolding.
3. Espah Borujeni, Channarasappa & Salis, *Nucleic Acids Res.* **42**, 2646
   (2013) — v2.0, standby-site module geometry.
4. Farasat *et al.*, *Mol. Syst. Biol.* **10**, 731 (2014) — multi-host use
   via anti-SD swapping.
5. Espah Borujeni & Salis, *J. Am. Chem. Soc.* **138**, 7016 (2016) —
   ribosome drafting (kinetic override; not implemented here).
6. Espah Borujeni *et al.*, *Nucleic Acids Res.* **45**, 5437 (2017) —
   ribosome footprint definition, ΔG_mRNA:rRNA decomposition.
7. Reis & Salis, *ACS Synth. Biol.* **9**, 3145 (2020) — v2.1 refit
   (stacking term, Gram-positive spacing/standby, start-codon table,
   Andronescu 2007 NN parameters).

v1.0 source (Python 2 + NuPACK):
[github.com/hsalis/Ribosome-Binding-Site-Calculator-v1.0](https://github.com/hsalis/Ribosome-Binding-Site-Calculator-v1.0).
v2.x source is not public; every coefficient beyond v1.0/v2.0 used here is
the published number, not a decompiled implementation.

## What is implemented

The six-term master equation:

```
delta_G_total = delta_G_standby + delta_G_mRNA:rRNA + delta_G_spacing
                + delta_G_start + delta_G_stacking - delta_G_mRNA

r = exp(-beta * delta_G_total)
```

| Term | Status | Notes |
|---|---|---|
| `delta_G_mRNA` | implemented (v1.1 two-fold) | Bound-state unfolding cost: `E_bound - E_initial`, where `E_initial` is the unconstrained MFE of the folded window and `E_bound` is the MFE of the same window with the start codon plus `HostPack.footprint_cds` (13 nt) of downstream CDS forced unpaired -- the 30S occupancy footprint. Structure entirely outside that occupied region costs nothing via this term (it doesn't have to be disrupted for the ribosome to bind); structure entirely inside it costs its full stability. Via a from-scratch simplified nearest-neighbor folder, or ViennaRNA `RNAfold` if installed on PATH. Exposed as `RBSCalculator.predict_one(mrna, start, extra_unpaired_global=None)`, the one constrained-fold hook translational coupling (in `operon.coupling`) also uses. |
| `delta_G_mRNA:rRNA` | implemented | SD:anti-SD duplex search over all registers, scored with temperature-scaled nearest-neighbor stacking. Uses the v1.0 selection rule: keep the alignment that minimizes `delta_G_hybrid + delta_G_spacing(aligned_spacing)` jointly, not the strongest duplex alone. |
| `delta_G_spacing` | implemented | Exact v1.0 push/pull formulas (logistic compression penalty, quadratic extension penalty), evaluated on *aligned* spacing (Chen/Salis-corrected for anti-SD registers that don't reach the tail's 3' terminus), not raw nucleotide count. |
| `delta_G_start` | implemented | v1.0 codon table (AUG/GUG/UUG/CUG apparent fMet-tRNA pairing energies). |
| `delta_G_standby` | implemented (v1.0 form) | 4-nt forced-unpaired method: the energetic cost of forcing the 4 nt upstream of the SD duplex to stay single-stranded. The richer v2.0 module geometry (`StandbyParams` in `hostpack.py`) is scaffolded but not wired up. |
| `delta_G_stacking` | **not implemented** | v2.1 homopolymer-spacer correction; its coefficient was never published (Reis & Salis 2020 supplement). Reported as an explicit `0.0` term rather than silently omitted. |
| Ribosome drafting | **not implemented** | Kinetic override of the unfolding penalty (Espah Borujeni & Salis 2016) requires folding-time estimates this package does not compute. |

Temperature enters through the standard van't Hoff decomposition applied
to nearest-neighbor stacking energies:

```
delta_G(T) = delta_H - T * delta_S,  delta_S = (delta_H - delta_G_37) / 310.15
```

using the unified RNA nearest-neighbor parameters (Xia et al. 1998,
consistent with Turner rules). Heat-capacity (ΔCp) effects are neglected —
a good approximation to ~60 C, progressively worse above that (see the
uploaded model summary's worked example, reproduced qualitatively in
`tests/test_nn_params.py`).

## Beta is not portable across temperature or organism without refitting

`r = exp(-beta * delta_G_total)` is a Boltzmann map, but **beta is an
empirical fit constant, not 1/RT.** At 37 C, the published beta (0.45
mol/kcal, equivalent to RT_eff = 2.22 kcal/mol) is ~3.6x smaller than the
ideal-solution 1/RT — the gap is attributed to molecular crowding and
unmodeled interactions (IFs, S1, Mg2+, supercoiling), and it was fit
almost entirely on *E. coli* near 37 C.

**Do not treat absolute `translation_initiation_rate` values from a
non-`ecoli` HostPack, or from `ecoli` at a temperature far from 37 C, as
calibrated.** Every HostPack in `hostpack.py` besides `ecoli` inherits
beta/spacing/standby coefficients from the *E. coli* fit as an explicit,
documented *prior* — this is enough to get correct **rankings** among
sequences scored with the same HostPack (aSD swapping alone recovers most
of the useful signal per Farasat et al. 2014), but not a validated
absolute rate. Getting real absolute numbers for a new organism requires
a calibration library: several dozen insulated reporter constructs
spanning SD strength, spacer length, start codon, and structure, run in
that host at its growth temperature, fit against `ln(protein/mRNA) =
const - beta * delta_G_total(theta, T)`.

## HostPack: adding a new organism

See the module docstring in `rbsforge/hostpack.py`. In short:

1. Get the *mature* 16S rRNA 3' tail — not a raw GenBank annotation, which
   is frequently truncated before the anti-SD region. Use RNA-seq 3'-end
   mapping when available.
2. Record growth temperature and Gram stain.
3. Leave `beta`, `spacing_push`/`spacing_pull`, and `standby` at the
   *E. coli* values until you have calibration data (above).

## Known limitations (inherited from the original model, not fixable by a better implementation)

- Equilibrium 30S occupancy is not protein yield when elongation,
  termination, or burden dominate.
- No initiation-factor dynamics, no 50S joining, no tmRNA.
- No RNA-binding proteins (CsrA, Hfq, ...).
- No pseudoknots or G-quadruplexes.
- Leaderless and SD-independent initiation (common in some phyla,
  including many archaea and some Bacteroidetes) are different physics
  this model does not represent; such start codons are reported as
  "leaderless" and excluded from ranked results rather than mis-scored.
- Proportional TIR is not an absolute, cross-lab-comparable rate.

## What this package deliberately omits relative to the full product

- **Design mode** (inverse design of a synthetic RBS for a target TIR via
  simulated annealing) — not implemented; this package is Predict-mode
  only.
- **mRNA-stability sidecar** (Cetnar & Salis 2020) — organism-specific
  RNase biology (RNase E/G vs J/Y); explicitly out of scope.
- **Promoter strength prediction** — out of scope for this package; the
  uploaded technical summary this implementation is built from covers
  translation initiation only.
