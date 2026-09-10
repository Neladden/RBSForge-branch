# Operon Calculator Internals

A specification for replicating the Salis Lab Operon Calculator around RBSForge, including the missing inverse RBS designer, CDS-footprint unfolding, terminator scans, and every periphery module.

This document merges: the published Operon Calculator product and papers; RBSForge as it actually exists (Predict-mode TIR engine); the gaps between them; and implementable algorithms for `design_rbs`, footprint unfolding, coupling, and terminator analysis.

Full URLs are written out as text, not markdown hyperlinks.

---

## 0. How to read this

The Operon Calculator is not one model. It is a pipeline of independent sequence-to-function models plus a multi-objective search. The only translation-initiation physics that is not already in an RBS calculator is translational coupling. Everything else is a scanner, a recoder, or an optimizer sitting on top of assembled DNA/mRNA.

RBSForge is a from-scratch Predict-mode reimplementation of the Salis RBS Calculator (v1.0–v2.1 papers plus public v1.0 source). It is a legitimate TIR engine for E. coli near 37 C. It is not a black box, it is not a designer, and it does not implement CDS-footprint unfolding even though `HostPack.footprint_cds = 13` exists.

Use this document as:

1. Architecture of the Operon product.
2. Adapter from RBSForge’s real API to the Operon TIR oracle.
3. Specs for every periphery module, independently.
4. Specs for the two holes that block a faithful Operon Design mode: inverse RBS design, and footprint / constrained-fold.
5. Implementation order and failure modes.

Do not reimplement the RBS free-energy model. Call RBSForge. Do not invent the unpublished ΔG_stacking coefficient. Do not copy Tian & Salis occupancy constants onto RBSForge’s au scale without re-fitting.

---

## 1. What the Operon Calculator is

Primary product pages:

- https://docs.denovodna.com/docs/operon-calculator
- https://docs.denovodna.com/docs/operon-calculator.md
- https://github.com/hsalis/salis-lab-protocol-book/blob/master/design/operon-calculator.md
- https://www.denovodna.com/software/predict_operon_calculator
- https://www.denovodna.com/software/design_operon_calculator
- https://salislab.net/software/design_operon_calculator

Core papers:

- Tian & Salis, Nucleic Acids Research 43:7137–7151 (2015), translational coupling. https://doi.org/10.1093/nar/gkv635
- Cetnar, SEED 2016 abstract, 15 models + GA. https://proceedings.aiche.org/sbe/conferences/synthetic-biology-engineering-evolution-design-seed/2016/proceeding/paper/operon-calculator-automated-design-synthetic-operon-sequences-using-15-models-and-design-rules
- Cetnar & Salis, ACS Synth. Biol. 2021, mRNA stability. https://doi.org/10.1021/acssynbio.0c00471
- Cetnar et al., Nature Communications 2024, isoform-aware decay. https://doi.org/10.1038/s41467-024-54059-7
- LaFleur, Hossain & Salis, Nature Communications 2022, Promoter Calculator. https://doi.org/10.1038/s41467-022-32829-5
- Hossain et al., Nature Biotechnology 2020, Nonrepetitive Parts Calculator. https://doi.org/10.1038/s41587-020-0584-2
- Halper, Hossain & Salis, ACS Synth. Biol. 2020, Synthesis Success Calculator. https://doi.org/10.1021/acssynbio.9b00460

### 1.1 Two modes

**Predict / Evaluate.** Concatenate promoter + RBS_i + CDS_i + terminator into one DNA/mRNA molecule. Then:

1. Predict translation initiation rate (TIR) of every annotated CDS, including translational coupling.
2. Predict translation elongation rate (TER) of each CDS.
3. Predict mRNA decay rate.
4. Scan for undesired genetic elements: cryptic promoters, intrinsic and rho-dependent terminators, highly translated internal start codons (HTISC), ribosomal pause sites, repeats ≥ 12 bp, transposon/phage insertion sites, forbidden restriction sites, RNase sites.
5. Emit a GenBank-annotated operon.

**Design.** Search over RBS sequences and synonymous CDS recodings (and optionally intergenic spacers) so the assembled operon simultaneously:

1. Hits target TIRs for each CDS.
2. Uses organism-specific codon tables (Highly Translated or Balanced).
3. Raises mRNA stability.
4. Minimizes internal promoters, internal terminators, HTISC, repeats > 12 bp, ribosomal pause sites, transposon/phage sites, synthesis-hard motifs (homopolymers, extreme GC).

The 2016 SEED talk stated this as 15 models/rules combined with multi-objective genetic algorithm optimization. The public web tool exposes those as selectable design rules and returns several equally optimal designs (a Pareto set), not one winner.

The first CDS is an ordinary RBS-calculator problem. Every CDS after that is not: its TIR is a function of the upstream CDS’s TIR, the intergenic sequence, and RNA structure that overlapping ribosomes unfold. That coupling model is the Operon Calculator’s original core. Everything else is a scan or a recoder that can be swapped independently.

### 1.2 The 15 rules (SEED 2016)

As listed in the 2016 abstract, paraphrased into implementable boxes:

1. Control TIR of each CDS (RBS calculator + coupling).
2. Control TER / synonymous codon choice.
3. Raise mRNA stability (RNase + ribosome protection).
4. Remove internal promoters.
5. Remove internal transcriptional terminators (intrinsic and rho-dependent).
6. Remove highly translated internal start codons.
7. Remove long repeats (published cutoff: 12 bp).
8. Harmonize codon usage to prevent elongation pile-ups.
9. Reduce ribosomal pause sites (unless desired).
10. Avoid transposon insertion sites.
11. Avoid phage att sites.
12. Reduce DNA synthesis complexity (homopolymers, extreme GC).
13. Exclude forbidden restriction sites.
14. Reduce RNase sites (especially 5' UTR AU-ssRNA).
15. Bound system-level RNAP and ribosomal load.

---

## 2. Architecture

```
                    +-----------------------------------------+
                    |           OperonAssembler               |
                    |  promoter + [RBS_i + CDS_i] + terminator|
                    +------------------+----------------------+
                                       | assembled DNA / mRNA
          +---------------+------------+------------+-----------------+
          v               v            v            v                 v
   Translation      mRNAStability   CrypticTX    CrypticTerm     Instability
   (RBSForge +      (RNase +        (Promoter     (intrinsic +    (repeats,
    coupling +       ribosome        Calculator    rho + pause)    IS, att, RE,
    TIR + TER +      protection)     both strands)                 synthesis)
    HTISC)
          |               |            |            |                 |
          +---------------+------------+------------+-----------------+
                                       |
                                       v
                            ObjectiveVector / Scores
                                       |
                    Design mode only:  MultiObjectiveSearch
                    mutates RBS IUPAC windows and synonymous
                    codon choices, re-assembles, re-scores
```

Every box is specified below as its own function with a closed interface. They share a host-organism record and a sequence coordinate system, and nothing else.

### 2.1 Shared types

```
HostPack:                     # RBSForge, do not fork
  name, phylogeny, gram_stain
  t_growth_c
  asd                         # mature 16S 3' tail, RNA alphabet
  s_opt, spacing_push, spacing_pull
  standby: StandbyParams      # v2.0 coeffs present, v1.0 4-nt used
  start_codon_dg              # AUG/GUG/UUG/CUG
  beta                        # 0.45 mol/kcal for ecoli
  footprint_cds               # 13; currently UNUSED by RBSForge
  cutoff_post                 # 35
  cutoff_pre                  # None = whole available 5' UTR, cap 200
  rbs_max_len                 # 35

OperonHost:                   # wrapper, Operon-side
  pack: HostPack
  codon_table
  codon_weights_high          # CAI / tAI / measured
  codon_weights_balanced
  r_elong_nominal_nt_s        # ~60 in E. coli
  ribosome_footprint_nt       # 30 (coupling occupancy; NOT footprint_cds)
  C_unfold, k_P               # coupling; RE-FIT on RBSForge au
  is_motifs[], att_motifs[]
  forbidden_re_sites[]
  rnase_e_like                # True for Gram– (RNase E)

CDS:
  id
  aa_or_nt
  target_tir                  # Design only, au
  rbs_constraint               # IUPAC, optional
  annotated_start              # 0-based on assembled mRNA
  annotated_stop

Operon:
  promoter_dna
  cds_list: [CDS]
  rbs_list: [str]
  terminator_dna
  intergenic_policy           # 'overlap-4' | 'overlap-1' | 'abut' | 'spacer-N' | 'free'
  host: OperonHost
```

Coordinate convention: assembled DNA is 5' to 3' on the coding strand. mRNA is the same sequence with T to U from the promoter’s TSS to the terminator’s 3' end. All scanners report `[start, end)` on that molecule plus strand.

Do not smash `footprint_cds = 13` (30S occupancy on the CDS in the bound-state fold) and `ribosome_footprint_nt = 30` (70S elongating ribosome used in coupling occupancy). They are different physics.

---

## 3. Module A — Operon assembly

Independent of every biophysical model. This is the only place sequence is concatenated.

```
assemble(operon) -> Assembled
    dna: str
    mrna: str                  # T→U from TSS to terminator end
    features: [{id, type, start, end, strand}]
    junctions: [{between, stop_i, start_i+1, d_nt}]
```

Layout:

```
[promoter][5' UTR / RBS_1][CDS_1][intergenic_1][RBS_2][CDS_2] ... [terminator]
```

RBS_1 is a 5' UTR. RBS_{i>1} usually lives inside CDS_{i−1} or in a short intergenic, because that is how coupling is encoded.

### 3.1 Intergenic distance d

Tian & Salis define d from the upstream stop codon to the downstream start codon:

```
d = start_{i+1} − (stop_i_end)
```

Negative means overlap.

| Configuration | Example | d |
|---|---|---|
| Overlap 25 nt | CDS2 start 25 nt upstream of STOP | −25 |
| Canonical 4 nt overlap | AUGA | −4 |
| 1 nt overlap | UAUG | −1 |
| Abutting | UAA AUG | 0 |
| Spacer of N nt | UAA + N nt + AUG | +N |

This single integer is the entire input to the re-initiation coefficient. Compute it once at assembly and store it on `junctions`.

### 3.2 Default junction policy (Design)

- Independent expression (CDS_{i+1} TIR ≈ its own RBS, not a function of CDS_i): overlap ≥ 25 nt (d ≤ −25). Re-initiation drops ~11.6-fold; coupling slope goes to 0.
- Hard-coded ratios via coupling: d = −4 (AUGA) plus a designed inhibitory hairpin that overlaps the upstream CDS. Most sensitive coupling architecture.
- Simple synthetic operons that just need decent downstream TIR: d in {−4, −1, +3} with a strong SD and no inhibitory structure (ΔG_coupling ≈ 0). Re-initiation then contributes a small linear leak (k_reinitiation ≈ 0.007–0.022).

Do not put long unstructured spacers between CDSs unless you also re-run the mRNA-stability module. 5' UTR ssRNA is highly destabilizing; intergenic ssRNA is much less so, but still an RNase E substrate. Prefer 25 nt coding overlap for insulation.

### 3.3 What is transcribed

The mRNA submitted to translation / coupling / stability starts at the TSS and ends at the terminator’s U-tract. Cryptic-promoter and cryptic-terminator scans run on the full DNA, both strands, including the promoter and terminator themselves (so you can see if the intended promoter is the strongest TSS).

---

## 4. RBSForge as the TIR engine

Source of this replica:

- https://github.com/Neladden/RBSForge-branch
- https://github.com/Neladden/RBSForge-branch/pull/1
- Package layout: `rbsforge/` with `rbs/calculator.py`, `hostpack.py`, `thermo/fold.py`, `thermo/duplex.py`, `rbs/standby.py`, `rbs/spacing.py`, `scoring.py`, `docs/MODEL.md`

It reimplements Predict mode of the Salis RBS Calculator from the published v1.0–v2.1 papers plus public v1.0 GPL source, not from the non-public v2.x source.

### 4.1 Master equation

```
ΔG_total = ΔG_standby + ΔG_mRNA:rRNA + ΔG_spacing
         + ΔG_start + ΔG_stacking − ΔG_mRNA

r_proportional = exp(−β · ΔG_total)                 # β = 0.45
r_v1           = 2500 · exp(−ΔG_total / 2.222)      # K = 2500, RT_eff = 1/β
```

| Term | Status in RBSForge |
|---|---|
| ΔG_mRNA | Implemented. MFE of a window around the start. Builtin nested-helix folder, or ViennaRNA RNAfold if on PATH. |
| ΔG_mRNA:rRNA | Implemented. Full duplex search, jointly minimizing hybridization + spacing (v1.0 selection rule). |
| ΔG_spacing | Implemented. v1.0 logistic push / quadratic pull on aligned (Chen/Salis-corrected) spacing. |
| ΔG_start | Implemented. v1.0 table: AUG −1.194, GUG −0.0748, UUG −0.0435, CUG −0.03406 kcal/mol. |
| ΔG_standby | Implemented in v1.0 form only: 4-nt forced-unpaired penalty. v2.0 module geometry is scaffolded on HostPack (`StandbyParams`) and not wired. |
| ΔG_stacking | Not implemented. v2.1 homopolymer-spacer coefficient was never published (Reis & Salis 2020 SI). Reported as explicit 0.0. |
| Ribosome drafting | Not implemented. 2016 kinetic override; needs folding times. |
| CDS-footprint unfolding | Listed in MODEL.md as a v1.1 feature. `footprint_cds = 13` sits on HostPack. Calculator never reads it. |
| Design mode | Entirely absent. |

Temperature: van’t Hoff ΔH/ΔS on Xia et al. 1998 NN stacks. ΔCp neglected. Good to ~60 C, worse above.

β is an empirical fit, not 1/RT. At 37 C, 0.45 mol/kcal is ~3.6× smaller than ideal-solution 1/RT. It was fit almost entirely on E. coli near 37 C. Non-ecoli HostPacks inherit β, spacing curvature, and standby as a documented prior. Rankings inside one HostPack are defensible; absolute TIR is not.

### 4.2 Public API (as it exists)

```python
from rbsforge import predict, RBSCalculator, HostPack, get_hostpack

result = predict(mrna_sequence, species="ecoli", temperature_c=37.0)
# or
calc = RBSCalculator.for_species("ecoli", temperature_c=37.0)
result = calc.predict(mrna_sequence)
```

`PredictResult`:

- `sequence`, `hostpack`, `temperature_c`
- `results: List[StartCodonResult]` — one per AUG/GUG/UUG (CUG is in the energy table, not in the default scan)
- `.best` — highest TIR among non-leaderless
- `.ranked()` — non-leaderless, descending TIR

`StartCodonResult`:

| Field | Meaning |
|---|---|
| `index` | 0-based start in the input mRNA |
| `codon` | AUG / GUG / UUG |
| `leaderless` | True if no SD duplex of length ≥ 4 — then TIR is None |
| `translation_initiation_rate` | exp(−β ΔG_total) |
| `v1_style_rate` | 2500 · exp(−ΔG_total / 2.222) |
| `delta_g_total` | `breakdown.total` |
| `breakdown.terms` | standby, mRNA:rRNA, spacing, start, stacking (=0), mRNA (= −MFE) |
| `duplex.mrna_start/end` | SD span in window coordinates, not full-mRNA |
| `aligned_spacing`, `raw_spacing` | Chen/Salis-corrected spacer |
| `standby_window` | 4-nt forced-unpaired slice in window coords |
| `window_start`, `window_end` | slice of full mRNA that was folded |

Windowing in `_score_start_codon`:

```
pre_len = min(cutoff_pre or available_UTR, available, MAX_PRE_WINDOW=200)
window  = sequence[index - pre_len : index + 3 + cutoff_post]
```

`cutoff_pre = None` means the whole available 5' UTR, capped at 200 nt. For an assembled operon mRNA, CDS_i’s window includes up to 200 nt of upstream sequence, which is enough to see the junction. You do not need CDS1’s start inside CDS3’s window.

Folding: unconstrained MFE of that window. Standby: second fold with 4 nt upstream of the SD forced unpaired. Duplex search is sequence-only (no structure). Footprint 13 is not applied.

Builtin folder (`thermo/fold.py`): nested Watson–Crick/wobble stacks and hairpin loops only. No bulges, internal loops, or multiloops. Conservative (underestimates stability). Fine for “is this isolated RBS strong vs weak.” Not fine for coupling hairpins.

Vienna folder: `RNAfold` wrapper, used automatically if the binary is on PATH (`use_vienna_if_available=True`).

Forced-unpaired already exists on `folder.fold(sequence, forced_unpaired=frozenset of 0-based indices)`. Standby uses it. The calculator does not expose it.

### 4.3 Built-in HostPacks

```
ecoli / e_coli / e_coli_k12     asd=ACCUCCUUA, T=37, s_opt=5, β=0.45  (calibrated)
b_subtilis                      asd=GAUCACCUCCUUUCU (15 nt, RNA-seq mapped), s_opt=8
                                β/spacing/standby = E. coli prior
geobacillus                     B. subtilis aSD borrowed, T=60, uncalibrated
thermus                         aSD placeholder ACCUCCUUA, Ψ flagged not applied, T=70
generic                         E. coli prior throughout
```

For Operon work: use `ecoli` as the calibrated engine. `b_subtilis` is the only other pack whose aSD is real; treat its absolute TIR as uncalibrated. Ignore Thermus until Layer-2 calibration exists.

### 4.4 Adapter: spec `predict_tir` → RBSForge

The Operon TIR oracle the rest of this document calls is:

```
predict_tir(mrna, start, calc) -> TIRResult | None
```

Today:

```python
def predict_tir(mrna, start, calc: RBSCalculator):
    out = calc.predict(mrna)
    hits = [r for r in out.results if r.index == start and not r.leaderless]
    return hits[0] if hits else None
```

Map:

| Operon field | RBSForge |
|---|---|
| tir | Freeze **one** of `v1_style_rate` or `translation_initiation_rate`. See §4.5. |
| dG_total | breakdown.total |
| dG_terms | breakdown.terms / notes |
| sd_span | (window_start + duplex.mrna_start, window_start + duplex.mrna_end) |
| standby_span | (window_start + standby_window[0], window_start + standby_window[1]) |
| footprint_span | Not computed. After §6, [start, start+footprint_cds). |
| inhibitory_helices | Not returned. Folder.pairs exist internally and are discarded. Use constrained fold instead of helix split. |
| design_rbs | Absent. Build per §5. |

HTISC is a filter over `result.results`. Coupling consumes two scores (folded / unfolded). Inverse design is a search around this oracle.

### 4.5 TIR scale (load-bearing choice)

Worked example: ΔG_total = −10 kcal/mol

- proportional ≈ 90
- v1-style ≈ 2.2 × 10^5

The Operon docs and Tian 2015 talk in 1 … 100,000+ au. C = 0.81 and k_P = 10 were fit on that scale.

Do not copy C = 0.81 or k_P = 10 onto `translation_initiation_rate`. Occupancy min(1, C · r) saturates for almost every real RBS if r is the v1-style 10^5 number, and never saturates if r is the ~90 proportional number.

Until C, k_P are re-fit on RBSForge outputs against a bi-cistronic reporter panel (Tian’s 22 + 76 constructs are the template):

- Freeze **v1_style_rate** as the au shown to users and used in HTISC thresholds (τ ≈ 1000 au on the Salis scale).
- Treat published C, k_P as order-of-magnitude priors.
- k_reinitiation(d) is a **ratio** of rates and does transfer (0.022 at d = −4, 0.0072 for d in [0, 25], reverse-overlap factor 11.6 at 25 nt).

Put the conversion in one function so HTISC, coupling, and Design never mix units.

### 4.6 What already works without new RBS physics

- Module F HTISC: `calc.predict(full_mrna)`, drop annotated starts, drop leaderless, threshold the rest.
- Module B first CDS, and insulated downstream CDS (d ≤ −25, no overlapping hairpin): one predict on assembled mRNA, pick the annotated index.
- Module G internal-SD pauses (partial): `thermo.duplex.best_hybridization` is a sliding anti-SD search. Not wrapped, engine is there.
- Modules A, C, E, H, I, J: no TIR calls.
- 1-cistron Predict: RBSForge plus scanners. No coupling.

### 4.7 What RBSForge will not do

No mRNA-stability sidecar, no promoter prediction, no terminator model, no codon recoding, no Design mode. Those were never this package’s job. MODEL.md is explicit.

Biology the original Salis model also does not capture: IF dynamics, 50S joining, tmRNA, CsrA/Hfq, pseudoknots, G-quadruplexes. Leaderless / SD-independent initiation is flagged and excluded from ranked results rather than mis-scored.

---

## 5. Missing piece: `design_rbs`

RBSForge is Predict-only. The original inverse designer is public:

- https://github.com/hsalis/Ribosome-Binding-Site-Calculator-v1.0
- https://github.com/hsalis/Ribosome-Binding-Site-Calculator-v1.0/blob/master/RBS_MC_Design.py
- Salis, Mirsky & Voigt, Nature Biotechnology 27:946 (2009), Fig. 2D. https://doi.org/10.1038/nbt.1568
- Salis, Methods in Enzymology 498:19 (2011) §4.4, including IUPAC constraints. https://doi.org/10.1016/B978-0-12-385120-8.00002-4

The Operon Calculator’s Design mode, and the 2011 chapter, are that simulated annealer plus IUPAC / constant-upstream. You do not need Salis’s private v2.x designer. You need this loop, calling RBSForge as the energy oracle.

### 5.1 Physics inversion

TIR and ΔG_total are 1–1 under the Boltzmann map. Design is “hit a target ΔG,” not “hit a TIR.”

Using the frozen v1 map:

```
ΔG_target = RT_eff * (logK − log(TIR_target))
          = 2.222 * (7.824 − ln(TIR_target))
```

This is exactly v1.0:

```
dG_target = RBS_Calculator.RT_eff * (RBS_Calculator.logK - math.log(float(TIR_target)))
```

Objective:

```
O = |ΔG_total − ΔG_target|
stop when O ≤ 0.25 kcal/mol
```

0.25 kcal/mol is about 1.13× in TIR at β = 0.45. Do not tighten it until Vienna is the folder; builtin MFE noise is larger than that.

### 5.2 Interface

```
design_rbs(
  cds_nt,                 # from start codon, at least cutoff_post (35) nt
  target_tir,             # au, same scale as predict
  hostpack,
  constraint_iupac=None,  # same length as designed RBS, or None
  constant_upstream=None, # pre_seq; NOT mutated
                          # cistron 1: promoter ITR / 5' UTR prefix
                          # cistron i>1: 3' tail of CDS_{i-1} including the junction
  rbs_init=None,
  rbs_len_range=(20, 35), # HostPack.rbs_max_len = 35
  max_iter=10000,
  temperature_c=None,
) -> DesignResult(rbs, predicted_tir, dG_total, breakdown, n_evals, hit_tol)
```

Full scored mRNA every evaluation:

```
mRNA = constant_upstream + RBS + cds_nt
start = len(constant_upstream) + len(RBS)
hit   = predict_tir(mRNA, start)
```

If that start is leaderless, treat as infeasible (infinite energy), not as TIR = 0.

### 5.3 Search (copy v1.0, add IUPAC)

From RBS_MC_Design.py / Salis 2009 / Salis 2011:

| Piece | v1.0 value | Keep? |
|---|---|---|
| Moves | replace 0.80, insert 0.10, delete 0.10 | Keep replace. Disable insert/delete if the IUPAC mask has fixed length (Operon constraints are fixed-length degenerate strings). |
| Init | random 35-mer, SD-biased when target is strong | Keep. Plant a core complementary to hostpack.asd (E. coli ACCUCCUUA → UAAGGAGGU-like) at aligned spacing s_opt (5 E. coli, 8 B. subtilis); randomize the rest. |
| Metropolis | accept if ΔO < 0; else P = exp(−ΔO / RT) | Keep |
| RT_init | 0.6 kcal/mol | Keep |
| Anneal | every 50 moves: if accept ratio > 0.20, RT /= 2; if < 0.01, RT *= 2 | Keep |
| MaxIter | 10,000 | Default. Operon inner loop should cap lower (~2,000) because it will call this many times. |
| tol | 0.25 kcal/mol | Keep |
| Max RBS length | cutoff = 35 | HostPack.rbs_max_len |
| Reject new start codon in RBS | RemoveStartCodons (ATG/GTG/TTG) | Keep |
| kinetic_score ≤ 0.50, three_state_indicator ≤ 6 | v1.0 folding-kinetics guards | Skip. RBSForge does not compute them. Do not fake them. |
| helical loop 4–12 | v1.0 | Optional. Vienna already penalizes pathological folds via ΔG_mRNA. |

IUPAC (not in the public Monte_Carlo_Design; is in the 2011 chapter and the Operon web tool):

- Mask length = designed RBS length.
- Replace: new base ∈ allowed set at that position (N → ATCG, R → AG, …).
- Literal A/C/G/T positions are frozen.
- If mask is absent: all 35 positions are N; length may vary with insert/delete.

`constant_upstream` is outside the mutatable RBS. For overlapping cistrons it is the thing coupling cares about; the designer must not nibble it.

### 5.4 Initialization details (v1.0 GetInitialRBS)

v1.0 maps ΔG_target onto [0, 1] using dG_range_high = 25.0 and dG_range_low = −18.0 kcal/mol, then sets P(choose SD nucleotides), core length, and allowed spacing deviation. After drawing a random RBS it strips start codons. Keep the spirit: strong targets get a longer, better-placed SD; weak targets get a worse SD or worse spacing. Then let Metropolis move.

max_init_energy = 10.0 kcal/mol: if the random init is worse than that, redraw.

### 5.5 Hard constraints on a candidate (reject, do not Metropolis)

- IUPAC satisfied
- No start codon in the RBS
- Intended start is not leaderless
- Length in rbs_len_range and ≤ rbs_max_len
- Optional Operon filters here: no forbidden RE site in upstream+RBS+N-term, no 12-mer clash with used_kmers — or leave those to the outer Operon GA

### 5.6 What this is not

- Not coupling-aware. design_rbs is mono-cistronic. For CDS_{i>1} always re-score with Module B after assembly. Pass constant_upstream = CDS_{i−1} tail so the folded window at least sees the junction; do not expect the returned TIR to be the coupled TIR.
- Not a Pareto search. One target ΔG. Operon-level tradeoffs (repeats, HTISC, promoters) stay in the outer NSGA-II.
- Not the RBS Library Calculator (Farasat et al., Mol. Syst. Biol. 2014, https://doi.org/10.15252/msb.20134955). That is degenerate-library design; a later expansion.

### 5.7 Minimal version to ship first

1. Fixed-length RBS, replace-only, optional IUPAC.
2. SD-biased random init at s_opt.
3. O = |ΔG_total − ΔG_target| from RBSForge predict on upstream+RBS+cds.
4. Metropolis + the RT schedule above.
5. Return best-so-far even if tol is missed, with hit_tol: bool.

That is enough for Operon Design as an inner proposal. Do not port kinetic_score / helical-loop repair from v1.0 until Predict is v1.1-complete (footprint + Vienna).

### 5.8 Architecture consequence

The earlier “Stage A = call design_rbs, Stage B = couple” split is still useful as a **proposal**, not as truth. Operon Design is sample-and-score from the start: mutate RBS nucleotides inside the IUPAC mask, assemble, predict, couple, scan. design_rbs is a cheap initializer for CDS1 and for insulated CDS_{i>1}. A fallback initializer if you have not built design_rbs yet: plant AGGAGG at aligned spacing s_opt and let the outer GA move it. That replaces a large fraction of what people use Design mode for.

---

## 6. Missing piece: CDS-footprint unfolding (`footprint_cds = 13` unused)

### 6.1 What the constant was supposed to mean

HostPack.footprint_cds = 13 is the v1.1/v2 CDS occupancy: once 30S is bound, about 13 nt past the start codon cannot form RNA structure.

- Salis 2011 (Methods Enzymol.): final-state spacing minimization and CDS-footprint unfolding. https://doi.org/10.1016/B978-0-12-385120-8.00002-4
- Espah Borujeni et al., Nucleic Acids Research 2017: later footprint definition. https://doi.org/10.1093/nar/gkx262
- RBSForge docs/MODEL.md lists this as a v1.1 feature of the model being reimplemented.

v1.0 source (RBS_Calculator.py) did something cruder: cutoff = 35 window, and a huge footprint = 1000 meaning “ignore almost all post-start structure in the bound state.”

https://github.com/hsalis/Ribosome-Binding-Site-Calculator-v1.0/blob/master/RBS_Calculator.py

RBSForge today:

- Folds [start − pre, start + 3 + cutoff_post) unconstrained (cutoff_post = 35).
- Puts −MFE into ΔG_total as the entire unfolding penalty.
- Uses forced-unpaired only for the 4-nt v1.0 standby slice.
- Never reads footprint_cds.

So structure in amino acids 1–5 is penalized as if the ribosome had to unfold it and also as if it still existed in the bound state. N-terminal recoding looks TIR-neutral when it is not. Coupling’s “inhibitory helix overlapping the footprint” is not a real term.

### 6.2 Correct two-fold (v1.1)

For each start:

```
window = mrna[window_start : window_end]          # already implemented

# initial state: empty 30S, full structure
E_initial = folder.fold(window).delta_g           # today's delta_g_mrna

# final / bound state: ribosome occupies SD + spacer + start + footprint
occupied = SD_span ∪ spacer ∪ start_codon ∪ [start, start + footprint_cds)
         ∪ standby_4nt     # see double-count note
E_bound = folder.fold(window, forced_unpaired=occupied_local).delta_g

unfolding = E_initial − E_bound                   # typically ≥ 0
```

Then hybridization, spacing, and start score the bound complex.

Standby’s own 4-nt forced-unpaired is a subset of occupied. Either drop the separate v1.0 standby term once occupied includes those 4 nt, or keep it and do not double-count them.

Recommended until v2.0 standby is wired: occupied = start codon + footprint_cds; keep the existing 4-nt standby term for the nts immediately upstream of the SD.

forced_unpaired already exists on BuiltinFolder.fold. For Vienna, RNAfold --constraint with `x` at occupied positions. This is the same hook coupling needs. Implement **one** `predict_one(mrna, start, extra_unpaired=None)`:

- default occupied = footprint 13 + start codon
- coupling passes additional unpaired = upstream-CDS nucleotides that fall in the window

Do not implement two parallel constrain systems.

### 6.3 Sign / bookkeeping

Today:

```
breakdown.add("mRNA", -delta_g_mrna)   # MFE is negative for stable folds → positive penalty
```

After the two-fold:

```
unfolding = E_initial - E_bound_constrained    # typically ≥ 0
breakdown.add("mRNA", unfolding,
              "unfolding cost, footprint+start constrained in bound state")
```

If you keep adding −MFE_unconstrained you have not implemented the footprint; you have renamed it.

### 6.4 What this changes downstream

- First ~5 codons of every CDS become TIR-sensitive. Operon recoding must re-predict after N-terminal synonym swaps. Freezing the first 15 amino acids after RBS design is the faster, slightly wrong option; re-running predict_tir is correct.
- HTISC inside a CDS: the footprint of a cryptic start is 13 nt of whatever protein it would make; constrained fold still applies.
- Coupling ΔG_coupling helices that sit on the upstream stop / first footprint of the downstream start are now actually unfoldable in the bound state of CDS_i. The two-call mixture in Module B stays; the inner predict_one just gets more accurate.
- Test: a hairpin entirely inside nt +4…+13 of the CDS should hurt TIR less after this change than before (ribosome pays to open it, then sits on it). A hairpin entirely in the SD should still hurt.

### 6.5 Shared `predict_one` (the one RBSForge patch Operon actually needs)

```
predict_one(mrna, start, extra_unpaired_global=None) -> StartCodonResult
    window = mrna[window_start:window_end]
    unpaired_local = {i - window_start for i in extra_unpaired_global
                      if window_start <= i < window_end}
    unpaired_local |= footprint and start-codon indices in window coords
    ΔG_mRNA from constrained vs unconstrained fold as in §6.2
    duplex / spacing / start = same as today (sequence, not structure)
    standby = recompute on the constrained fold
    rest of ScoreBreakdown = unchanged
```

Coupling, footprint, and standby then share one code path so the constrained folds cannot drift.

---

## 7. Module B — Translational coupling

This is the only translation physics that is not already in RBSForge. Source: Tian & Salis, NAR 2015, Eqs. 1–5. https://doi.org/10.1093/nar/gkv635  PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC4538824/

Five fitted parameters. Nested iteration from 5' to 3'.

### 7.1 First CDS (no coupling)

```
r_1 = predict_tir(mrna, start_1).tir
```

β and K live inside RBSForge. Do not refit them here.

### 7.2 Two sources of downstream translation

```
r_i = r_reinit_i + r_denovo_i
```

#### Re-initiation

```
r_reinit_i = k_P * k_reinitiation(d_{i-1,i}) * r_{i-1}
```

Reported values (Tian 2015):

| Quantity | Value | Notes |
|---|---|---|
| k_P | 10 (best-fit 14; 95% CI [8, 33]) | Converts the ratiometric coefficient into TIR au. Treat as 10 on **their** au. Re-fit on RBSForge. |
| k_reinitiation(d) at d = −4 | 0.022 | AUGA: tRNAfMet anticodon 3′-UACU-5′ pairs the 4 nt, no scanning. Transfers. |
| k_reinitiation(d) for d in [0, 25] | 0.0072 ± 0.0018 | Forward scan, almost flat. Transfers. |
| k_reinitiation at d = −25 | ~11.6-fold below the d = 0 value (~0.00062) | Reverse scan collides with elongating 70S. Transfers. |
| d in (−25, 0) except −4 | linear interpolation | Three-parameter piecewise fit, R² = 0.91. |
| d > 25 | ~0.5–0.7 per +100 nt | From Levin-Karp et al.; not refit by Tian. |

```python
def k_reinitiation(d: int) -> float:
    if d == -4:
        return 0.022
    if 0 <= d <= 25:
        return 0.0072
    if d > 25:
        return 0.0072 * (0.6 ** ((d - 25) / 100.0))
    if d <= -25:
        return 0.0072 / 11.6
    return 0.0072 / 11.6 + (0.0072 - 0.0072 / 11.6) * (d + 25) / 25.0
```

No Shine–Dalgarno is required for re-initiation. An SD only raises the de novo (basal) rate.

#### De novo initiation with ribosome-assisted unfolding

An elongating ribosome occupies the junction with probability

```
f_unfold = min(1, C * r_{i-1}^{phys})
```

Theoretically C = footprint / r_elong ≈ 30 nt / 60 nt s⁻¹ = 0.50 s when r is in s⁻¹. Empirically Tian fitted C = 0.81 ± 0.17 (95% CI [0.67, 1]) on the au scale they used. The conversion between RBS-calculator au and s⁻¹ is absorbed into C. If you reuse 0.81 you must use their au. On RBSForge, re-fit C.

The de novo rate is a two-state mixture of the fully folded and fully unfolded RBS-calculator evaluations:

```
r_folded   = predict_one(mrna, start_i).tir
r_unfolded = predict_one(mrna, start_i,
                 extra_unpaired=nts_of_CDS_{i-1} that fall in this window).tir

r_denovo_i = (1 - f_unfold) * r_folded + f_unfold * r_unfolded
```

which is sigmoidal in log r_{i−1}. At low upstream TIR the downstream CDS sees the inhibitory hairpin; at high upstream TIR the hairpin is held open and downstream TIR saturates at the unfolded (basal) rate.

A helix is coupling-sensitive if it overlaps the upstream CDS, including the last nucleotide of the stop. Those are the nucleotides to force unpaired. You do not need a helix list from RBSForge; the constrained fold is the split.

Do not interpolate a single ΔG_total in log-space unless you have verified it against the two-state mixture. The mixture is what the paper’s sigmoid is.

If the downstream start is leaderless (no SD), r_folded is undefined. Treat de novo as 0 and keep only re-initiation. That is a real operon architecture (SD-less overlapping starts); RBSForge is correctly refusing to invent an SD score.

### 7.3 Nested iteration for n cistrons

```python
def coupled_tirs(assembled, host, calc) -> list[float]:
    n = len(assembled.starts)
    r = [None] * n
    r[0] = predict_one(assembled.mrna, assembled.starts[0], calc=calc).tir
    for i in range(1, n):
        d = assembled.junctions[i-1].d_nt
        unpaired = range(assembled.cds_end[i-1] - overlap_into_window,
                         assembled.cds_end[i-1])
        tir_fold   = predict_one(assembled.mrna, assembled.starts[i], calc=calc)
        tir_unf    = predict_one(assembled.mrna, assembled.starts[i],
                                 extra_unpaired=unpaired, calc=calc)
        f = min(1.0, host.C_unfold * to_physical(r[i-1], host))
        r_denovo = (1-f)*tir_fold.tir + f*tir_unf.tir
        r_reinit = host.k_P * k_reinitiation(d) * r[i-1]
        r[i] = r_reinit + r_denovo
    return r
```

Each step uses the already-coupled upstream rate, not the mono-cistronic prediction. That is why changing the promoter (or RBS_1) can make a downstream RBS “disappear” in Evaluate mode.

ViennaRNA is required when this module is on. Builtin folder cannot represent bulge-containing intergenic hairpins. Refuse or warn if folder.backend != "vienna".

### 7.4 What coupling does not include

- Transcription–translation coupling (~40% more protein per kb of distance from the TSS). Residual error, not in the model.
- Ribosome crowding that would decrease downstream TIR at very high r_{i−1}. They looked for it and did not see it.
- mRNA-level polarity (RNAP fall-off, RNase entry at intergenic sites). Module D, not this module.
- Elongation-limited translation. The paper assumes CDSs are codon-optimized enough that initiation is rate-limiting. If you skip codon optimization, TER can become the bottleneck and coupling will over-estimate protein.

### 7.5 Design use of coupling

Three regimes, from Tian’s design-criteria section:

1. ΔG_coupling = 0. r_i vs r_{i−1} is a straight line with slope k_P * k_reinitiation(d) ∈ [0, 0.022]. Small slope, large effect when basal rates differ a lot. Example: 100-fold basal difference at d = −4 collapses to 31-fold once re-initiation is included. Kill this by overlapping ≥ 25 nt.
2. Inhibitory overlapping hairpin. Ratio r_{i−1}/r_i is set by ΔG_coupling and the basal unfolded rate. At basal −20 kcal/mol, varying ΔG_coupling from −5 to −30 kcal/mol tunes the ratio from ~2-fold to ~70-fold. Weak hairpin + high basal: ~10-fold swing as upstream TIR is titrated. In principle 1000-fold if the hairpin is very stable and fully unfolded.
3. Do not try to encode ratios with re-initiation alone. The dynamic range is tiny. Encode ratios with SD strength (de novo, unfolded) and use coupling either as insulation (d ≤ −25) or as a hairpin-gated dependency.

A Design-mode sub-search for a junction: enumerate short IUPAC intergenic sequences subject to d and a target ΔG_coupling, then re-score with this module.

Do not invert coupling. There is no closed form.

---

## 8. Module C — Translation elongation rate and codon recoding

Independent of coupling. Operates on one CDS at a time. Used as a score (Predict) and as a mutator (Design).

### 8.1 TER prediction

The web tool reports TER as an average synthesis rate in nucleotides per second.

```
r_elong_codon[c] = r_max * w[c] / w_max          # w from tAI or measured dwell times
TER_cds = 3 / mean(1/r_elong_codon[c] for c in cds)   # nt/s, harmonic mean
```

Use a harmonic mean, not an arithmetic mean: one slow codon dominates dwell time. Typical E. coli: r_max ≈ 20 aa/s ≈ 60 nt/s. Pause sites (Module G) are extra dwells added to the mean.

Next increment: a TASEP / ribosome-flow model with initiation rate r_i from Module B, codon-specific hops, and footprint 10 codons. That is how you detect pile-ups. The 2016 abstract lists “harmonize synonymous codon usage to prevent translation elongation pile-ups” as a rule.

```
density[j] ≈ r_init / r_elong[j]     # if this approaches 1, the ribosome jams
```

When r_init is small, rare codons are harmless. When r_init is large, they are not. That is the entire reason for two codon tables.

### 8.2 Two codon tables (Design mode input)

**Highly Translated Synonymous Codons.** For every amino acid, always pick the codon with the highest w[c]. Use when a small number of CDSs will be overexpressed.

**Balanced Synonymous Codons.** Match local elongation to the CDS’s target TIR:

```python
def choose_codon(aa, target_tir, table_high, table_genomic):
    mix = clip(log10(target_tir) / log10(tir_max), 0, 1)
    p = softmax(mix * log(table_high[aa]) + (1-mix) * log(table_genomic[aa]))
    return sample(p)   # or argmax for deterministic recoding
```

Use when many CDSs are expressed at varied levels. Fast initiation + slow elongation = pile-up, aborted elongation, and a TER that no longer matches the TIR you designed.

Host-specific tables: codon usage of highly expressed genes (ribosomal proteins) for table_high; whole-genome usage for table_genomic. Do not use E. coli tables in Bacillus.

### 8.3 Recoding as a mutator, not a one-shot rewrite

A single greedy CAI rewrite will fight every other objective (repeats, cryptic promoters, internal starts, RNase sites, synthesis). Recoding is a neighborhood operator inside the genetic algorithm:

- Pick a random codon whose amino acid has synonyms.
- Swap to another synonym (biased toward the active table).
- Re-score the whole operon.

Preserve: start codon, stop codon, user-locked motifs, and — until footprint unfolding is on — the 5'-most ~10–15 codons, because they change TIR. After footprint is wired, re-run predict_tir after every N-terminal recode. That is the correct option.

### 8.4 Internal start-codon suppression is a recoding objective

Every in-frame and out-of-frame AUG/GUG/UUG inside a CDS is a potential HTISC. Recoding can:

- Remove the start codon (AUG → AUC, etc., if the amino acid allows).
- Destroy its SD (synonym swaps in the 10–20 nt upstream, possibly in another frame).
- Leave it but drive its TIR below a threshold (~1000 au on the v1 scale).

HTISC count is a separate objective, not folded into TIR error.

---

## 9. Module D — mRNA stability

Independent of coupling except that it consumes the coupled TIRs as ribosome-protection inputs. Two generations of Salis models exist. The Operon Calculator docs describe the 2021 biophysical version. The 2024 calculator is a hybrid LightGBM on top of those features plus isoforms. Replicate 2021 first.

RBSForge explicitly does not include this.

### 9.1 Physics the 2021 paper measured

Cetnar & Salis, ACS Synth. Biol. 2021. https://doi.org/10.1021/acssynbio.0c00471
bioRxiv: https://www.biorxiv.org/content/10.1101/2020.07.22.216051v1

82 operons, RT-qPCR, exponential E. coli:

| Perturbation | Effect on mRNA level |
|---|---|
| Long ssRNA in the 5' UTR | up to 9.4-fold down |
| Lowering translation rate | up to 11.8-fold down |
| RNase sites in intergenic regions | much smaller |
| Terminator efficiency | no effect on upstream mRNA |
| Long ssRNA in the 3' UTR | no effect on upstream mRNA |

A polarity-from-the-3'-end model is the wrong expansion. Decay in E. coli for these synthetic operons is 5'-end / ribosome-protection dominated. Do not put terminator strength into the half-life formula for this host and architecture.

### 9.2 2021 biophysical model (implement first)

**Ribosome protection.** Average unprotected distance between ribosomes:

```
Δx = r_elong / r_init − L_footprint     # nt of naked mRNA between 70S’s
# r_init in s⁻¹, r_elong in nt/s, L_footprint = 30
```

Convert TIR au → s⁻¹ with the same factor you used for C in Module B, or treat Δx as a feature and fit one linear coefficient on log-mRNA.

More protection (smaller Δx) → slower decay. At low TIR, Δx is large, RNase E has a long internal substrate, decay is fast. This is why “design a weak RBS” also silently designs an unstable mRNA.

For multi-cistronic operons use each CDS’s coupled TIR and a length-weighted average Δx, or the minimum protection (the least-translated cistron is the RNase entry point). Weight toward the 5' cistron.

**RNase binding probability.** RNase E (E. coli) cuts RNA that is single-stranded, AU-rich, and accessible.

```
fold the 5' UTR + ~30 nt of CDS1
ss[i] = 1 if unpaired in MFE (or 1 − pairing probability)
motif_w['A']/['U'] > ['G']/['C']
RNase_score = sum(ss[i] * motif_w[nt[i]]) over 5' UTR
```

Long unpaired AU tracts are the 9.4-fold effect. G/C hairpins at the 5' end are protective (they also block RppH in the 2024 model).

**Combine (minimal):**

```
log(mRNA) = a0 + a1 * (−RNase_score) + a2 * (−Δx) + a3 * 1[5p_hairpin]
k_decay   = k0 * exp(b1 * RNase_score + b2 * Δx)
t_half    = ln(2) / k_decay
```

Fit {a, b} on any paired sequence/qPCR set; the 2021 paper is the calibration.

### 9.3 2024 expansion

Cetnar, Hossain, Vezeau & Salis, Nat. Commun. 15:9601 (2024). https://doi.org/10.1038/s41467-024-54059-7

Public code: https://github.com/hsalis/SalisLabCode (mRNA_Stability_Calculator)

62,120 barcoded 5' UTRs, rifampicin decay, half-lives 0.31–25.4 min.

Feature pipeline:

1. Isoforms. Run the Promoter Calculator (Module E) on the promoter+UTR. Keep the top 5 TSS isoforms.
2. Per isoform: rppH (first 4 nt), TX_rate, TL_rate, initial_sum_external_ssRNA, initial_dsRNA_sum, initial_external_ssRNA_{A,T,G,C}, GQUAD, iMOTIF, CsrA, nt_0..nt_79.
3. LightGBM, one model per timepoint (T0, T2, T4, T8, T16), target log(RNA_counts/DNA_counts).
4. Test R²: 0.69 (T0) falling to 0.43 (T16).

Hyperparameters from the public training script: num_leaves=100, max_depth=5, n_estimators=119, learning_rate=0.1, bagging_fraction=0.50, feature_fraction=0.25, min_data_in_leaf=50.

Design rules the 2024 model recovered, usable even without the GBDT:

- Stable RppH 5' end (not a fast-dephosphorylation tetramer).
- A 5' hairpin with a short unpaired 5' extension.
- High TIR on CDS1 (protection), unless you also introduce long ssRNA — at high TIR, extra ssRNA increases decay (R = 0.92 in that slice). At low TIR the mRNA is already unstable and extra ssRNA does little.
- Avoid long polyA loops in the 5' UTR.
- G-quadruplexes help; i-motifs do not.
- Do not waste degrees of freedom on terminator strength or 3' UTR ssRNA for this host.

### 9.4 Design-mode actions that raise mRNA stability

1. Recode / redesign RBS_1 so CDS1 TIR is not accidentally tiny.
2. Insert a small 5' UTR hairpin (6–8 bp stem, 4–6 nt loop) that does not occlude the standby site.
3. Shorten unpaired AU tracts in the 5' UTR.
4. Pick a TSS-proximal tetramer with a slow RppH rate (once you have the 256-tetramer table from the 2024 SI).
5. Do not try to stabilize by massaging the terminator.

---

## 10. Module E — Internal promoters (cryptic transcription)

Independent scanner. Source: LaFleur, Hossain & Salis, Nat. Commun. 13:5159 (2022). https://doi.org/10.1038/s41467-022-32829-5

Public code: https://github.com/hsalis/SalisLabCode (Promoter_Calculator)

This is how Design mode removes internal promoters and how Predict mode lists them. It is also how the 2024 mRNA model finds isoforms.

### 10.1 State of a σ70 promoter

For every candidate TSS:

| Element | Length | Notes |
|---|---|---|
| UP | 24 bp | distal 12 + proximal 12; groove width + rigidity |
| 1 nt gap | 1 | between UP and −35 |
| −35 hexamer | 6 | split into two 3-mers, energy matrix |
| spacer | 15–20 bp | quadratic penalty in length; last 2 nt are the extended −10 (TGN) |
| −10 hexamer | 6 | split into two 3-mers |
| discriminator | 6–10 bp | first 3-mer scored |
| ITR (R-loop) | 20 bp | DNA:RNA hybrid energy of the initial transcribed region |

Minimum promoter ~78 bp, maximum ~89 bp. Scan both strands. For circular plasmids, wrap the origin.

### 10.2 Free energy and transcription rate

```
dG_10     = matrix10_left[hex10[0:3]] + matrix10_right[hex10[3:6]]
dG_35     = matrix35_left[hex35[0:3]] + matrix35_right[hex35[3:6]]
dG_ext10  = matrix2mer[spacer[-3:-1]]
dG_disc   = matrix3mer[disc[0:3]]
dG_spacer = 0.1463*s**2 − 4.9113*s + 41.119    # s = spacer length
dG_ITR    = coef_ITR * (DNA:RNA_hybrid / 4.30)
dG_UP     = coef_dist*(groove_dist/256) + coef_prox*(groove_prox/255)
            + coef_rig*(rigidity/25.78)
dG_total  = dG_10 + dG_35 + dG_disc + dG_ITR + dG_ext10 + dG_spacer + dG_UP + intercept

Tx_rate   = K * exp(−β * dG_total)
```

From the public implementation:

- K = 42 for E. coli MG1655 and in vitro.
- β_vivo = 1.636217, β_vitro = 0.816326.
- LOGK = −2.80271176 appears in the apparent-dG conversion used during training.

Coefficients live in free_energy_coeffs.npy (346 parameters). Do not refit unless you have a new MPRA.

### 10.3 Scanning algorithm (Predict)

```
for strand in (+1, −1):
    for each TSS:
        for disc_len in 6..10:
            for spacer_len in 15..20:
                score configuration
                keep min dG_total per TSS
```

A site is internal / cryptic if its TSS is not the intended TSS and Tx_rate > τ_tx (start with a few percent of the intended promoter’s Tx_rate, or an absolute cutoff around the weak-constitutive range).

Report count, positions, strands, and Tx_rates. Design mode’s objective is n_internal_promoters, optionally weighted by log(Tx_rate).

### 10.4 Design-mode kill

Cryptic promoters inside CDSs: synonymous swaps that break the −10 or −35 hexamer. Inside RBSs / spacers there is more nucleotide freedom — one reason RBS constraints are IUPAC rather than a fixed sequence.

Cheap prefilter: PWM-scan for TATAAT-like and TTGACA-like hexamers at 15–20 bp spacing, then rescore hits with the full model. The full model should gate the objective, because non-canonical motifs collectively matter.

σ70 only. Other sigma factors and T7 RNAP are different models (the lab has a T7 Promoter Calculator in SalisLabCode). Expanding = swapping this module, not changing coupling.

---

## 11. Module F — Highly translated internal start codons (HTISC)

This is RBSForge run as a scanner, plus a threshold. Truncated proteins are a common operon-failure mode (2016 abstract).

```python
def scan_htisc(mrna, annotated_starts, calc, tau_htisc=1000.0):
    hits = []
    out = calc.predict(mrna)
    for r in out.results:
        if r.leaderless:
            continue
        if r.index in annotated_starts:
            continue
        tir = r.v1_style_rate   # same au as everywhere else
        if tir >= tau_htisc:
            hits.append(r)
    return hits
```

The “Translated Open Reading Frames” plot in the web tool is exactly this: TIR vs position for every start, with annotated CDSs highlighted. Check that no internal start outranks the annotated one.

Design-mode kill: synonym swaps that remove the start or its SD. Out-of-frame AUGs that land on a codon you cannot change without changing the protein are killed by destroying the SD, or by introducing a nearby stop in that frame (only legal if it doesn’t create a pause / RE / promoter).

tau_htisc is a knob. 1000 au is a reasonable default on the v1 scale. For toxic truncated products, lower it.

Leaderless starts have tir=None. They are not HTISC.

---

## 12. Module G — Ribosomal pause sites

Independent motif + structure scanner inside CDSs. The Operon Calculator lists “fewer ribosomal pause sites” as a selectable objective. The lab never published a standalone Pause Calculator. The operational definition used in related Salis design tools is the union of:

1. Internal Shine–Dalgarno / anti-SD hybridization during elongation. Hexamers complementary to the 16S tail (E. coli 3′-AUCACCUCCUUA-5′, motifs like AGGAGG, GGAG, GAGG and weaker neighbors) occurring in-frame or out-of-frame inside a CDS. Score as n_internal_SD weighted by predicted SD duplex ΔG. Reuse RBSForge’s `thermo.duplex.best_hybridization` on a sliding window. Host-dependent: use HostPack.asd, not a hardcoded E. coli tail.
2. Slow codon runs. Consecutive low-w codons (Pro-Pro, Pro-Glu, Arg-Arg rare, Ile AUA, Leu CUA). A run of ≥2–3 is a pause.
3. Polyproline / stall motifs. PPP, PPG. SecM / TnaC class if you care about arrest peptides. For generic operons, PPP is enough.
4. Stable mRNA hairpins in the ribosomal E/P site window (optional, weaker evidence for bacteria than the internal SD).

Predict: count and positions. Design: synonym swaps that break the SD-like hexamer or replace the slow codon. The 2016 abstract says you may want to keep a pause (e.g. to allow folding). That is a per-CDS whitelist.

---

## 13. Module H — Internal transcriptional terminators

There is no Salis Terminator Calculator in SalisLabCode (the repo has Promoter, T7 Promoter, mRNA Stability, NRP, Synthesis Success, ELSA, Oligopool, ModelTestSystem). Terminators in the Operon Calculator are a scan plus a parts-database append, not a fitted rate model of their own.

https://github.com/hsalis/SalisLabCode

### 13.1 What the Operon product actually does

From the De Novo DNA Operon page:

- Input: optional terminator name + sequence from their genetic parts database, placed at the 3' end.
- Predict output: “internal transcriptional terminators” listed under Undesired Genetic Elements.
- Design objective 5: fewer internal transcriptional terminators.
- Design extras: “Non-repetitive promoters and transcriptional terminators may be inserted at the beginning and end of the operon.”
- mRNA stability: Cetnar & Salis 2021 found changing terminator efficiency or 3' UTR ssRNA had **no effect on upstream mRNA** in their 82-operon E. coli set. The terminator is not a Module D lever. It is a “don’t die in the middle of the operon” scan plus a 3' part.

2016 SEED abstract, rule (v): remove internal terminators, both intrinsic and rho-dependent.

### 13.2 Intrinsic (rho-independent)

Motif: GC-rich hairpin + U-tract. Two layers.

**Presence (TransTermHP / Kingsford 2007 style)** — enough for the Operon count.

Kingsford, Ayanbule, Salzberg, Genome Biology 2007. https://doi.org/10.1186/gb-2007-8-2-r22

- Stem 4–15 bp, loop 3–8, T-tail ~15 nt
- Score = hairpin MFE + U-tract RNA:DNA hybrid energy
- Scan both strands
- Hit is “internal” if it is not the user-supplied 3' terminator (allow a window of ~50 bp at the annotated terminator)

**Strength (Chen et al. 2013, Nature Methods)** — optional, better as a weight than a binary.

https://doi.org/10.1038/nmeth.2515

- ΔG_U = RNA:DNA hybrid of the 8-nt U-tract vs template (Sugimoto RNA:DNA stacks)
- ΔG_L = loop-closure energy (this tracked termination strength; hairpin MFE ΔG_H did not)
- Threshold: e.g. Chen TS ≳ 10, or crude gate ≥5 U in 8 and hairpin MFE < −7 kcal/mol

Cambray et al. 2013 is a second calibrated strength set if you want a cross-check.

NRP Maker can generate non-repetitive hairpin+U-tract parts when Design needs to insert a terminator. That is part design, not internal-hit detection.

Design-mode kill: synonym-swap the U-tract (T→C/A/G) or break GC pairs in the stem. U-tract is the easier lever.

### 13.3 Rho-dependent

No Salis model. Use published motif detectors:

- RhoTermPredict (Di Salvo et al. 2019). https://doi.org/10.1186/s12859-019-2704-x
  78-nt Rut, C/G > 1, C every 11–13 nt, pick max C/G in a 128-nt neighborhood, then a pause (intrinsic-like hairpin or pyrimidine tract) within ~50–100 nt.
- Nadiras et al., NAR 2018 OPLS-DA. https://doi.org/10.1093/nar/gky705
  C>G skew, regularly spaced CC/UC, low structure (~85% classification).

Internal rho hits are common in high-CAI C-rich recoding. Balanced codon table is a side-effect mitigation. Design kill: drop C/G in the window or put structure on the Rut.

### 13.4 Design-mode use of a good terminator

If the user did not pass one: append a Chen-class intrinsic terminator from a non-repetitive toolbox (NRP Maker, L_max = 12 vs the rest of the operon). Score it with the intrinsic detector to confirm it is strong, and with Module I to confirm it doesn’t clone a 12-mer. That matches “non-repetitive … terminators may be inserted at the end.”

### 13.5 What not to do

- Do not put terminator strength into the mRNA half-life formula for E. coli (Cetnar 2021).
- Do not treat the intended 3' terminator as an internal hit.
- Do not skip the reverse strand.
- Do not wait for a Salis terminator energy model; it was not in the published Operon stack as a numbered calculator, only as a rule.

---

## 14. Module I — Repeats and genetic stability

Independent. Sources:

- Hossain et al., Nat. Biotechnol. 2020. https://doi.org/10.1038/s41587-020-0584-2
- NRP Calculator in https://github.com/hsalis/SalisLabCode
- RepeatFinder used inside the Synthesis Success Calculator

The Operon Calculator’s published objective is “fewer repetitive DNA sequences above 12 bp.”

### 14.1 What to report (Predict)

From the web tool: longest repeat, repeat-length histogram, chord diagram, full list with type/length/sequence.

| Type | Detection |
|---|---|
| Direct | identical k-mers, k ≥ 12, two or more loci, either strand |
| Reverse-complement (inverted) | k-mer and its RC |
| Tandem | consecutive copies, period ≥ 1, span ≥ 12 |
| Terminal | tandem of length ≥ 5 at either end of the fragment (synthesis + recombination) |

Algorithm: seed with 12-mers (canonical strand + RC as the key), extend matches, merge overlapping, keep unique pairs. O(L) expected with a dict of 4^12 keys. For L ~ 10 kb operons it is trivial.

L_max = 12 is the published cutoff. Homologous recombination in E. coli becomes efficient around 20–25 bp; 12 bp is a conservative design rule that also helps DNA synthesis (Module J).

### 14.2 Design: constrain the mutators

The Nonrepetitive Parts Calculator’s Maker mode generates parts that share no k-mer of length `homology` with each other or with a background (genome, plasmid backbone, already-accepted cistrons). Lighter version inside the GA:

- Maintain used_kmers of all 12-mers in the current operon plus optional background.
- A codon swap or RBS mutation is illegal if it introduces a 12-mer already in used_kmers (unless it only duplicates within the same local tandem you are already breaking).
- When designing multiple RBSs / promoters / terminators for a multi-operon system, run Maker-style generation once to get a toolbox of non-repetitive parts, then pick from the toolbox.

Finder mode of NRP (graph + vertex cover on the homology graph) is for selecting a subset of already-designed parts. Useful if you generate 100 RBS candidates and want the largest non-repetitive subset.

### 14.3 Transposon insertion sites and phage att sites

Independent motif scan. 2016 abstract rules (xii)/(xiii).

Ship a table, per host, of:

- IS element inverted-repeat cores (IS1, IS2, IS3, IS5, IS10, IS50, IS186, Tn3, etc. for E. coli K-12; different list for B. subtilis, P. putida).
- attB / attP cores for relevant phages (λ attB GCTTTTTTATACTAAC family, 186, P4, etc.).
- Chi sites if you care about RecBCD (GCTGGTGG in E. coli) — not in the original 15, easy expansion.

Scan both strands, exact or 1-mismatch. Design-mode kill: a single synonym swap in the motif. These motifs are rare enough that they almost never fight TIR.

---

## 15. Module J — Synthesis complexity and restriction sites

Independent. Sources:

- Halper, Hossain & Salis, ACS Synth. Biol. 2020. https://doi.org/10.1021/acssynbio.9b00460
- Synthesis class in NRP Maker; seq_eval.py in Synthesis Success Calculator
- Web tool: “reduce DNA synthesis complexity (homopolymers and extreme GC%)”

### 15.1 Hard rules (cheap, use as GA constraints not just scores)

| Rule | Default |
|---|---|
| Homopolymer C/G | no run of 9 |
| Homopolymer A/T | no run of 13 |
| Trinucleotide tandem | no 6× period-3 |
| Dinucleotide tandem | no 10× AT/AG/AC/… |
| Windowed GC, 20 bp | 15–90% |
| Windowed GC, 100 bp | 28–76% |
| Terminal 30 bp GC | 24–76% if L > 60 |
| Hairpins | no long / GC-rich / terminal hairpins |
| i-motif / G-quadruplex motifs | flag; G-quads are also a stability feature (Module D) |
| Local repeat density | 90% of a 70-bp window in some repeat is a synthesis fail; 60% of a 500-bp window; 69% of the whole fragment; a single repeat occupying 40% |

SSC trains a random forest on these features (F1 0.928 on 251 unseen fragments). You do not need the forest to start: the hard rules already define a feasible set. Embed the forest later as a smooth objective P(synthesis success) if you are designing >1 kb fragments.

SSC’s most important feature class was repeats. Module I and Module J are coupled in practice: killing 12-mer repeats does most of the synthesis work.

### 15.2 Restriction sites

User-supplied enzyme list. Exclude recognition sequences on both strands (unless the enzyme is not palindromic, then both orientations). Hard constraint: n_forbidden_RE > 0 is invalid. Optionally insert requested RE sites at part boundaries by making those nucleotides constant in the IUPAC RBS / intergenic constraint.

### 15.3 Design-mode use

Halper et al. showed you only need to recode ~5% of a hard CDS, targeted at the nucleotides the forest says matter, to push 73% of unsynthesizable E. coli CDSs into the easy bin. Inside the GA: if SSC flags a region, raise the mutation rate in that window.

ΔG_stacking = 0 in RBSForge means polyA/polyU spacers (common when you also want mRNA stability / RNase rules) will be mis-ranked on TIR. Leave stacking at 0; treat TIR as ordinal if spacer homopolymers dominate a Pareto set.

---

## 16. Module K — System-level RNAP and ribosomal load (optional, original rule xv)

Not a big control on the public web form, but in the 2016 list. Independent bookkeeping:

```
RNAP_load   = sum_over_promoters( Tx_rate * copy_number )
ribo_load   = sum_over_CDS( TIR * mRNA_level * L_codons / r_elong )
mRNA_level  = Tx_rate / (k_decay + μ)     # μ = growth rate
```

Use as a soft cap when the operon is one of many on a plasmid, or when you are overexpressing several CDSs with the Highly Translated table. It does not change sequence locally; it changes which Pareto designs you accept.

---

## 17. Predict mode: the whole pipeline in order

```python
def predict(operon, calc) -> PredictReport:
    A = assemble(operon)
    r_coupled = coupled_tirs(A, operon.host, calc)          # Module B
    r_mono    = [predict_tir(A.mrna, s, calc) for s in annotated_starts]
    ter       = [ter_of(cds) for cds in operon.cds_list]    # Module C
    htisc     = scan_htisc(A.mrna, annotated_starts, calc)  # Module F
    pauses    = scan_pauses(A.dna, operon.host, calc)       # Module G
    stability = mrna_stability(A, r_coupled, operon.host)   # Module D
    promoters = scan_promoters(A.dna, operon.host)          # Module E
    terms_i   = scan_intrinsic_terminators(A.dna)           # Module H
    terms_r   = scan_rho_terminators(A.dna)
    repeats   = find_repeats(A.dna, k=12)                   # Module I
    is_hits   = scan_motifs(A.dna, host.is_motifs + host.att_motifs)
    re_hits   = scan_motifs(A.dna, operon.forbidden_re)
    synth     = synthesis_features(A.dna)                   # Module J
    return PredictReport(...)
```

No search. The “physics-only” property Salis described: you can run this on an unannotated mRNA by taking every start codon as a candidate CDS and letting HTISC + coupling produce the ORF plot.

Host organism is an input even in Predict, because 16S anti-SD, codon weights, IS tables, and Promoter Calculator β all depend on it.

Evaluate mode takes mRNA. Promoter sequence is not in the TIR calculation; it only changes isoforms (Module E) and therefore which 5' UTR the RBS calculator sees. If a user changes the promoter and the RBS strengths move, either the TSS moved and the 5' UTR changed, or they were looking at coupled downstream TIRs while upstream TIR also changed. Both are real.

---

## 18. Design mode: search, not a closed-form inverse

There is no inverse of Module B. Design is sample sequences, assemble, run Predict, score, mutate.

### 18.1 Decision variables

| Variable | Domain | Notes |
|---|---|---|
| RBS_i nucleotides | IUPAC constraint, typically 15–30 nt | CDS1 is a 5' UTR; CDS_{i>1} may overlap CDS_{i−1} |
| Synonymous codon at each CDS position | amino-acid-preserving | frozen positions allowed |
| Intergenic extra nucleotides | optional short IUPAC spacer | determines d and ΔG_coupling |
| Promoter / terminator | usually fixed user parts, or picked from a non-repetitive toolbox | not recoded by the GA unless you opt in |

Protein sequence is invariant. DNA sequence is not.

### 18.2 Objectives (all independently toggleable, matching the web tool)

```
f1  = sum_i |log y.tir[i] − log target_tir[i]|     # TIR error, log-space
f2  = −sum_i y.ter[i]                              # or |ter − ter_target| for balanced
f3  = −y.stability                                 # higher mRNA half-life is better
f4  = n_internal_promoters
f5  = n_internal_terminators
f6  = n_htisc
f7  = max(0, longest_repeat − 12)
f8  = n_pause_sites
f9  = n_is_att_sites
f10 = n_rnase_sites                                # 5' UTR AU-ssRNA tracts
f11 = n_forbidden_re                               # hard: must be 0
f12 = synthesis_fail_score
```

The web tool returns several designs with “equally optimal characteristics,” which is Pareto language. Use NSGA-II (or any Pareto GA). Scalarizing with weighted sums will hide tradeoffs the user is supposed to see (TIR vs repeats vs HTISC).

Hard constraints (repair or reject, do not put on the Pareto front):

- n_forbidden_re == 0
- protein sequence unchanged
- IUPAC RBS constraints satisfied
- no homopolymer beyond Module J limits (or let them sit on f12)

### 18.3 Inner loop for one candidate

```
1. write RBS sequences and recoded CDSs into the template
2. assemble
3. optionally call design_rbs() as a proposal for RBS_i
   (mono-cistronic; do not trust its TIR for i>1)
4. run Predict, including coupling
5. return objective vector
```

Practical 2-stage scheme:

**Stage A — decoupled proposal.** For each CDS, run design_rbs(cds, target_tir, constraint) on the recoded CDS as if it were mono-cistronic. Choose junctions with d ≤ −25 (insulation) or a designed AUGA plus no overlapping hairpin if you want a little re-initiation. Recode CDSs with the selected codon table, applying Modules E–J as local repair.

**Stage B — coupled re-evaluation and GA polish.** Assemble, run full Predict with Module B. If CDS_{i>1} TIR drifted because of residual coupling or because recoding changed the footprint, mutate RBS_i / nearby codons with NSGA-II until the Pareto front is acceptable.

Stage A alone is already useful. Stage B is what the coin cost (50 coins vs 1 for an RBS prediction) is paying for.

If design_rbs is not built yet, Stage A is: plant AGGAGG at s_opt, randomize the rest of a 35-mer inside the IUPAC mask, go straight to Stage B.

### 18.4 Mutation operators (keep them local)

1. RBS nucleotide flip inside the IUPAC mask. Biased toward SD-like at 5–9 nt upstream of start for CDS1; less biased for overlapping RBS_i (the SD may sit on the upstream protein).
2. Codon swap to a synonym, biased by the active table, with a repair pass for RE / 12-mer / homopolymer.
3. N-terminal recode + RBS redesign as a paired move (they interact through the RBS calculator footprint).
4. Junction recode of the last 5–10 codons of CDS_i plus the first 5 of CDS_{i+1}, to hit a d and a ΔG_coupling target.
5. Targeted repair: if Predict flags a cryptic promoter at position p, only mutate the −10/−35 of that hit.

Crossover: uniform over cistrons (swap CDS_2 of parent A with CDS_2 of parent B), not over nucleotides. Nucleotide-level crossover wrecks SDs.

### 18.5 Why multiple designs are returned

The feasible set is large. Several sequences hit the same TIR vector with different repeat / HTISC / promoter leftovers. Returning a set lets the user pick by cloning convenience (RE sites, length) without re-running the search. Keep 5–20 non-dominated points, de-duplicated by 12-mer identity so you don’t return the same operon with one silent change.

---

## 19. Other RBSForge holes that bias Operon scores

These are smaller than design_rbs / footprint but belong in any faithful replica.

**v2.0 standby not wired.** StandbyParams (c1=0.038, c2=−1.629, c3=17.359, c_slide=0.20, h_cap=15, a0=15) sit on HostPack unused. Tian 2015’s ΔG_standby is the 2013 module-geometry term (Espah Borujeni, Channarasappa & Salis, NAR 2014, https://doi.org/10.1093/nar/gkt1139). Structured intergenic RBSs (the coupling case) are exactly where v1.0 vs v2.0 standby diverges. Highest-value RBSForge upgrade after footprint. Do not ship both v1.0 4-nt and v2.0 as additive.

Formula scaffolded on HostPack:

```
A_s = A0 + P + D − min(H, H_cap)
ΔG_distortion = C1*A_s^2 + C2*A_s + C3
```

plus c_slide per nt of downstream hairpin height.

**ΔG_stacking = 0.0.** Unpublished v2.1 homopolymer-spacer term (Reis & Salis 2020, https://doi.org/10.1021/acssynbio.0c00394). Honestly zeroed. Do not invent a coefficient.

**No ribosome drafting (2016).** Espah Borujeni & Salis, JACS 2016. https://doi.org/10.1021/jacs.6b00216. Spec already said not to put this in coupling. Fast-folding 5' UTRs may be over-penalized.

**NN parameters are Xia 1998 / Turner, not Andronescu 2007 (v2.1).** Relative ranking inside one engine is fine. Do not compare RBSForge ΔG to a Salis web-tool ΔG term-by-term.

**Leaderless starts:** tir=None. Coupling: de novo = 0, re-init only. Design: infeasible.

**Window coords:** duplex is window-local; GenBank needs window_start + offset.

**HostPack ≠ OperonHost.** Wrap; do not fork; do not put codon usage or IS motifs into RBSForge.

**Non-E. coli HostPacks are Layer-1 priors.** aSD swap is real (b_subtilis 15-nt tail is the important one). β, spacing curvature, standby, start table, C, k_P, k_reinitiation(d) are E. coli 37 C. A B. subtilis operon run through this stack produces rankings, not au. Gram-positive RNase J/Y vs E. coli RNase E is a Module D problem.

**asd_modifications (Thermus Ψ)** are flagged and not applied. Ignore the Thermus pack for Operon work.

**MAX_PRE_WINDOW = 200** is enough for a junction (tens of nt upstream of start_i). Fine.

---

## 20. Implementation order

If RBSForge already works, implement in this order. Each step is a complete, testable program. Do not start with the GA.

| Step | What | You can ship | Test |
|---|---|---|---|
| 1 | Assembler + d | Concatenate and annotate GenBank | Round-trip a known bi-cistronic (AUGA → d = −4) |
| 2 | TIR scale freeze (v1_style_rate) + predict_tir adapter | 1-cistron TIR | Known strong vs weak RBS ranking |
| 3 | HTISC | ORF plot | Any mRNA with an internal AUG+SD |
| 4 | predict_one + forced_unpaired + footprint 13 | Bound-state unfolding | Hairpin in +4..+13 hurts less than hairpin in SD |
| 5 | Coupling, two predict_one calls | Coupled TIR for n CDSs | Tian 2015: re-init constructs; hairpin constructs. Re-fit C, k_P. Vienna required. |
| 6 | design_rbs (replace-only, IUPAC) | Inverse RBS for CDS1 | Hit a target TIR within 0.25 kcal/mol on a reporter CDS |
| 7 | Codon recode + TER | Recoded CDS, TER in nt/s | Protein unchanged; CAI up under Highly Translated |
| 8 | 12-mer repeat finder | Longest repeat, chord list | Plant a 15-mer twice |
| 9 | RE + homopolymer + GC | Hard-constraint filter | Known synthesis-fail sequences from SSC SI |
| 10 | Promoter Calculator scan | Internal promoter list | J23100 as a strong TSS; TATAAT…17 bp…TTGACA flags |
| 11 | Intrinsic terminator scan | Internal terminator list | Known rrnB T1 / Chen strong terminator |
| 12 | Pause / internal SD | Pause list | Plant AGGAGG in-frame |
| 13 | 2021 stability | Half-life / decay score | Weak RBS vs strong RBS on the same CDS: stability drops with TIR |
| 14 | NSGA-II over RBS+codons | Design mode | Recover a target TIR vector on 2–3 fluorescent reporters |
| 15 | 2024 GBDT, rho, IS tables, load, v2.0 standby | Current web-tool parity | Cetnar 2024 test split if you retrain |

Do not wait for ΔG_stacking or v2.0 standby to start. Ship insulated (d ≤ −25) multi-cistronic Predict with HTISC + scanners. Add constrained-fold coupling next. Inverse design after that. Stacking never, unless Salis publishes the coefficient.

---

## 21. Worked information flow (2-cistron Design)

User inputs: host = E. coli MG1655, promoter = J23100, terminator = BBa_B1006, CDS1 = enzyme A (AA sequence), CDS2 = enzyme B (AA sequence), target TIR = (10^4, 10^3) on the v1 au scale, table = Balanced, exclude EcoRI/XbaI, rules = all on, RBS constraints = 20 N's for both.

1. Recode A and B to DNA with the Balanced table mixed at log(1e4)/log(tir_max) and log(1e3)/log(tir_max). Until footprint is on, freeze the first ~5 codons or accept that they will be jointly redesigned with the RBS.
2. Choose junction d = −25 (insulation — the targets differ by 10× and you probably want them independent). Overlap the last 25 nt of recoded A with the start of B; this forces a joint recode of that window so both proteins stay correct. If a joint recode is impossible, drop to d = −4 (AUGA) and accept a known re-initiation leak of ~0.022 × r_1 into r_2.
3. design_rbs(A_nt, 1e4) → RBS1. design_rbs(B_nt, 1e3, constant_upstream=A_tail) → RBS2 (the overlapping tail is the constant upstream).
4. Assemble, Predict:
   - Module B may show r_2 ≠ 1e3 because the overlap changed ΔG_mRNA. Adjust RBS2 / last codons of A.
   - Modules E–J flag leftover hits. Local synonym swaps.
   - Module D: if half-life is poor, raise r_1 or add a 5' hairpin that doesn't steal from RBS1.
5. If flags remain, NSGA-II with the mutation operators in §18.4 for a few hundred generations, population ~50, return the Pareto set.
6. Emit GenBank with features: promoter, TSS, RBS1, CDS1, overlap, RBS2, CDS2, terminator, plus every cryptic hit as misc_feature. SD spans use window_start + duplex offsets.

If instead the user wants B's expression to track A (a coupled ratio), skip insulation. Set d = −4, design an overlapping hairpin with ΔG_coupling from Tian's Figure 6 for the desired ratio, and let RBS2 be weak in the folded state. Then r_2(r_1) is the sigmoid of Module B, and Design's TIR objective is on the pair, not on r_2 alone.

Vienna must be on PATH for step 4.

---

## 22. How to expand each module without touching the others

Every expansion is a drop-in replacement of one box in §2.

| Want | Replace / extend | Do not touch |
|---|---|---|
| New host | HostPack aSD + T, plus OperonHost codon tables, IS/att, Promoter Calculator β (or disable E for non-σ70), RNase E vs Y | Coupling equations (C and k_reinitiation should be re-fit; they were measured in E. coli) |
| Re-fit coupling | New reporter library → new C, k_P, optionally k_reinitiation(d) | RBSForge internals |
| T7 / σS / σ32 promoters | Swap Module E | Translation |
| Better decay | Swap Module D's scorer for the 2024 GBDT, or add RNase III / RppH tables | Coupling, except that D consumes TIRs |
| Rho-only hosts | Weight Module H toward Rut scanning; intrinsic still exists | — |
| Eukaryotic IRES / Kozak | This architecture does not apply. Coupling via 70S re-initiation is bacterial. | Everything |
| Riboswitches / toeholds in a UTR | Freeze those nts in the RBS IUPAC mask / constant_upstream. Re-run Predict. | GA will respect the mask |
| Multi-operon plasmids | Run Design per operon with a shared used_kmers background (plasmid + already designed operons). NRP Maker's background argument. | — |
| Libraries rather than single sequences | Replace design_rbs with an RBS Library Calculator (Farasat 2014) per cistron, still re-score coupling on sampled members. Keep junctions non-degenerate. | — |
| Load / burden | Turn on Module K as a cap | Local sequence scores |
| ML TIR instead of RBSForge | As long as it returns a scalar TIR and you can evaluate two sequences (folded / unfolded) you can keep Module B's mixture. You lose the ΔG_coupling split-by-helix unless the ML model explains structure. | — |
| v2.0 standby | Wire StandbyParams; drop v1.0 4-nt | Coupling mixture |
| Inverse RBS | §5 design_rbs | Do not put coupling inside it |

The clean expansion surface is the OperonHost record plus Predict's list of scanners. Add a scanner, add an objective, add a mutator that knows how to kill that scanner's hits. That is the whole design-rule plugin.

---

## 23. Parameters to copy, and parameters to re-fit

**Copy (E. coli, do not invent):**

- Coupling ratios: k_reinitiation(−4) = 0.022, k_reinitiation([0,25]) = 0.0072, reverse-overlap factor 11.6 at 25 nt.
- Ribosome: occupancy footprint 30 nt, r_elong ≈ 60 nt/s, CDS unfold footprint 13 nt.
- Promoter Calculator: K = 42, β_vivo = 1.636, matrices from free_energy_coeffs.npy.
- Repeat cutoff: 12 bp.
- Homopolymer / GC windows: §15.1.
- RBSForge β = 0.45, K = 2500, start-codon table, spacing push/pull, s_opt = 5, aSD = ACCUCCUUA.
- design_rbs: tol = 0.25 kcal/mol, RT_init = 0.6, accept-ratio band [0.01, 0.20], MaxIter = 10000, replace/insert/delete weights 0.80/0.10/0.10, rbs_max_len = 35.

**Re-fit if the host or the TIR au scale changes:**

- C and k_P. They mix physical occupancy with the au scale. A different RBS calculator (different K, different β, LinearFold vs Vienna) changes C. Recalibrate on a small bi-cistronic reporter panel.
- Stability coefficients {a, b} in §9.2, or the GBDT: E. coli RNase E / RppH specific.
- Codon tables, always.
- β itself for non-E. coli HostPacks (Layer 2).

**Do not copy from RBSForge into coupling:** standby-site sliding, RNA folding kinetics (ribosome drafting), and NEQ warnings. Those already affected r_folded / r_unfolded. Coupling only mixes those two numbers.

**Do not invent:** ΔG_stacking coefficient, Thermus Ψ stacking bonus, Gram-positive v2.1 spacing/standby tables (unpublished).

---

## 24. Known failure modes

1. Evaluating coupling on the promoter DNA. Evaluate mode takes mRNA. Promoter sequence is not in the TIR calculation; it only changes isoforms.
2. Trusting design_rbs for cistron 2+. Mono-cistronic design ignores ΔG_coupling and re-initiation. Always re-predict on the assembled mRNA.
3. Recoding the N-terminus without re-predicting TIR. The footprint is ~13 nt / ~5 aa, and the RBS calculator window includes ~35 nt of CDS. Silent recoding is not silent there.
4. Insulating with a long unstructured spacer. Kills coupling, creates an RNase E site, may create a promoter. Prefer 25 nt coding overlap.
5. Highly Translated table on a 10-gene operon. Ribosomal load and synthesis complexity explode; cryptic promoters appear because high-CAI sequence is hexamer-biased. Use Balanced.
6. Treating terminator strength as an mRNA-stability lever in E. coli. The 2021 paper says it is not, for upstream mRNA.
7. σ70 Promoter Calculator on a T7 cassette. Wrong polymerase, confident cryptic-promoter calls that do not exist (and missed T7 cryptics that do).
8. Forgetting the reverse strand. Promoters, terminators, IS motifs.
9. Pareto vs weighted sum. If you scalarize, you will never surface the design that missed TIR by 1.4-fold but has zero internal promoters. That is often the one to clone.
10. Mixing proportional_rate and v1_style_rate. Occupancy and HTISC thresholds become junk. One au, one τ, one C.
11. Running coupling on the builtin folder. Intergenic hairpins will be under-stabilized and the sigmoid will flatten. Require Vienna.
12. Copying C = 0.81 onto RBSForge au without a reporter re-fit.
13. Double-counting standby 4-nt inside occupied when footprint unfolding is added.
14. Returning leaderless starts as TIR = 0 instead of infeasible / de novo = 0.
15. Comparing RBSForge ΔG term-by-term to the Salis web tool (Xia vs Andronescu, v1.0 standby vs v2.0, stacking = 0, no drafting, no footprint until you add it).

---

## 25. Minimal reference interfaces

```python
# rbsforge (exists)
class RBSCalculator:
    def predict(self, mrna: str) -> PredictResult: ...
    # add:
    def predict_one(self, mrna, start, extra_unpaired=None) -> StartCodonResult: ...

def design_rbs(cds_nt, target_tir, hostpack, constraint_iupac=None,
               constant_upstream=None, **kw) -> DesignResult: ...

# coupling.py
def k_reinitiation(d: int) -> float: ...
def coupled_tirs(assembled, host, calc) -> list[float]: ...

# elongation.py
def recode_cds(aa: str, table, rng) -> str: ...
def ter_nt_s(cds_nt: str, host) -> float: ...

# stability.py
def mrna_stability(assembled, tirs: list[float], host) -> StabilityResult: ...

# scans.py
def scan_promoters(dna, host) -> list[PromoterHit]: ...
def scan_intrinsic_terminators(dna) -> list[TermHit]: ...
def scan_rho_terminators(dna) -> list[TermHit]: ...
def scan_htisc(mrna, annotated, calc, tau) -> list[StartHit]: ...
def scan_pauses(cds_nt, host, calc) -> list[PauseHit]: ...
def find_repeats(dna, k=12) -> RepeatReport: ...
def scan_motifs(dna, motifs) -> list[MotifHit]: ...
def synthesis_score(dna) -> SynthReport: ...

# assemble.py
def assemble(operon) -> Assembled: ...

# predict.py
def predict(operon, calc) -> PredictReport: ...

# design.py
def design(operon, calc, objectives, constraints) -> list[Design]: ...
```

Each scans.py function is a self-contained file in a real codebase, with its own tests and its own paper. That is the unit of replication. Coupling is the only function that must call the RBS calculator more than once per start codon. Everything else either calls it once per start (HTISC) or not at all.

---

## Appendix A. design_rbs algorithm (compact)

```
dG_target = 2.222 * (7.824 - ln(TIR_target))
O(seq)    = |predict(upstream + RBS + cds).dG_total - dG_target|

init:
  RBS = random 35-mer
  plant aSD-complementary core at spacing s_opt, length biased by how negative dG_target is
  strip start codons
  if O > 10 kcal/mol: redraw

loop up to MaxIter:
  move ~ {replace: 0.8, insert: 0.1, delete: 0.1}
  if IUPAC fixed-length: replace only, new base in mask[pos]
  strip start codons; clamp length to rbs_max_len
  if leaderless or IUPAC broken: reject
  O_new = O(RBS_new)
  if O_new < O or rand() < exp((O - O_new)/RT): accept
  every 50 moves:
    if accept_ratio > 0.20: RT /= 2
    if accept_ratio < 0.01: RT *= 2
  if O <= 0.25: return hit_tol=True

return best-so-far, hit_tol=False
```

Constants from v1.0: MaxIter=10000, tol=0.25, RT_init=0.6, annealing_min_moves=50, annealing_accept_ratios=[0.01, 0.20], max_init_energy=10.0, dG_range_high=25.0, dG_range_low=−18.0, Max_RBS_Length=35.

Skip v1.0 kinetic_score and three_state_indicator until those quantities exist.

IUPAC alphabet: N=ACGT, R=AG, Y=CT, S=GC, W=AT, K=GT, M=AC, B=CGT, D=AGT, H=ACT, V=ACG, and the four literals.

---

## Appendix B. Constrained-fold / footprint API

```
folder.fold(seq, forced_unpaired: FrozenSet[int]) -> FoldResult(delta_g, pairs, backend)

occupied_default(start_in_window, duplex, footprint_cds=13):
    return set(range(start_in_window, start_in_window + 3 + footprint_cds))
    # optionally union SD span and 4-nt standby; do not double-count standby energy

unfolding = fold(window).delta_g - fold(window, occupied ∪ extra).delta_g
```

Vienna: RNAfold --constraint, constraint string with `x` at occupied positions, `.` elsewhere.

Coupling extra_unpaired: every index < start_i that is part of CDS_{i−1} including the last nucleotide of the stop, clipped to the window.

Tests:

- No structure in window → unfolding ≈ 0, TIR unchanged vs today.
- Hairpin wholly in footprint → TIR increases relative to unconstrained-MFE penalty.
- Hairpin wholly in SD → TIR still low.
- Forcing upstream-CDS nts unpaired on a known coupling hairpin → r_unfolded >> r_folded.

---

## Appendix C. Coupling equations

```
r_1 = k * exp(−β ΔG_total^(1))                 # RBSForge, first CDS

k_reinit(d) as piecewise in §7.2               # transfers across au scales

r_reinit_i = k_P * k_reinit(d) * r_{i-1}       # k_P re-fit

f = min(1, C * r_{i-1,phys})                   # C re-fit
r_denovo_i = (1-f) * r_folded + f * r_unfolded
r_i = r_reinit_i + r_denovo_i
```

Tian published: C = 0.81 ± 0.17, k_P = 10 (best-fit 14, CI [8, 33]). These are not portable onto RBSForge without a reporter library.

Regimes: ΔG_coupling = 0 → linear leak; overlapping hairpin → sigmoid in log r_{i−1}; d ≤ −25 → insulation.

---

## Appendix D. Promoter Calculator scan (σ70)

Window ~78–89 bp. Scan both strands. Keep min dG_total per TSS. Tx_rate = 42 * exp(−1.636217 * dG_total) in vivo.

Spacer penalty 0.1463 s^2 − 4.9113 s + 41.119, s in 15..20.

346-parameter matrices: free_energy_coeffs.npy in SalisLabCode/Promoter_Calculator.

Internal if TSS ≠ intended and Tx_rate > τ_tx.

---

## Appendix E. Terminator detectors

Intrinsic presence: stem 4–15, loop 3–8, U-tract ~15 nt, both strands. Strength: Chen ΔG_U (8-nt RNA:DNA) and ΔG_L (loop closure), not ΔG_hairpin. Internal if not the annotated 3' terminator.

Rho: 78-nt Rut, C/G > 1, C every 11–13 nt, max C/G in 128-nt neighborhood, pause within 50–100 nt. Or Nadiras C>G + CC/UC periodicity + low structure.

Do not feed terminator strength into E. coli mRNA half-life.

If user omitted a terminator: append a non-repetitive Chen-class intrinsic part (NRP Maker, L_max=12 vs the operon).

---

## Appendix F. HostPack vs OperonHost

HostPack (RBSForge, organism-agnostic TIR machinery):

- asd, t_growth_c, gram_stain
- s_opt, spacing_push, spacing_pull
- standby StandbyParams (v2.0 numbers unused)
- start_codon_dg, beta
- footprint_cds=13 (unused until §6), cutoff_post=35, cutoff_pre=None, rbs_max_len=35
- nn_param_set, asd_modifications (flagged, unused)

OperonHost adds:

- codon_table, codon_weights_high, codon_weights_balanced
- r_elong_nominal_nt_s ≈ 60, ribosome_footprint_nt = 30
- C_unfold, k_P (re-fit)
- is_motifs, att_motifs, forbidden_re_sites
- rnase_e_like
- Promoter Calculator β / K if you cache them per host

Do not put OperonHost fields into RBSForge.

---

## Appendix G. Gap checklist

| Gap | Blocks | Fix |
|---|---|---|
| No design_rbs | Spec Stage A | Build Appendix A around RBSForge.predict |
| Constrained fold not public | Hairpin-mediated coupling and footprint | predict_one + forced_unpaired |
| footprint_cds unused | N-terminal recode vs TIR; bound-state ΔG_mRNA | Two-fold, Appendix B |
| C, k_P on a different au | Occupancy f_unfold | Re-fit on v1_style_rate; k_reinitiation(d) keeps |
| Builtin folder too simple | Coupling ΔG | Require Vienna when coupling or footprint is on |
| Standby v1.0 only | Structured intergenic RBS | Wire existing StandbyParams; do not add to v1.0 |
| ΔG_stacking = 0 | Homopolymer spacers | Leave at 0 |
| No drafting | Fast-folding UTRs | Out of scope |
| No helices on StartCodonResult | Optional helix-split of ΔG_mRNA | Unnecessary if constrained fold exists |
| duplex coords window-local | SD annotation in GenBank | Add window_start (already there) |
| Leaderless tir=None | SD-less overlapping CDS | de novo = 0, re-init only; design infeasible |
| HostPack ≠ OperonHost | Codons, IS, r_elong | Wrapper |
| Non-E. coli β/C uncalibrated | Absolute TIR / coupling in other hosts | Rankings only |
| No stability / promoter / terminator in RBSForge | Rest of Operon | Still yours; terminator is a scan, not a Salis calculator |
| HTISC threshold in mixed units | False truncated-protein calls | One au, one τ |

---

## Appendix H. Full source list (URLs as text)

Operon Calculator product and protocol:

- https://docs.denovodna.com/docs/operon-calculator
- https://docs.denovodna.com/docs/operon-calculator.md
- https://github.com/hsalis/salis-lab-protocol-book/blob/master/design/operon-calculator.md
- https://raw.githubusercontent.com/hsalis/salis-lab-protocol-book/master/design/operon-calculator.md
- https://www.denovodna.com/software/predict_operon_calculator
- https://www.denovodna.com/software/design_operon_calculator
- https://salislab.net/software/design_operon_calculator
- https://www.denovodna.com/software
- https://salislab.net/software

RBS Calculator product:

- https://docs.denovodna.com/docs/rbs-calculator
- https://docs.denovodna.com/docs/rbs-calculator/rbs-calculator-free-energy-model.md
- https://www.denovodna.com/software/predict_rbs_calculator
- https://www.denovodna.com/software/design_rbs_calculator

RBSForge:

- https://github.com/Neladden/RBSForge-branch
- https://github.com/Neladden/RBSForge-branch/pull/1

Salis public source:

- https://github.com/hsalis/Ribosome-Binding-Site-Calculator-v1.0
- https://github.com/hsalis/Ribosome-Binding-Site-Calculator-v1.0/blob/master/RBS_MC_Design.py
- https://github.com/hsalis/Ribosome-Binding-Site-Calculator-v1.0/blob/master/RBS_Calculator.py
- https://github.com/hsalis/SalisLabCode

Papers, RBS model:

- Salis, Mirsky & Voigt, Nat. Biotechnol. 2009. https://doi.org/10.1038/nbt.1568
- Salis, Methods Enzymol. 2011. https://doi.org/10.1016/B978-0-12-385120-8.00002-4
- Espah Borujeni, Channarasappa & Salis, NAR 2014 (v2.0 standby). https://doi.org/10.1093/nar/gkt1139
- Farasat et al., Mol. Syst. Biol. 2014 (RBS Library Calculator). https://doi.org/10.15252/msb.20134955
- Espah Borujeni & Salis, JACS 2016 (drafting). https://doi.org/10.1021/jacs.6b00216
- Espah Borujeni et al., NAR 2017 (footprint). https://doi.org/10.1093/nar/gkx262
- Reis & Salis, ACS Synth. Biol. 2020 (v2.1 / model test). https://doi.org/10.1021/acssynbio.0c00394

Papers, Operon / coupling / stability / promoters / parts / synthesis:

- Tian & Salis, NAR 2015 (coupling). https://doi.org/10.1093/nar/gkv635
- Tian & Salis PMC. https://pmc.ncbi.nlm.nih.gov/articles/PMC4538824/
- Cetnar, SEED 2016 abstract. https://proceedings.aiche.org/sbe/conferences/synthetic-biology-engineering-evolution-design-seed/2016/proceeding/paper/operon-calculator-automated-design-synthetic-operon-sequences-using-15-models-and-design-rules
- Cetnar & Salis, ACS Synth. Biol. 2021 (mRNA stability). https://doi.org/10.1021/acssynbio.0c00471
- Cetnar & Salis 2021 bioRxiv. https://www.biorxiv.org/content/10.1101/2020.07.22.216051v1
- Cetnar et al., Nat. Commun. 2024 (isoform-aware decay). https://doi.org/10.1038/s41467-024-54059-7
- LaFleur, Hossain & Salis, Nat. Commun. 2022 (Promoter Calculator). https://doi.org/10.1038/s41467-022-32829-5
- Hossain et al., Nat. Biotechnol. 2020 (NRP). https://doi.org/10.1038/s41587-020-0584-2
- Halper, Hossain & Salis, ACS Synth. Biol. 2020 (Synthesis Success). https://doi.org/10.1021/acssynbio.9b00460

Papers, terminators:

- Chen et al., Nat. Methods 2013. https://doi.org/10.1038/nmeth.2515
- Kingsford, Ayanbule & Salzberg, Genome Biology 2007 (TransTermHP). https://doi.org/10.1186/gb-2007-8-2-r22
- Di Salvo et al., BMC Bioinformatics 2019 (RhoTermPredict). https://doi.org/10.1186/s12859-019-2704-x
- Nadiras et al., NAR 2018 (rho). https://doi.org/10.1093/nar/gky705

NN parameters / related:

- Xia et al., Biochemistry 1998 (RNA NN stacks). https://doi.org/10.1021/bi9809425

Other RBS reimplementation (not this stack, cited only as contrast):

- https://github.com/allyourbasepair/rbscalculator

---

## Bottom line

The Operon Calculator is those papers composed. Replicate them as independent functions, feed them an assembled operon, and put a Pareto search over RBS nucleotides and synonymous codons.

RBSForge is a legitimate predict_tir for E. coli ~37 C, transparent enough to build coupling on, and missing exactly the two Operon-facing features a naive spec assumes: inverse design, and constrained-unfold / footprint. Inverse design you copy from public v1.0 RBS_MC_Design around RBSForge.predict. Constrained unfold you add as one predict_one hook, then re-fit C. Vienna is required for coupling. ΔG_stacking stays zero. Terminators are scanners (Chen / TransTermHP / RhoTermPredict), not a missing Salis calculator. Everything else in the Operon Calculator was never RBSForge's job.
