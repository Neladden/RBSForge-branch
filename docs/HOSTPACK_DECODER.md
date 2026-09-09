# HostPack Decoder — Design Document v0.1

**Status: draft, for discussion.** This is a first design pass for a network
that consumes `genome_encoder`'s embeddings (Architecture Document v1.2,
hereafter "the encoder doc") and predicts RBSForge `HostPack` parameters —
the organism-specific constants in `rbsforge/hostpack.py` (anti-SD tail,
optimal spacing, β, standby coefficients, start-codon energies, growth
temperature, Gram stain) — de novo from genome sequence, rather than by hand
curation. It follows the encoder doc's own conventions (numbered sections,
falsifiable hypotheses, staged decisions register, open questions) because
that structure earned its keep on a harder problem than this one, and
because a shared vocabulary makes it easy to reason about where the two
documents' decisions touch each other — several do, and are called out
explicitly in §10.

**Relationship to the two existing systems.** `genome_encoder` is upstream
and independent: it does not know RBSForge exists, and nothing here should
require reopening its locked decisions. RBSForge is downstream and already
built: its `HostPack` dataclass and calibration-tier badge system
(`calibrated` / `prior` / `speculative`) are the actual consumer of this
network's output, not a hypothetical future integration — §9 makes this
concrete. This document sits between them.

---

## 1. Goals and Non-Goals

**Goal.** Given a `genome_encoder` embedding (and, where available, its
auxiliary outputs — MoE routing histogram, the learned codon-table
correction Δ, transcription-unit-scale operon/promoter/terminator
predictions) for a genome the RBS Calculator model has no hand-curated
profile for, produce a `HostPack` — every field RBSForge's scoring engine
needs — with each field carrying an honest confidence tier, not a single
opaque number. The network should get a genuinely uncharacterized organism
closer to right than RBSForge's current fallback (borrow the nearest
phylogenetic relative's numbers by hand), and should know when it can't.

**Non-goals.**

- **Genome-wide, gene-level TIR prediction.** A `HostPack` is a compact
  parameter set for RBSForge's existing thermodynamic model to consume, not
  a replacement for that model. Nothing here tries to predict individual
  gene expression directly from sequence.
- **Producing a `calibrated`-tier output from sequence alone, ever.**
  β, spacing curvature, and standby coefficients are empirical fits against
  real translation-rate measurements (see §2's Tier 3). No amount of
  embedding quality changes what that tier means; this network's job for
  those fields is a well-calibrated `prior`, not a substitute for
  wet-lab calibration. Stated as a hard non-goal because it's the mistake
  most worth naming explicitly before anyone builds against this system.
- **Replicon-level or gene-level HostPacks in v0.1.** Output is one
  HostPack per genome. §10 discusses why this is a real simplification
  (a plasmid can carry its own thermotolerance or regulatory machinery)
  and why it's deferred rather than solved now.
  A promoter-consensus/sigma-factor-strength decoder — genome_encoder's own
  transcription-unit-scale auxiliary heads (encoder doc §3.2/§4.4) are
  already targeting exactly this territory. If/when RBSForge's descoped
  promoter-strength module (see `docs/MODEL.md`, "What this package
  deliberately omits") is revisited, it should consume those heads' output
  rather than duplicate the effort with a separate motif detector.
- **Training on or requiring metagenome-assembled genomes (MAGs).**
  Inherits the encoder doc's own chimerism/contamination caveat (§1 there)
  directly and for the same reason, but the failure mode is worse here: a
  chimeric contig doesn't just degrade an embedding, it can confidently
  blend two organisms' physiology into one internally-consistent-looking
  but wrong HostPack. Out of scope until the encoder doc's own §9 entry 16
  trigger fires.
- **Joint or end-to-end training with `genome_encoder`.** Discussed and
  rejected in §7 — this network is a strictly downstream consumer of a
  frozen embedding, never a gradient source back into the encoder.

---

## 2. Output Schema: Three Tiers by Epistemic Status

This is the central design decision, and everything else in this document
follows from it. `HostPack`'s fields are not one homogeneous prediction
target — they differ in *kind*, and treating them uniformly (one network,
one loss, one confidence story) would either overclaim on the fields that
need real data or waste modeling effort on fields that don't need learning
at all.

### Tier 1 — Directly sequence-determined

Fields where the honest task is extraction (possibly with a learned
refinement step for a genuinely ambiguous boundary), not regression against
sparse ground truth.

| Field | What determines it | Learned component |
|---|---|---|
| `asd` (anti-SD tail) | Locate the 16S rRNA gene; read its 3′ end. | A boundary-refinement head — see below. |
| `gram_stain` | Cell-envelope gene content (peptidoglycan synthesis vs. outer-membrane/LPS biogenesis operons). | Classifier, abundant labels. |
| Standby-platform architecture (proxy for which `standby` regime applies) | Presence/copy number of ribosomal protein S1 and related platform proteins. | Classifier or direct HMM search, no real learning needed. |

**Why `asd` needs a learned step at all, given it's "just" extraction:**
RBSForge's own docs already flag this — GenBank 16S annotations are
frequently truncated before the anti-SD region (Nakov et al. 2018
re-annotated >12,000 3′ ends for exactly this reason). The task is not
"find the 16S gene" (solved, cheap, deterministic) but "given the region
around the annotated 3′ end, predict the true mature terminus," which is a
short-range sequence-to-boundary regression problem with real, if scarce,
training signal (§6). This is the one Tier-1 field that benefits from a
neural component; the other two are close to fully deterministic and are
included here mainly so the network learns to reconstruct them from the
*embedding alone* — useful for the degraded/incomplete-fragment case the
encoder doc's own Goal statement (§1 there) is designed to tolerate, where
a direct HMM search on a partial assembly might miss the relevant gene
entirely but the whole-genome embedding may still carry the signal.

### Tier 2 — Statistically inferable, self-supervisable

Fields with real comparative-genomic signal, refinable without wet-lab
calibration, but where the *scale* still ultimately traces back to a small
number of measured organisms.

| Field | Signal | Learned component |
|---|---|---|
| `t_growth_c` | Proteome amino-acid composition (IVYWREL-type indices), rRNA stem GC%, codon/tRNA GC3, purine-loading asymmetry. | Regression, real precedent (§6). |
| `s_opt` (optimal spacing) | The organism's own native SD-strength-vs-spacing distribution, computed classically. | A *correction* on top of a classical estimate — see "prior-anchored regression" below. |

**Prior-anchored regression for `s_opt`.** This is worth naming as a
pattern because it recurs and because it deliberately mirrors something the
encoder doc already does. §3.1 of the encoder doc initializes the
DNA-tower's ORF hypothesis from a heuristic caller as a *soft prior*, not a
hard constraint, and rewards the loop for diverging from it when the
protein-tower's own signal supports divergence. The `s_opt` head does the
same thing with a *classical* prior instead of a heuristic-caller one:

1. **Native-RBS survey (classical, no learning).** Once `genome_encoder`
   has called high-confidence ORFs and their upstream regions for a genome
   (its own ORF-discovery loop, encoder doc §4.3), run RBSForge's existing
   `thermo/duplex.py` search against each one, using the genome's own
   (Tier-1-predicted) anti-SD tail. This produces an empirical distribution
   of aligned SD:aSD spacing across the whole genome, weighted by duplex
   strength — exactly the "pull native 5′ UTRs genome-wide; compute SD:aSD
   ΔG at growth temperature — this is your prior on what 'strong' means in
   that species" step already described, for a human doing it by hand, in
   `docs/MODEL.md`'s Layer-0/Layer-1 guidance. The survey's mode is a
   genuinely strong `s_opt` estimate requiring zero wet-lab data — SD-usage
   statistics genome-wide are self-evidently a large-N signal within a
   single genome.
2. **Learned correction.** `s_opt_predicted = survey_mode +
   correction(E_genome, survey_summary_stats)`. The network only has to
   learn the *residual* — a much lower-variance target than the raw value
   — and degrades gracefully to the classical estimate alone when the
   correction network has no useful signal (e.g. at initialization, or for
   an embedding region far from anything seen in training).
3. Correction is needed because the raw survey can be systematically
   biased: under-called ORFs skew it, and some lineages (Deinococcus-
   Thermus, some Bacteroidetes) are leaderless-heavy enough that "the
   genome-wide mode" isn't a meaningful concept at all — which the
   correction network should learn to recognize (ideally by predicting a
   wide/degenerate output, not a confident wrong one) from the embedding.

**Operon-awareness the survey should have and currently doesn't, in
RBSForge itself.** Internal genes within a polycistronic operon are
subject to translational coupling (re-initiation after the upstream gene,
different SD-strength/spacing statistics than a genuine transcript-leading
gene) — RBSForge's own docs already flag this by name ("need the Operon
Calculator, not this model alone," `docs/MODEL.md`, Known Limitations).
`genome_encoder`'s transcription-unit-scale operon-boundary output (encoder
doc §3.2/§4.4, contingent on its own H6 validating) lets the native-RBS
survey stratify by first-gene-in-operon vs. internal-gene, rather than
averaging them into one number that describes neither well. This is not
required for a first version of the survey, but it's a concrete, cheap
win available *because* the two systems exist together, and worth flagging
as a reason to prefer waiting for validated H6 output over shipping the
survey blind to operon structure.

### Tier 3 — Requires real calibration data

| Field | Why it can't be shortcut |
|---|---|
| `beta` | An empirical ΔG→rate conversion fit against real reporter/expression data — not 1/RT, not derivable from sequence (`docs/MODEL.md` makes this point at length for the existing hand-curated HostPacks; it applies with equal force to a learned one). |
| `spacing_push` / `spacing_pull` (7 coefficients) | Curve *shape* assumed fixed by convention; organism-specific scale needs spacer-length reporter libraries. |
| `standby` `{C1,C2,C3,c_slide,H_cap,A0}` | Needs structured hairpin-occlusion reporter libraries. |
| `start_codon_dG` `{AUG,GUG,UUG,CUG}` | Apparent fMet-tRNA pairing energies, folding in IF3 discrimination effects — fit constants, not pure duplex thermodynamics, by the original model's own design choice. |

As of this writing, the only large, high-quality calibration panel that
exists is E. coli's (the Salis-lab 1014IC/8848FS-class datasets already
cited in `docs/MODEL.md`). Effectively every other organism in the world
has **zero** real Tier-3 ground truth. This is not a data-cleaning problem
or a modeling problem in the usual sense — it's a **N≈1 problem**, and it
is the single hardest and most important thing this design document has to
be honest about. §5 and §6 describe the approach; §11 opens the road to
actually growing N past 1.

---

## 3. Explicit Hypotheses

**H-P1 — Taxonomy-as-lookup-table is the composition-shortcut risk's
analogue here.** A network given a raw taxonomy string as an input feature
for Tier 2/3 prediction will, if allowed, use it as a lookup key ("if
Firmicutes, emit B. subtilis's numbers") rather than learning genuine
embedding-to-physiology structure — structurally the same failure the
encoder doc's H1 names for GC%/genus, one level downstream. *Falsification:*
ablate the raw taxonomy feature and check whether held-out Tier 2 accuracy
degrades marginally (embedding is doing real work) or catastrophically
(the network had memorized a lookup table). Unlike H1, this one has an easy
mitigation if confirmed — drop the feature — so it's a check, not a
structural risk requiring its own audit mechanism.

**H-P2 — The encoder's embedding space is close to linearly separable
along the physiological axes RBSForge needs.** If `genome_encoder`'s own
H6/§4.6 goals succeed (temperature/sigma-factor regime emerging as
structure without supervision), a *linear* probe on `E_genome` should
capture most of Tier 1/2's achievable accuracy, and a small MLP shouldn't
meaningfully beat it. *Falsification:* compare a linear head against a
small-MLP head on identical held-out data. *Why this matters beyond
architecture choice:* RBSForge's whole reason for existing is transparent,
term-by-term scoring, never a black box. A HostPack Decoder that needs deep
nonlinear machinery to work is in tension with the project it feeds — so
this hypothesis isn't just an efficiency preference, it's a values
commitment, and a negative result (linear isn't enough) should prompt
real reconsideration of scope, not just a bigger MLP.

**H-P3 — For `s_opt`, the classical native-RBS-survey estimate alone
already captures most of the achievable accuracy, and the learned
correction's marginal contribution is small.** *Falsification:* compare
survey-only vs. survey-plus-correction against every organism with any
real characterization (even partial/uncalibrated ones, e.g. the
Geobacillus "50K RBS" case already discussed in `docs/MODEL.md`). If the
delta is small, Head B should stay mostly classical rather than absorbing
more modeling investment — the same "cheapest mechanism that clears
validation" discipline the encoder doc applies throughout.

**H-P4 — Tier 3 prediction error for a held-out organism is monotonic in
its embedding distance to the nearest calibrated anchor.** This is the
load-bearing assumption underneath the entire Tier-3 approach in §5, and
it is currently *untested*, not just unconfirmed, because there's only one
real anchor (E. coli) to measure distance from. *Falsification (partial,
until N grows):* while real Tier-3 held-out data doesn't exist, check
whether embedding distance correlates with an independent physiological-
similarity proxy — 16S percent identity, shared Gram-stain/S1-architecture,
shared optimal-spacing regime from Tier 2 — as a weaker but currently
available substitute. A real test of H-P4 is one of the strongest arguments
for §11's active-learning proposal: it can't be properly checked without
at least one more real calibration point.

---

## 4. Component Specification

**Inputs consumed from `genome_encoder`** (per the encoder doc's own
component list):

- `E_genome` — the genome-summary embedding (post-loop, post-audit).
- `R_genome` — MoE expert-routing histogram (encoder doc §4.6). Brought in
  explicitly, not incidentally — see §10 for why this is a deliberate,
  slightly uncomfortable choice.
- `Δ(genome_embedding)` — the learned per-genome codon-table correction
  (encoder doc §4.7). A plausible weak feature for `start_codon_dG`; the
  same apparatus that corrects for genuine codon-table deviation may carry
  some signal about non-canonical start-codon usage, though this is
  speculative and untested (not asserted as a strong signal, just a cheap
  input to include and probe).
- Per-locus embedding for the identified 16S rRNA gene and its 3′ flank —
  needed for the Tier-1 boundary-refinement head.
- Transcription-unit-scale auxiliary output (operon boundaries, promoter/
  terminator positions) where validated — used by the native-RBS survey
  (§2) and, per §10, for TSS-anchored 5′ UTR windows.
- External metadata, never required, used when present, same "given not
  inferred" pattern the encoder doc applies to topology (§4.8 there): 16S
  copy number, genome size, coarse taxonomic rank (phylum/class only — see
  H-P1 on why organism-level taxonomy is deliberately excluded).

**Heads, one cluster per tier:**

- **Head A (Tier 1).** 16S 3′-boundary regression (a small window of
  candidate offsets around the annotated end) + Gram-stain classifier + S1/
  platform-architecture classifier. Trained on abundant, cheap, mostly-
  real-not-approximate labels (§6).
- **Head B (Tier 2).** OGT regression; the `s_opt` prior-anchored
  correction described in §2. Both are low-capacity (linear-or-shallow,
  per H-P2) regressors over `E_genome` (+ `R_genome`, + survey summary
  stats for `s_opt`).
- **Head C (Tier 3).** Not a standard regressor. Given N≈1, a deep
  multi-task regression head would simply memorize E. coli. Instead: a
  kernel/distance-weighted estimator over the (currently tiny) set of
  calibrated anchor organisms in embedding space — output is a posterior
  (mean, variance) per field, collapsing toward the nearest anchor's
  values, with variance that *only* shrinks near genuinely-calibrated
  embedding-space regions. This is deliberately closer to k-nearest-
  neighbors-with-uncertainty than to a learned deep function, because
  that's what N≈1–5 honestly supports. See §5 stage γ and §9 entry 3 for
  the gate on when this graduates to something more expressive.

**Explicit abstention.** Every head can output "no usable prediction" for
a given field, not just a wide-uncertainty number — mirroring RBSForge's
own existing "leaderless, skipped" pattern for start codons with no SD
match (`rbs/calculator.py`). For Head C specifically, this is the expected
common case early on: most of embedding space is not near E. coli, and a
number in that region is not "uncertain," it is not there yet.

---

## 5. Training Objectives (Staged)

Mirrors the encoder doc's staged-composite approach (§4.9 there) — bring
objectives online in order of how much they depend on infrastructure or
data that doesn't exist yet, rather than assuming everything is trainable
from day one.

**Stage α — Head A + Head B, supervised, no Tier-3 objective at all.**
Standard classification/regression losses against the label sources in §6
(rRNA maturation data, taxonomy-derived Gram stain, BacDive OGT metadata,
the native-RBS survey's own self-generated statistics). This stage is
cheap, needs no wet-lab RBS data whatsoever, and doubles as an *external*
validation of `genome_encoder`'s own H6/§4.6 claims — if Stage α can't get
reasonable OGT/Gram-stain accuracy from `E_genome`, that is itself evidence
worth feeding back to the encoder project, not just a failure of this one.

**Stage β — Self-consistency via native-RBS ranking.** A label-free
objective, computable for any genome with no wet-lab data at all: fold the
current predicted parameter set back through RBSForge's own thermodynamic
model and check whether it correctly ranks the genome's own high-confidence
native RBSs above shuffled/scrambled controls (and, where available, above
lower-expression genes by an independent proxy — codon adaptation index as
a weak stand-in, real Ribo-seq density where it exists, see §6). This is
the direct analogue of the encoder doc's own self-consistency stage (§4.9
item 3, protein-tower coherence as reward) — "coherence" here means "the
predicted HostPack makes this genome's own biology look sensible under
RBSForge's model," which is checkable at the full genome-panel scale
(thousands of genomes), not just the handful with Tier-3 ground truth.

**Stage γ — Few-shot calibration transfer for Head C.** Leave-one-
organism-out cross-validation is the only honest evaluation at this N.
Explicitly gated: Head C should not ship a point-estimate output for any
organism beyond the calibrated anchor set until at least some minimum
number of independent real calibration points exist — placeholder **N ≥ 5**
distinct organisms, chosen the same way the encoder doc treats its own
phase-transition gate (§4.7 there): a tunable number, set for real once a
second and third anchor exist to calibrate against, not asserted here as
final. Below that gate, Head C's only honest behavior is "same as nearest
calibrated relative, with an embedding-distance confidence penalty" — which
is a learned, embedding-distance-weighted version of exactly what
RBSForge's hand-curated HostPacks already do today (Geobacillus borrowing
B. subtilis's tail, explicitly marked as an unfit prior). Stage γ's job in
v0.1 is to make that existing, honest heuristic slightly less arbitrary,
not to leapfrog past it.

---

## 6. Omics and Data Sources

This is the practical bottleneck for the whole system, and the tiers above
map onto genuinely different data types with very different availability.
Ordered here by leverage — highest-value, most concretely actionable
first.

**1. Culture-collection phenotype metadata — the highest-leverage resource
for Tier 2, and probably the single best next step.** Databases such as
BacDive (DSMZ), and structured fields in NCBI BioSample/GOLD/IMG-M, hold
optimal growth temperature, pH range, salinity tolerance, oxygen
requirement, Gram stain, motility, and sporulation for tens of thousands of
strains, a meaningful fraction already genome-linked. This is large-N,
already-collected, and entirely independent of any RBS-specific
measurement — it is what makes Stage α's OGT regression realistic at all,
and there is real published precedent for genome-to-OGT prediction at this
scale reaching useful accuracy (proteome-composition and rRNA-GC-content-
based predictors, e.g. Sauer & Wang 2019-class approaches) that this
project should benchmark against rather than reinvent from scratch.

**2. rRNA 3′/5′-end mapping data — needed specifically for Tier 1's
anti-SD boundary head.** Term-seq (3′ end mapping) and dRNA-seq/Cappable-
seq (5′/TSS mapping) datasets, growing in number across both model and
non-model bacteria, plus curated structural-RNA resources (Rfam, SILVA,
RNAcentral) for comparative mature-end evidence across close relatives.
This is a narrow, well-defined data need, not a large one — the boundary-
refinement task only needs a local window of sequence and a true endpoint
per example.

**3. Ribosome profiling (Ribo-seq) — the highest-value data type for
breaking the Tier-3 bottleneck, and probably the headline recommendation
of this document.** Ribo-seq gives genome-wide, per-gene, in-vivo
translation-initiation signal, which is the closest thing to a native TIR
ground truth without a designed reporter library — and critically, public
Ribo-seq datasets already exist for a real and growing set of non-model
bacteria (thermophiles, Bacillus, Pseudomonas, cyanobacteria, among
others), reanalyzable today at zero new wet-lab cost. Two distinct uses:

- **At genome-panel scale:** feeds Stage β's self-consistency ranking with
  real signal instead of a codon-adaptation-index proxy, for every genome
  that has public Ribo-seq — likely dozens to low hundreds of organisms
  already, an order of magnitude beyond any hand-curated reporter panel.
- **At single-organism scale, where density is high enough:** a genome's
  own native genes, spanning a real range of SD strength and spacing
  (computed via RBSForge's own duplex search), become an implicit,
  free reporter library — in principle enough to *directly* fit an
  organism-specific β and spacing curve via regression against measured
  ribosome density, without designing or ordering a single construct. This
  is speculative and untested (§11 opens it as the concrete next step, not
  claimed as working here) but is the most plausible path to growing
  Tier-3's N past the handful reachable by commissioning new reporter
  libraries one at a time.

RNA-seq (bulk) is a useful, low-cost companion here, not a primary target
in its own right — ribosome-density needs mRNA abundance as a denominator
to become translation *efficiency* (closer to what TIR represents) rather
than a signal conflating transcription and translation, and the standard
Ribo-seq/RNA-seq ratio already does exactly this in the literature.

**4. Designed reporter libraries — the gold standard for Tier 3, and the
actual scarce resource.** The Salis-lab 1014IC/8848FS-class panels already
cited in `docs/MODEL.md` remain the best real calibration data that
exists, and are E. coli-only. New organism-specific libraries (what
`docs/MODEL.md` already calls "Layer 2" calibration, in language borrowed
directly from the source technical summary this whole project started
from) are the actual bottleneck this system exists to reduce *reliance on*,
not eliminate — see §11's active-learning framing for how to spend a small
number of new libraries where they matter most, rather than needing many.

**5. Proteomics and metabolomics/growth-condition metadata — secondary,
validation-tier.** Mass-spec protein abundance, combined with Ribo-seq, can
validate translation-efficiency estimates at the protein level; growth-
condition metadata matters because a single organism's OGT is a range, not
a point, and Ribo-seq/reporter data collected under a mismatched condition
is a real, easy-to-miss confound when comparing across studies. Worth
tracking, not worth blocking on.

**6. What NOT to use: MAGs, as already stated in §1.** Named again here
specifically because MAG-derived data is exactly the kind of thing that
looks like "more training data, for free" and is the most likely place
someone building this system later reaches for it without re-reading the
non-goal.

---

## 7. Compute and Scale

The inverse concern from `genome_encoder`. This network is a handful of
linear-or-shallow heads over a fixed, low-dimensional, already-computed
embedding — well under 1M parameters total, trivially inside the encoder
doc's own 4GB ceiling as a post-hoc addition, plausibly CPU-trainable given
both the model size and the realistic dataset sizes in §6. **This system
is data-constrained, specifically Tier-3-data-constrained, not compute-
constrained** — worth stating plainly so effort isn't misallocated toward
scaling a network that doesn't need it.

**Frozen-embedding, strictly downstream — not a nice-to-have, a
requirement.** `genome_encoder`'s protein and DNA towers are already
architecturally independent, coupled only through the back-feed interface
(encoder doc §4.2). This network should be looser still: a pure consumer
of `E_genome` (and the other frozen outputs listed in §4), with zero
gradient path back into the encoder. Joint training would risk destabilizing
an already carefully-staged training process (encoder doc §6) for a
downstream task the encoder was never designed around, and would make it
impossible to cleanly attribute a failure to "the embedding doesn't carry
this signal" versus "joint training broke something upstream." If `E_genome`
turns out not to carry enough signal (H-P2 failing), the right response is
feeding that back as a validation result for the encoder project to weigh,
not quietly wiring in a gradient to compensate.

---

## 8. Validation Strategy

- **Linear vs. shallow probing (H-P2).** Standard practice already, given
  the encoder doc's own validation culture.
- **Leave-one-organism/leave-one-genus-out cross-validation for Tier 3.**
  The only honest evaluation at N≈1–5. No held-out split of a single
  E. coli-anchored dataset can substitute for this.
- **Self-consistency ranking (§5 Stage β) as a broadly-applicable check.**
  Runs on any genome, calibrated or not — the closest thing this system has
  to a large-N validation axis.
- **Beat the naive baseline, explicitly.** The naive baseline is not a
  strawman — it's literally what RBSForge does today: borrow the nearest
  hand-picked phylogenetic relative's numbers, flagged `prior`. If Head C
  can't beat "just use B. subtilis's numbers for a Firmicute" on held-out
  organisms, it isn't earning its complexity, and Stage γ should say so
  rather than ship anyway.
- **Upstream leading indicator.** `genome_encoder`'s own linear probes for
  sigma-factor/RNAP-subunit presence and its MoE-specialization validation
  (encoder doc §6/§7) are a cheap early signal for this project. If those
  probes fail, that's evidence Tier 2 will also struggle, and Stage α
  work should be gated on the encoder project's own Stage C results rather
  than started blind.

---

## 9. Integration with RBSForge

Concrete, not aspirational — this is how the output actually lands in the
existing codebase.

- A HostPack Decoder prediction becomes a `HostPack` instance the same
  shape as any hand-curated one in `hostpack.py`, with one addition: a
  per-field confidence tier alongside the value, not just one
  whole-profile `calibration`/`uncertainty` pair as the current hand-
  curated entries have. Tier 1 fields can be genuinely `calibrated`-grade
  confident even for a brand-new organism; Tier 3 fields should almost
  never claim better than `prior`, and should default to `speculative` or
  outright abstain (§4) below Stage γ's gate.
- The RBSForge Console's existing badge system (`calibrated` / `prior` /
  `speculative`, `rbsforge-console.html`) is the right place to surface
  this once it exists — a decoder-generated HostPack should render with
  the same visual honesty already built for the hand-curated ones, not a
  separate "AI-generated" treatment that implies either more or less trust
  than the badges already communicate per-field.
- The native-RBS survey (§2) is new shared infrastructure, not decoder-
  specific — it's a direct reuse of `thermo/duplex.py` against
  genome-wide called ORFs, and is independently useful to RBSForge even
  without the rest of this network (e.g. as a "does this organism look
  SD-dependent at all" check before running Predict mode on it).

---

## 10. Dependencies on `genome_encoder`'s Open Questions

The encoder doc is explicitly WIP with several unresolved
`[QUESTION FOR NATHAN]` items. Most don't touch this design at all; two do,
directly enough to flag before either project commits further.

**H1a / routing-leakage resolution (encoder doc §4.10, open question 2).**
This is the sharpest dependency, and it's an uncomfortable one:
`R_genome` is useful to *this* network specifically because it may carry
composition/genus-correlated signal — which is precisely what the encoder
project's own audit (§4.10 there) is trying to scrub, and what H1a worries
is leaking there *unaudited*. If routing-leakage is confirmed and "fixed"
by folding routing into the live adversarial audit, the fix is explicitly
designed to destroy the correlation this network currently exploits.
**This is not an argument against fixing it upstream** — the encoder
project's goals are correct to prioritize over this one's convenience — but
it means Head B's use of `R_genome` should be re-validated (not assumed to
keep working) if and when that fix lands, and this dependency is worth
naming to whoever makes that call so it's a known tradeoff, not a surprise
regression discovered later.

**Replicon-provenance metadata (encoder doc §4.8, open question 5) and the
one-HostPack-per-genome simplification (§1 above).** If replicon-provenance
tagging is adopted upstream, it opens a real question this document
currently ducks: should a HostPack be per-genome or per-replicon? A
plasmid carrying its own heat-shock or alternative-sigma-factor machinery
is a real case where "the organism's" growth-temperature tolerance isn't
a single number. Not resolved here — deferred, consistent with H7's own
"local signal, real biology, not a shortcut to suppress" framing in the
encoder doc (§2 there) — but the two documents' open questions are the same
shape, and it's worth deciding them with awareness of each other rather
than independently.

---

## 11. Staged Decisions Register

1. **Minimum anchor count (N) gating Head C's point-estimate output.**
   Placeholder N≥5 (§5). Trigger to revisit: a second and third real
   calibration organism actually exist to set this empirically rather than
   by placeholder.
2. **Ribo-seq reanalysis as an in-scope Stage-β/Tier-3 data source.**
   Recommend adopting early — real data, zero new wet-lab cost, directly
   addresses the sharpest bottleneck in the whole design (§6). Low cost,
   high plausible value, in the same spirit as the encoder doc's own
   "adopt now rather than wait to discover the gap empirically" calls
   (e.g. its audit target-statistic widening, §4.10 there).
3. **Whether/which organism to commission a second real Tier-3 reporter
   library for.** Not decidable from this document alone — depends on
   what H-P4 (§3) and the Ribo-seq self-calibration idea (§6, item 3) show
   once tried. Candidates worth naming now: a Gram-positive with partial
   existing characterization, or a moderate thermophile, since either
   would break the current "everything is a distance from E. coli"
   degeneracy in a way a second mesophilic Gram-negative wouldn't.
4. **Per-replicon HostPacks.** Deferred per §10, contingent on the
   encoder doc's own replicon-provenance decision.
5. **`R_genome` re-validation trigger.** Tied to the encoder doc's H1a
   resolution — see §10. Not a decision to make now; a tripwire to make
   sure isn't missed later.

---

## 12. Open Questions

1. **(§5, §11 entry 1)** What's an acceptable placeholder for N in Head
   C's phase-transition gate before a second real anchor organism exists —
   is 5 too conservative, too permissive, or fine as a placeholder pending
   real data?
2. **(§6, §11 entry 2)** Is Ribo-seq reanalysis worth prioritizing now, or
   does it compete too directly with other near-term RBSForge/encoder work
   for it to be worth picking up yet?
3. **(§11 entry 3)** If a single new organism-specific reporter library
   were commissioned next, which organism would do the most to validate
   (or break) H-P4, given what's already partially characterized
   (Geobacillus, per `docs/MODEL.md`) versus genuinely novel?
4. **(§10)** Per-replicon HostPacks — worth scoping now alongside the
   encoder project's replicon-provenance question, or genuinely fine to
   defer until it's forced by a real case?
5. **(§4)** Is including `Δ(genome_embedding)` as a `start_codon_dG`
   feature worth the coupling to a mechanism (encoder doc §4.7) that is
   itself still gated near-zero in the encoder's own phase-one training —
   or should this wait until Δ is actually un-gated there, rather than
   depending on an input that's currently close to inert by design?
